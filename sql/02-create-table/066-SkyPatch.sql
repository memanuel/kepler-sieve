-- ************************************************************************************************
-- Set the variable N
SET @N = 4;

-- The six faces of a SkyPatch (correspond to six sides of an inscribed cube in the unit sphere)
CREATE OR REPLACE TABLE KS.SkyPatchFace(
    f TINYINT NOT NULL PRIMARY KEY
        COMMENT "Face counter from 1 to 6; formula is argmax(z, y, x, -z, -y, -x), e.g. 1 for the XY+ plane",
    SkyPatchFaceCD CHAR(3) NOT NULL UNIQUE,
    alpha CHAR(1) NOT NULL
    	COMMENT "The first of two axes that vary on this face, e.g. x; discrete index i",
    beta CHAR(1) NOT NULL
    	COMMENT "The second of two axes that vary on this face, e.g. y; discrete index j",
    gamma CHAR(1) NOT NULL
    	COMMENT "The third axis; solved given the other two, e.g. z",
    c DOUBLE NOT NULL
    	COMMENT "Value of gamma axis on this face; one of -1 and +1"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "The six faces of a SkyPatch, which correspond to six sides of an inscribed cube in the unit sphere.";

-- Manually populate the six major faces
INSERT INTO KS.SkyPatchFace
(f, SkyPatchFaceCD, alpha, beta, gamma, c)
VALUES
(1, 'XY+', 'X', 'Y', 'Z', 1.0),
(2, 'XZ+', 'X', 'Z', 'Y', 1.0),
(3, 'YZ+', 'Y', 'Z', 'X', 1.0),
(4, 'XY-', 'X', 'Y', 'Z', -1.0),
(5, 'XZ-', 'X', 'Z', 'Y', -1.0),
(6, 'YZ-', 'Y', 'Z', 'X', -1.0);

-- ************************************************************************************************
-- SkyPatchGridCell describes the 4N^2 square grid cells on one major face
CREATE OR REPLACE TABLE KS.SkyPatchGridCell(
    i SMALLINT UNSIGNED NOT NULL
        COMMENT "Counter for the first Cartesian coordinate on this face",
    j SMALLINT UNSIGNED NOT NULL
        COMMENT "Counter for the second Cartesian coordinate on this face",
    -- Center point
    u DOUBLE NOT NULL COMMENT "Midpoint of u (coordinate corresponding to i)",
    v DOUBLE NOT NULL COMMENT "Midpoint of v (coordinate corresponding to j)",
    w DOUBLE NOT NULL COMMENT "Midpoint of w (third coordinate, not i or j)",
	-- Lower left corner (u00, v00, w00)
    u00 DOUBLE NOT NULL COMMENT "Lower left corner, u",
    v00 DOUBLE NOT NULL COMMENT "Lower left corner, v",
    w00 DOUBLE NOT NULL COMMENT "Lower left corner, w",
	-- Upper left corner (u01, v01, w01)
    u01 DOUBLE NOT NULL COMMENT "Upper left corner, u",
    v01 DOUBLE NOT NULL COMMENT "Upper left corner, v",
    w01 DOUBLE NOT NULL COMMENT "Upper left corner, w",
	-- Lower right corner (u10, v10, w10)
    u10 DOUBLE NOT NULL COMMENT "Lower right corner, u",
    v10 DOUBLE NOT NULL COMMENT "Lower right corner, v",
    w10 DOUBLE NOT NULL COMMENT "Lower right corner, w",
	-- Upper right corner (u11, v11, w11)
    u11 DOUBLE NOT NULL COMMENT "Upper right corner, u",
    v11 DOUBLE NOT NULL COMMENT "Upper right corner, v",
    w11 DOUBLE NOT NULL COMMENT "Upper right corner, w",       
    -- Unique key
    PRIMARY KEY (i,j)
    	COMMENT "The pair (i,j) uniquely determines one grid cell of a major face"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "SkyPatchGridCell describes the 4N^2 square grid cells on one major face";

-- Populate the 4N^2 grid cells on each major face
INSERT INTO KS.SkyPatchGridCell
(i, j, u, v, w, u00, v00, w00, u01, v01, w01, u10, v10, w10, u11, v11, w11)
WITH t1 AS (
SELECT
	-- Integer face coordinates (i,j)
	i._ AS i,
	j._ AS j,
	-- Midpoint of (a, b)
	CAST(( (i._+0.5)/N._ - 1.0) AS DOUBLE) AS a,
	CAST(( (j._+0.5)/N._ - 1.0) AS DOUBLE) AS b,
	CAST(1.0 AS DOUBLE) AS c,
	-- Bound on a
	CAST((i._ /N._ - 1.0) AS DOUBLE) AS a0,
	CAST(( (i._+1)/N._ - 1.0) AS DOUBLE) AS a1,
	-- Bound on b
	CAST((j._ /N._ - 1.0) AS DOUBLE) AS b0,
	CAST(( (j._+1)/N._ - 1.0) AS DOUBLE) AS b1
FROM
	KS.Counter AS N
	INNER JOIN KS.Counter AS i ON i._ < 2*N._
	INNER JOIN KS.Counter AS j ON j._ < 2*N._
WHERE
	N._ = 4
), t2 AS(
SELECT
	-- Integer face coordinates (i,j)
	t1.i,
	t1.j,
	-- Midpoint (a, b, c)
	t1.a,
	t1.b,
	t1.c,
	-- Radius r of point (a, b, c) on cube
	SQRT(POW(t1.a, 2) + POW(t1.b, 2) + POW(t1.c, 2)) AS r,
	-- Bound on a
	t1.a0,
	t1.a1,
	-- Bound on b
	t1.b0,
	t1.b1,
	-- Radius at 4 corners
	SQRT(POW(t1.a0, 2) + POW(t1.b0, 2) + POW(t1.c, 2)) AS r00,
	SQRT(POW(t1.a0, 2) + POW(t1.b1, 2) + POW(t1.c, 2)) AS r01,
	SQRT(POW(t1.a1, 2) + POW(t1.b0, 2) + POW(t1.c, 2)) AS r10,
	SQRT(POW(t1.a1, 2) + POW(t1.b1, 2) + POW(t1.c, 2)) AS r11
FROM
	t1
)
SELECT
	-- Integer face coordinates (i,j)
	t2.i,
	t2.j,
	-- Midpoint (u, v, w)
	t2.a/r AS u,
	t2.b/r AS v,
	t2.c/r AS w,
	-- Lower left corner (u00, v00, w00)
	t2.a0/t2.r00 AS u00, 
	t2.b0/t2.r00 AS v00, 
	t2.c/t2.r00 AS w00,
	-- Upper left corner (u01, v01, w01)
	t2.a0/t2.r01 AS u01, 
	t2.b1/t2.r01 AS v01, 
	t2.c/t2.r01 AS w01,
	-- Lower right corner (u10, v10, w10)
	t2.a1/t2.r10 AS u10, 
	t2.b0/t2.r10 AS v10, 
	t2.c/t2.r10 AS w10,
	-- Upper right corner (u11, v11, w11)
	t2.a1/t2.r11 AS u11,
	t2.b1/t2.r11 AS v11,
	t2.c/t2.r11 AS w11
FROM
	t2;

-- ************************************************************************************************
-- The whole SkyPatch table; six major faces
CREATE OR REPLACE TABLE KS.SkyPatch(
	SkyPatchID INT UNSIGNED NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this one patch of sky",
    f TINYINT NOT NULL
        COMMENT "Face counter from 1 to 6",
    i SMALLINT UNSIGNED NOT NULL
        COMMENT "Counter for the first Cartesian coordinate on this face",
    j SMALLINT UNSIGNED NOT NULL
        COMMENT "Counter for the second Cartesian coordinate on this face",
    -- Bound on x
    x0 DOUBLE NOT NULL COMMENT "Lower bound on x",
    x1 DOUBLE NOT NULL COMMENT "Upper bound on x",
    -- Bound on y
    y0 DOUBLE NOT NULL COMMENT "Lower bound on y",
    y1 DOUBLE NOT NULL COMMENT "Upper bound on y",
    -- Bound on z
    z0 DOUBLE NOT NULL COMMENT "Lower bound on z",
    z1 DOUBLE NOT NULL COMMENT "Upper bound on z",
    -- Center point
    x DOUBLE NOT NULL COMMENT "Midpoint of x",
    y DOUBLE NOT NULL COMMENT "Midpoint of y",
    z DOUBLE NOT NULL COMMENT "Midpoint of z",
    -- Unique key
    UNIQUE KEY UNQ_SkyPatch_f_i_j(f, i, j)
    	COMMENT "The pair (i,j) determines one small patch on a major face; the trio (f,i,j) is unique",
    -- Foreign key
    CONSTRAINT FK_SkyPatch_f FOREIGN KEY (f) REFERENCES SkyPatchFace(f)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Time at which one or more detections were made an observatory.  Cartesian position and velocity includes topos adjustment.";
