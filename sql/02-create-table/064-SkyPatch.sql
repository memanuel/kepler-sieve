-- ************************************************************************************************
-- Set the variable N
SET @N = 4;

-- The six faces of a SkyPatch (correspond to six sides of an inscribed cube in the unit sphere)
CREATE OR REPLACE TABLE KS.SkyPatchFace(
    f TINYINT NOT NULL PRIMARY KEY
        COMMENT "Face counter from 1 to 6; formula is argmax(x, y, z, 1-x, 1-y, 1-z), e.g. 3 for the XY+ plane",
    SkyPatchFaceCD CHAR(3) NOT NULL UNIQUE,
    Axis_i CHAR(1) NOT NULL
    	COMMENT "The first of two axes that vary on this face, e.g. x",
    Axis_j CHAR(1) NOT NULL
    	COMMENT "The second of two axes that vary on this face, e.g. y",
    Axis_k CHAR(1) NOT NULL
    	COMMENT "The third axis; solved given the other two, e.g. z",
    Axis_k_sign DOUBLE NOT NULL
    	COMMENT "Sign of variables on this face; one of -1 and +1"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "The six faces of a SkyPatch, which correspond to six sides of an inscribed cube in the unit sphere.";

-- Manually populate the six major faces
INSERT INTO KS.SkyPatchFace
(f, SkyPatchFaceCD, Axis_i, Axis_j, Axis_k, Axis_k_sign)
VALUES
(1, 'YZ+', 'Y', 'Z', 'X', 1.0),
(2, 'XZ+', 'X', 'Z', 'Y', 1.0),
(3, 'XY+', 'X', 'Y', 'Z', 1.0),
(4, 'YZ-', 'Y', 'Z', 'X', -1.0),
(5, 'XZ-', 'X', 'Z', 'Y', -1.0),
(6, 'XY-', 'X', 'Y', 'Z', -1.0);

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
    -- Bound on u
    u0 DOUBLE NOT NULL COMMENT "Lower bound on u",
    u1 DOUBLE NOT NULL COMMENT "Upper bound on u",
    -- Bound on v
    v0 DOUBLE NOT NULL COMMENT "Lower bound on v",
    v1 DOUBLE NOT NULL COMMENT "Upper bound on v",
    -- Bound on w
    w0 DOUBLE NOT NULL COMMENT "Lower bound on w",
    w1 DOUBLE NOT NULL COMMENT "Upper bound on w",
    -- Unique key
    PRIMARY KEY (i,j)
    	COMMENT "The pair (i,j) uniquely determines one grid cell of a major face"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "SkyPatchGridCell describes the 4N^2 square grid cells on one major face";

-- Populate the 4N^2 grid cells on each major face
INSERT INTO KS.SkyPatchGridCell
(i, j, u, v, w, u0, u1, v0, v1, w0, w1)
WITH t1 AS (
SELECT
	-- Integer face coordinates (i,j)
	i._ AS i,
	j._ AS j,
	-- Midpoint of (u, v)
	CAST(( (i._+0.5)/N._ - 1.0) AS DOUBLE) AS u,
	CAST(( (j._+0.5)/N._ - 1.0) AS DOUBLE) AS v,
	-- Bound on u
	CAST((i._ /N._ - 1.0) AS DOUBLE) AS u0,
	CAST(( (i._+1)/N._ - 1.0) AS DOUBLE) AS u1,
	-- Bound on v
	CAST((j._ /N._ - 1.0) AS DOUBLE) AS v0,
	CAST(( (j._+1)/N._ - 1.0) AS DOUBLE) AS v1
FROM
	KS.Counter AS N
	INNER JOIN KS.Counter AS i ON i._ < 2*N._
	INNER JOIN KS.Counter AS j ON j._ < 2*N._
WHERE
	N._ = @N
), t2 AS (
SELECT
	-- Integer face coordinates (i,j)
	t1.i,
	t1.j,
	-- Midpoint (u,v)
	t1.u,
	t1.v,
	-- Bound on u
	t1.u0,
	t1.u1,
	-- Bound on v
	t1.v0,
	t1.v1,
	-- Squared Distance in (u,v) plane
	(POW(t1.u,2)+POW(t1.v,2)) AS rs,
	(POW(t1.u0,2)+POW(t1.v0,2)) AS rs0,
	(POW(t1.u1,2)+POW(t1.v1,2)) AS rs1
FROM
	t1
)
SELECT
	-- Integer face coordinates (i,j)
	t2.i,
	t2.j,
	-- Midpoint (u,v, w)
	t2.u,
	t2.v,
	SQRT(1.0-rs) AS w,
	-- Bound on u
	t2.u0,
	t2.u1,
	-- Bound on v
	t2.v0,
	t2.v1,
	-- Bound on w
	SQRT(1.0-rs0) AS w0,
	SQRT(1.0-rs1) AS w1
FROM
	t2
WHERE
	-- Only take points that fit inside the unit sphere!
	t2.rs0 <= 1.0 AND t2.rs1 <= 1.0;

-- Delete grid cells where w has no chance of being the larger than |u| and |v|
DELETE 
	KS.SkyPatchGridCell
FROM
	KS.SkyPatchGridCell
WHERE
	-- On these rows, |u| is always larger than w
	(pow(u0,2)>pow(w1,2) AND pow(u1,2)>pow(w1,2)) 
	OR
	-- On these rows, |v| is always larger than w
	(pow(v0,2)>pow(w1,2) AND pow(v1,2)>pow(w1,2));

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
