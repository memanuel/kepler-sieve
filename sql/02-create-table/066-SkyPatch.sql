-- ************************************************************************************************
-- Set the variable N
-- SET @N = POW(2,8);
SET @N = 4;
-- The size of each grid is 2Nx2N with 4N^2 entries
SET @M = 4 * @N * @N;

-- ************************************************************************************************
-- SkyPatchGridCell describes the 4N^2 square grid cells on one major face
CREATE OR REPLACE TABLE KS.SkyPatchGrid(
    i SMALLINT NOT NULL
        COMMENT "Counter for the first Cartesian coordinate on this face",
    j SMALLINT NOT NULL
        COMMENT "Counter for the second Cartesian coordinate on this face",
    k INT NOT NULL
    	COMMENT "Counter for the entire grid",
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
    -- Lower and upper bounds on u
    uMin double NOT NULL COMMENT "Lower bound on u",
    uMax double NOT NULL COMMENT "Upper bound on u",
    -- Lower and upper bounds on v
    vMin double NOT NULL COMMENT "Lower bound on v",
    vMax double NOT NULL COMMENT "Upper bound on v",
    -- Lower and upper bounds on w
    wMin double NOT NULL COMMENT "Lower bound on w",
    wMax double NOT NULL COMMENT "Upper bound on w",
    -- Unique key
    PRIMARY KEY (i,j)
    	COMMENT "The pair (i,j) uniquely determines one grid cell of a major face",
    UNIQUE KEY UNQ_SkyPatchGrid_k (k)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "SkyPatchGrid describes the 4N^2 square grid cells on one major face";

-- Populate the 4N^2 grid cells on each major face
INSERT INTO KS.SkyPatchGrid
(i, j, k, u, v, w, u00, v00, w00, u01, v01, w01, u10, v10, w10, u11, v11, w11,
 uMin, uMax, vMin, vMax, wMin, wMax)
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
	N._ = @N
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
	(2*@N*CAST(t2.i AS INT) + t2.j) AS k,
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
	t2.c/t2.r11 AS w11,
	-- Placeholder for the lower and upper bounds
	0.0 AS uMin,
	0.0 AS uMax,
	0.0 AS vMin,
	0.0 AS vMax,
	0.0 AS wMin,
	0.0 AS wMax	
FROM
	t2;

-- Set the bounds on u, v, w
UPDATE KS.SkyPatchGrid
SET
	uMin = LEAST(u00, u01, u10, u11),
	uMax = GREATEST(u00, u01, u10, u11),
	vMin = LEAST(v00, v01, v10, v11),
	vMax = GREATEST(v00, v01, v10, v11),
	wMin = LEAST(w00, w01, w10, w11),
	wMax = GREATEST(w00, w01, w10, w11);

-- ************************************************************************************************
-- The whole SkyPatch table; six major faces
CREATE OR REPLACE TABLE KS.SkyPatch(
	SkyPatchID INT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this patch of sky; assigned by formula (2N)^2*(f-1) + (2N)i + j",
    CubeFaceID TINYINT NOT NULL
        COMMENT "The major face of the cube on which the 2Nx2N grid is inscribed; foreign key to CubeFace",
    i SMALLINT NOT NULL
        COMMENT "Counter for the first Cartesian coordinate on this CubeFace, which is named alpha on the CubeFace table",
    j SMALLINT NOT NULL
        COMMENT "Counter for the second Cartesian coordinate on this CubeFace, which is named beta on the CubeFace table",
    -- Center point
    x DOUBLE NOT NULL COMMENT "Midpoint of x",
    y DOUBLE NOT NULL COMMENT "Midpoint of y",
    z DOUBLE NOT NULL COMMENT "Midpoint of z",
    -- Bound on x
    x0 DOUBLE NOT NULL COMMENT "Lower bound on x",
    x1 DOUBLE NOT NULL COMMENT "Upper bound on x",
    -- Bound on y
    y0 DOUBLE NOT NULL COMMENT "Lower bound on y",
    y1 DOUBLE NOT NULL COMMENT "Upper bound on y",
    -- Bound on z
    z0 DOUBLE NOT NULL COMMENT "Lower bound on z",
    z1 DOUBLE NOT NULL COMMENT "Upper bound on z",
    -- Unique key
    UNIQUE KEY UNQ_SkyPatch_CubeFaceID_i_j(CubeFaceID, i, j)
    	COMMENT "The pair (i,j) determines one small patch on a major face; the trio (f,i,j) is unique"
    -- Foreign key
    -- CONSTRAINT FK_SkyPatch_CubeFaceID FOREIGN KEY (CubeFaceID) REFERENCES CubeFace(CubeFaceID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Collection of discrete patches of the sky corresponding to a cube in which the unit sphere is inscribed.";

-- SELECT * FROM KS.SkyPatch;
INSERT INTO KS.SkyPatch
(SkyPatchID, CubeFaceID, i, j, x, y, z, x0, x1, y0, y1, z0, z1)
SELECT
	-- Integer IDs
	(cf.CubeFaceID-1)*@M + 2*gr.i*@N + j AS SkyPatchID,
	cf.CubeFaceID,
	-- Local grid coordinates (i, j) on the major face
	gr.i,
	gr.j,
	-- Coordinates of midpoint (x, y, z)
	IF(cf.alpha='X', gr.u, 0.0) + IF(cf.beta='X', gr.v, 0.0) + IF(cf.gamma='X', gr.w, 0.0) AS x,
	IF(cf.alpha='Y', gr.u, 0.0) + IF(cf.beta='Y', gr.v, 0.0) + IF(cf.gamma='Y', gr.w, 0.0) AS y,
	IF(cf.alpha='Z', gr.u, 0.0) + IF(cf.beta='Z', gr.v, 0.0) + IF(cf.gamma='Z', gr.w, 0.0) AS z,
	-- Bounds on x
	IF(cf.alpha='X', gr.uMin, 0.0) + IF(cf.beta='X', gr.vMin, 0.0) + IF(cf.gamma='X', gr.wMin, 0.0) AS x0,	
	IF(cf.alpha='X', gr.uMax, 0.0) + IF(cf.beta='X', gr.vMax, 0.0) + IF(cf.gamma='X', gr.wMax, 0.0) AS x1,
	-- Bounds on y
	IF(cf.alpha='Y', gr.uMin, 0.0) + IF(cf.beta='Y', gr.vMin, 0.0) + IF(cf.gamma='Y', gr.wMin, 0.0) AS y0,	
	IF(cf.alpha='Y', gr.uMax, 0.0) + IF(cf.beta='Y', gr.vMax, 0.0) + IF(cf.gamma='Y', gr.wMax, 0.0) AS y1,
	-- Bounds on z
	IF(cf.alpha='Z', gr.uMin, 0.0) + IF(cf.beta='Z', gr.vMin, 0.0) + IF(cf.gamma='Z', gr.wMin, 0.0) AS z0,	
	IF(cf.alpha='Z', gr.uMax, 0.0) + IF(cf.beta='Z', gr.vMax, 0.0) + IF(cf.gamma='Z', gr.wMax, 0.0) AS z1
FROM
	KS.CubeFace AS cf,
	KS.SkyPatchGrid AS gr
ORDER BY SkyPatchID;

-- ************************************************************************************************
-- Neighbor distance on the SkyPatchGrid table
CREATE OR REPLACE TABLE KS.SkyPatchGridNeighbor(
	i1 INT NOT NULL
        COMMENT "i coordinate of the first SkyPatchGrid cell",
	j1 INT NOT NULL
        COMMENT "j coordinate of the first SkyPatchGrid cell",
	i2 INT NOT NULL
        COMMENT "i coordinate of the second SkyPatchGrid cell",
	j2 INT NOT NULL
        COMMENT "ji coordinate of the second SkyPatchGrid cell",
	dr DOUBLE NOT NULL
		COMMENT "Distance from center point to center point",
	dr_min DOUBLE NOT NULL
		COMMENT "Estimated minimum distance based on getting closer by up to half the width on each end",
	PRIMARY KEY (i1, j1, i2, j2)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Catalog neighbor interactions on SkyPatchGrid table.";

-- Width for populating the SkyPatchGridNeighbor table
SET @dw = 5;

INSERT INTO KS.SkyPatchGridNeighbor
(i1, j1, i2, j2, dr, dr_min)
WITH t1 as(
SELECT
	g1.i AS i1,
	g1.j AS j1,
	g2.i AS i2,
	g2.j AS j2,
	-- The distance along each coordinate from center to center
 	ABS(g2.u-g1.u) AS du,
 	ABS(g2.v-g1.v) AS dv,
 	ABS(g2.w-g1.w) AS dw,
 	-- Estimated distance improvement by traveling from the center to the nearest boundary
 	((g1.uMax-g1.uMin)+(g2.uMax-g2.uMin))/2 AS wu,
 	((g1.vMax-g1.vMin)+(g2.vMax-g2.vMin))/2 AS wv,
 	((g1.wMax-g1.wMin)+(g2.wMax-g2.wMin))/2 AS ww
FROM
	KS.SkyPatchGrid AS g1
	INNER JOIN KS.SkyPatchGrid AS g2 ON 
		(g2.i BETWEEN g1.i-@dw AND g1.i+@dw) AND
		(g2.j BETWEEN g1.j-@dw AND g1.j+@dw)
), t2 AS(
SELECT
	t1.i1,
	t1.j1,
	t1.i2,
	t1.j2,
	t1.du,
	t1.dv,
	t1.dw,
	GREATEST(t1.du - t1.wu, 0.0) AS du_min,
	GREATEST(t1.dv - t1.wv, 0.0) AS dv_min,
	GREATEST(t1.dw - t1.ww, 0.0) AS dw_min
FROM
	t1
)
SELECT
	t2.i1,
	t2.j1,
	t2.i2,
	t2.j2,
	SQRT(POW(du,2)+POW(dv,2) + POW(dw,2)) AS dr,
	SQRT(POW(du_min,2)+POW(dv_min,2) + POW(dw_min,2)) AS dr_min
FROM
	t2;

-- ************************************************************************************************
-- Neighbor distance on the SkyPatch table
CREATE OR REPLACE TABLE KS.SkyPatchNeighbor(
	SkyPatchID_1 INT NOT NULL
        COMMENT "The first SkyPatchID",
	SkyPatchID_2 INT NOT NULL
        COMMENT "The second SkyPatchID",
	dr DOUBLE NOT NULL
		COMMENT "Distance from center point to center point",
	dr_min DOUBLE NOT NULL
		COMMENT "Estimated minimum distance based on getting closer by up to half the width on each end",
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2)
)
COMMENT "Catalog neighbor interactions on SkyPatch table.";
