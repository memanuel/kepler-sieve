DELIMITER $$

/* *********************************************************************************
Procedures in this module:
KS.MakeTable_SkyPatchGrid(N INT, dr_max DOUBLE)
KS.MakeTable_SkyPatch()
KS.MakeTable_SkyPatchDistance_face(dr_max DOUBLE)
KS.MakeTable_SkyPatchDistance_extend(dr_max DOUBLE)
KS.MakeTable_SkyPatchDistance_min()
* **********************************************************************************/

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchGrid(
	IN N INT,
	IN dr_max DOUBLE)
COMMENT "Populate the KS.SkyPatchGrid and KS.SkyPatchGridDistance tables"
BEGIN 

-- Calculate M and N as a float
SET @M = 2*N;
SET @Nf = CAST(N AS DOUBLE);

-- Calculate grid_width based on dr_max
SET @gw = CAST(1 + CEILING(dr_max*N) AS INT);

-- Calculate maximum sum of squares of di^2 + dj^2 based on dr_max
-- SET @gw2 = CAST(CEILING(1+POW(dr_max*N,2)) AS INT);

-- Empty both tables
TRUNCATE TABLE KS.SkyPatchGrid;
TRUNCATE TABLE KS.SkyPatchGridDistance;

-- Populate the M^2 grid cells in SkyPatchGrid
INSERT INTO KS.SkyPatchGrid
(i, j, k, u, v, w, u00, v00, w00, u01, v01, w01, u10, v10, w10, u11, v11, w11)
WITH t1 AS (
SELECT
	-- Integer face coordinates (i,j)
	i._ AS i,
	j._ AS j,
	-- Midpoint of (a, b)
	CAST(( (i._+0.5)/@Nf - 1.0) AS DOUBLE) AS a,
	CAST(( (j._+0.5)/@Nf - 1.0) AS DOUBLE) AS b,
	CAST(1.0 AS DOUBLE) AS c,
	-- Bound on a
	CAST((i._ /@Nf - 1.0) AS DOUBLE) AS a0,
	CAST(( (i._+1)/@Nf - 1.0) AS DOUBLE) AS a1,
	-- Bound on b
	CAST((j._ /@Nf - 1.0) AS DOUBLE) AS b0,
	CAST(( (j._+1)/@Nf - 1.0) AS DOUBLE) AS b1
FROM
	KS.Counter AS i CROSS JOIN
	KS.Counter AS j 
WHERE	
	(i._ < @M) AND (j._ < @M)
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
	(2*N*CAST(t2.i AS INT) + t2.j) AS k,
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

-- Populate SkyPatchGridDistance; just pairs of points less than dr_max apart
-- Step 1: get the coordinates of each candidate pair of corners
CREATE OR REPLACE TEMPORARY TABLE t1(
	i1 INT NOT NULL,
	j1 INT NOT NULL,
	i2 INT NOT NULL,
	j2 INT NOT NULL,
	cr1 TINYINT NOT NULL,
	cr2 TINYINT NOT NULL,
	u1 DOUBLE NOT NULL,
	v1 DOUBLE NOT NULL,
	w1 DOUBLE NOT NULL,
	u2 DOUBLE NOT NULL,
	v2 DOUBLE NOT NULL,
	w2 DOUBLE NOT NULL,
	dr_mid DOUBLE NOT NULL,
	dr_corner DOUBLE NOT NULL,
	-- This temp table keyed by the grid points and choice of corners
	PRIMARY KEY (i1, j1, i2, j2, cr1, cr2)
);

INSERT INTO t1
(i1, j1, i2, j2, cr1, cr2, u1, v1, w1, u2, v2, w2, dr_mid, dr_corner)
SELECT
	g1.i AS i1,
	g1.j AS j1,
	g2.i AS i2,
	g2.j AS j2,
 	-- Selected corners
 	cr1._ AS cr1,
 	cr2._ AS cr2,
 	-- Selected corner from grid cell 1
 	CASE cr1._ WHEN 0 THEN g1.u00 WHEN 1 THEN g1.u01 WHEN 2 THEN g1.u10	WHEN 3 THEN g1.u11 END AS u1,
 	CASE cr1._ WHEN 0 THEN g1.v00 WHEN 1 THEN g1.v01 WHEN 2 THEN g1.v10	WHEN 3 THEN g1.v11 END AS v1,
 	CASE cr1._ WHEN 0 THEN g1.w00 WHEN 1 THEN g1.w01 WHEN 2 THEN g1.w10	WHEN 3 THEN g1.w11 END AS w1,
 	-- Selected corner from grid cell 2
 	CASE cr2._ WHEN 0 THEN g2.u00 WHEN 1 THEN g2.u01 WHEN 2 THEN g2.u10	WHEN 3 THEN g2.u11 END AS u2,
 	CASE cr2._ WHEN 0 THEN g2.v00 WHEN 1 THEN g2.v01 WHEN 2 THEN g2.v10	WHEN 3 THEN g2.v11 END AS v2,
 	CASE cr2._ WHEN 0 THEN g2.w00 WHEN 1 THEN g2.w01 WHEN 2 THEN g2.w10	WHEN 3 THEN g2.w11 END AS w2,
	-- Distance at midpoint
	SQRT(POW(g2.u-g1.u,2) + POW(g2.v-g1.v,2) + POW(g2.w-g1.w,2)) AS dr_mid,
	0.0 AS dr_corner
FROM
	-- The starting grid cell
	KS.SkyPatchGrid AS g1
	-- The change in i and j
	INNER JOIN KS.CounterSigned AS di ON di._ BETWEEN - @gw AND @gw
	INNER JOIN KS.CounterSigned AS dj ON dj._ BETWEEN - @gw AND @gw
	INNER JOIN KS.SkyPatchGrid AS g2 ON 
		g2.i = g1.i + di._ AND
		g2.j = g1.j + dj._
	-- Counter to choose the corner of grid cell 1
	INNER JOIN KS.Counter AS cr1 ON cr1._ < 2
	-- Counter to choose the corner of grid cell 2
	INNER JOIN KS.Counter AS cr2 ON cr2._ < 2;

-- Calculate distance on the corners
UPDATE t1 SET dr_corner = SQRT( POW(u2-u1,2) + POW(v2-v1,2) + POW(w2-w1,2));

-- Step 2: group by the pair of candidate grid points
CREATE OR REPLACE TEMPORARY TABLE t2(
	i1 INT NOT NULL,
	j1 INT NOT NULL,
	i2 INT NOT NULL,
	j2 INT NOT NULL,
	dr_mid DOUBLE NOT NULL,
	dr_min DOUBLE NOT NULL,
	PRIMARY KEY (i1, j1, i2, j2)
);

INSERT INTO t2
(i1, j1, i2, j2, dr_mid, dr_min)
SELECT
	t1.i1,
	t1.j1,
	t1.i2,
	t1.j2,
	MIN(t1.dr_mid) AS dr_mid,
	MIN(t1.dr_corner) AS dr_min
FROM
	t1
-- Group by the pair of points so we take only the minimum distance between grid points
GROUP BY t1.i1, t1.j1, t1.i2, t1.j2;

-- Now insert all the pairs of grid points that are less than dr_max apart their closest corners
INSERT INTO KS.SkyPatchGridDistance
(i1, j1, i2, j2, dr_mid, dr_min)
SELECT
	t2.i1,
	t2.j1,
	t2.i2,
	t2.j2,
	t2.dr_mid,
	t2.dr_min
FROM
	t2
WHERE
	t2.dr_min < dr_max;

-- Clean up temporary tables
-- DROP TEMPORARY TABLE t1;
-- DROP TEMPORARY TABLE t2;

END $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatch()
COMMENT "Populate the KS.SkyPatch table"
BEGIN 

-- Calculate @M and @N
SELECT 
	MAX(spg.i)+1 AS M 
INTO @M 
FROM KS.SkyPatchGrid AS spg;
SET @N = (@M DIV 2);	

TRUNCATE TABLE KS.SkyPatch;

-- Build SkyPatch by joining CubeFace to SkyPatchGrid
INSERT INTO KS.SkyPatch
(SkyPatchID, CubeFaceID, i, j, x, y, z, x00, y00, z00, x01, y01, z01, x10, y10, z10, x11, y11, z11)
SELECT
	-- Integer IDs
	(cf.CubeFaceID-1)*@M*@M + gr.i*@M + j AS SkyPatchID,
	cf.CubeFaceID,
	-- Local grid coordinates (i, j) on the major face
	gr.i,
	gr.j,
	-- Coordinates of midpoint (x, y, z)
	IF(cf.alpha='X', gr.u, 0.0) + IF(cf.beta='X', gr.v, 0.0) + IF(cf.gamma='X', gr.w*cf.ci, 0.0) AS x,
	IF(cf.alpha='Y', gr.u, 0.0) + IF(cf.beta='Y', gr.v, 0.0) + IF(cf.gamma='Y', gr.w*cf.ci, 0.0) AS y,
	IF(cf.alpha='Z', gr.u, 0.0) + IF(cf.beta='Z', gr.v, 0.0) + IF(cf.gamma='Z', gr.w*cf.ci, 0.0) AS z,
	-- Coordinates of lower left corner (x00, y00, z00)
	IF(cf.alpha='X', gr.u00, 0.0) + IF(cf.beta='X', gr.v00, 0.0) + IF(cf.gamma='X', gr.w00*cf.ci, 0.0) AS x00,
	IF(cf.alpha='Y', gr.u00, 0.0) + IF(cf.beta='Y', gr.v00, 0.0) + IF(cf.gamma='Y', gr.w00*cf.ci, 0.0) AS y00,
	IF(cf.alpha='Z', gr.u00, 0.0) + IF(cf.beta='Z', gr.v00, 0.0) + IF(cf.gamma='Z', gr.w00*cf.ci, 0.0) AS z00,
	-- Coordinates of upper left corner (x01, y01, z01)
	IF(cf.alpha='X', gr.u01, 0.0) + IF(cf.beta='X', gr.v01, 0.0) + IF(cf.gamma='X', gr.w01*cf.ci, 0.0) AS x01,
	IF(cf.alpha='Y', gr.u01, 0.0) + IF(cf.beta='Y', gr.v01, 0.0) + IF(cf.gamma='Y', gr.w01*cf.ci, 0.0) AS y01,
	IF(cf.alpha='Z', gr.u01, 0.0) + IF(cf.beta='Z', gr.v01, 0.0) + IF(cf.gamma='Z', gr.w01*cf.ci, 0.0) AS z01,
	-- Coordinates of lower right corner (x10, y10, z10)
	IF(cf.alpha='X', gr.u10, 0.0) + IF(cf.beta='X', gr.v10, 0.0) + IF(cf.gamma='X', gr.w10*cf.ci, 0.0) AS x10,
	IF(cf.alpha='Y', gr.u10, 0.0) + IF(cf.beta='Y', gr.v10, 0.0) + IF(cf.gamma='Y', gr.w10*cf.ci, 0.0) AS y10,
	IF(cf.alpha='Z', gr.u10, 0.0) + IF(cf.beta='Z', gr.v10, 0.0) + IF(cf.gamma='Z', gr.w10*cf.ci, 0.0) AS z10,
	-- Coordinates of upper right corner (x11, y11, z11)
	IF(cf.alpha='X', gr.u11, 0.0) + IF(cf.beta='X', gr.v11, 0.0) + IF(cf.gamma='X', gr.w11*cf.ci, 0.0) AS x11,
	IF(cf.alpha='Y', gr.u11, 0.0) + IF(cf.beta='Y', gr.v11, 0.0) + IF(cf.gamma='Y', gr.w11*cf.ci, 0.0) AS y11,
	IF(cf.alpha='Z', gr.u11, 0.0) + IF(cf.beta='Z', gr.v11, 0.0) + IF(cf.gamma='Z', gr.w11*cf.ci, 0.0) AS z11
FROM
	KS.CubeFace AS cf CROSS JOIN
	KS.SkyPatchGrid AS gr
ORDER BY SkyPatchID;

END $$
	
-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchDistance_face(
	IN dr_max DOUBLE)
COMMENT "Populate the KS.SkyPatchDistance table with records that are on the same face"
BEGIN 

-- Empty table first so this procedure can be re-run on demand
TRUNCATE TABLE KS.SkyPatchDistance;
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min, IsCrossFace)
SELECT
	-- The two SkyPatch cells in this pair
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	-- The midpoint and minimum distance
	gd.dr_mid,
	gd.dr_min,
	FALSE AS IsCrossFace
FROM
	-- Start with the cube faces
	KS.CubeFace AS cf
	-- Each cube face is crossed with the GridDistance table; only those within the maximum distance
	INNER JOIN KS.SkyPatchGridDistance AS gd ON gd.dr_min < dr_max
	-- The two SkyPatch cells sharing this cube face and pair of grid coordinates
	INNER JOIN KS.SkyPatch AS p1 ON	p1.CubeFaceID = cf.CubeFaceID AND p1.i=gd.i1 AND p1.j=gd.j1
	INNER JOIN KS.SkyPatch AS p2 ON	p2.CubeFaceID = cf.CubeFaceID AND p2.i=gd.i2 AND p2.j=gd.j2;

END $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchDistance_extend_bf(
	IN dr_max DOUBLE)
COMMENT "Populate records on the KS.SkyPatchDistance table spanning across faces using brute force"
BEGIN 

-- Brute force approach for distance on pairs that span two different faces
-- This does not scale, but useful for checking with small N
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min, IsCrossFace)
WITH t1 AS (
SELECT
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	SQRT(POW(p2.x-p1.x,2)+POW(p2.y-p1.y,2)+POW(p2.z-p1.z,2)) AS dr_mid
FROM
	KS.SkyPatch AS p1
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID <> p1.CubeFaceID
)
SELECT
	t1.SkyPatchID_1,
	t1.SkyPatchID_2,
	t1.dr_mid,
	-1.0 AS dr_min,
	TRUE AS IsCrossFace
FROM
	t1
WHERE 
	t1.dr_mid < dr_max;

-- Calculate the minimum distance
CALL KS.MakeTable_SkyPatchDistance_min();
	
END $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchDistance_extend(
	IN dr_max DOUBLE)
COMMENT "Extend the KS.SkyPatchDistance table to include pairs that span different cube faces"
BEGIN 

-- Get M and calculate N
SELECT
	MAX(i)+1 AS M
INTO @M
FROM KS.SkyPatchGrid;

SET @N = @M DIV 2;
	
-- Temporary table for direct neighbors that span a face
CREATE OR REPLACE TABLE KS.SkyPatchNeighbors(
	SkyPatchID_1 INT NOT NULL,
	SkyPatchID_2 INT NOT NULL,
	dr DOUBLE NOT NULL,
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2)
) 
ENGINE='Aria' TRANSACTIONAL=0;

-- Insert the first batch of neighbors, on the i direction (i0 or i1)
INSERT INTO KS.SkyPatchNeighbors
(SkyPatchID_1, SkyPatchID_2, dr)
SELECT
	-- The two SkyPatch cells
 	sp1.SkyPatchID AS SkyPatchID_1,
 	sp2.SkyPatchID AS SkyPatchID_2,
 	-- Distance between these points
	MIN(SQRT(POW(sp2.x-sp1.x,2)+POW(sp2.y-sp1.y,2)+POW(sp2.z-sp1.z,2))) AS dr
FROM 
	-- Start with one cube face
	KS.CubeFace AS cf1
	-- Neighbors of this CubeFace
	INNER JOIN KS.CubeFaceNeighbor AS cfn ON cfn.CubeFaceID = cf1.CubeFaceID
	-- The value of i whose neighbors are being found; 0 for left (i0) and 1 for right (i1)
	INNER JOIN KS.Counter AS i1 ON i1._ IN (0, @M-1)
	-- The neighboring cube face in the selected direction
	INNER JOIN KS.CubeFace AS cf2 ON cf2.CubeFaceID = 
	CASE i1._ WHEN 0 THEN cfn.CubeFaceID_i0 WHEN (@M-1) THEN cfn.CubeFaceID_i1 END
	-- The counter j1 for where we are on the perimeter
	INNER JOIN KS.Counter AS j1 ON j1._ < @M
	-- The SkyPatch of the first face
	INNER JOIN KS.SkyPatch AS sp1 ON
		sp1.CubeFaceID = cf1.CubeFaceID AND
		sp1.i=i1._ AND sp1.j=j1._
	-- The value of j2 will be either 0 or @M-1 depending on the sign of the face
	INNER JOIN KS.Counter AS j2 ON j2._ = (1 + cf1.ci)*@N
	-- The second SkyPatch
	INNER JOIN KS.SkyPatch AS sp2 ON
		sp2.CubeFaceID = cf2.CubeFaceID AND
		sp2.i = sp1.j AND sp2.j = j2._
GROUP BY sp1.SkyPatchID, sp2.SkyPatchID;

-- Insert the second batch of neighbors in the j direction
INSERT INTO KS.SkyPatchNeighbors
(SkyPatchID_1, SkyPatchID_2, dr)
SELECT 
	-- The two SkyPatch cells
 	sp1.SkyPatchID AS SkyPatchID_1,
 	sp2.SkyPatchID AS SkyPatchID_2,
 	-- Distance between these points
	MIN(SQRT(POW(sp2.x-sp1.x,2)+POW(sp2.y-sp1.y,2)+POW(sp2.z-sp1.z,2))) AS dr
FROM 
	-- Start with one cube face
	KS.CubeFace AS cf1
	-- Neighbors of this CubeFace
	INNER JOIN KS.CubeFaceNeighbor AS cfn ON cfn.CubeFaceID = cf1.CubeFaceID
	-- The value of j whose neighbors are being found; 0 for left (i0) and 1 for right (i1)
	INNER JOIN KS.Counter AS j1 ON j1._ IN (0, @M-1)
	-- The neighboring cube face in the selected direction
	INNER JOIN KS.CubeFace AS cf2 ON cf2.CubeFaceID = 
	CASE j1._ WHEN 0 THEN cfn.CubeFaceID_j0 WHEN (@M-1) THEN cfn.CubeFaceID_j1 END
	-- The counter i1 for where we are on the perimeter
	INNER JOIN KS.Counter AS i1 ON i1._ < @M
	-- The SkyPatch of the first face
	INNER JOIN KS.SkyPatch AS sp1 ON
		sp1.CubeFaceID = cf1.CubeFaceID AND
		sp1.i=i1._ AND sp1.j=j1._
	-- The value of i2 will be either 0 or @M-1 depending on the sign of the face
	INNER JOIN KS.Counter AS i2 ON i2._ = (1 + cf1.ci)*@N
	-- The second SkyPatch
	INNER JOIN KS.SkyPatch AS sp2 ON
		sp2.CubeFaceID = cf2.CubeFaceID AND
		sp2.i = i2._ AND sp2.j = i1._
GROUP BY sp1.SkyPatchID, sp2.SkyPatchID;

-- Temporary table for candidate paths spanning a neighbor
-- This preliminary table has no primary key, just an index; the same path can be added multiple times
CREATE OR REPLACE TABLE KS.SkyPatchDistance_batch(
	SkyPatchID_1 INT NOT NULL,
	SkyPatchID_2 INT NOT NULL,
	dr_mid DOUBLE NOT NULL,
	INDEX IDX_SkyPatchID_1_2 (SkyPatchID_1, SkyPatchID_2)
)
ENGINE='Aria' TRANSACTIONAL=0;

INSERT INTO KS.SkyPatchDistance_batch
(SkyPatchID_1, SkyPatchID_2, dr_mid)
SELECT
	-- The two SkyPatch cells identified
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	SQRT(POW(p2.x-p1.x,2) + POW(p2.y-p1.y,2) + POW(p2.z-p1.z,2)) AS dr_mid
FROM
	-- Start with neighbors spanning a cube face
	KS.SkyPatchNeighbors AS spn
	-- The two "anchor" SkyPatch cells
	INNER JOIN KS.SkyPatch AS ap1 ON ap1.SkyPatchID = spn.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS ap2 ON ap2.SkyPatchID = spn.SkyPatchID_2
	-- Search along the first cube face up to dr_max
	INNER JOIN KS.SkyPatchGridDistance AS spg1 ON 
		spg1.i1=ap1.i AND spg1.j1=ap1.j AND	spg1.dr_mid < dr_max
	-- Search along the second cube face up to dr_max
	INNER JOIN KS.SkyPatchGridDistance AS spg2 ON 
		spg2.i1=ap2.i AND spg2.j1=ap2.j AND	spg2.dr_mid < dr_max
	-- The two SkyPatch cells that have been connected
	INNER JOIN KS.SkyPatch AS p1 ON
		p1.CubeFaceID = ap1.CubeFaceID AND p1.i=spg1.i2 AND p1.j=spg1.j2
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = ap2.CubeFaceID AND p2.i=spg2.i2 AND p2.j=spg1.j2;
	
-- Temporary table for candidate paths; now each path only added once
CREATE OR REPLACE TABLE KS.SkyPatchDistance_cand(
	SkyPatchID_1 INT NOT NULL,
	SkyPatchID_2 INT NOT NULL,
	dr_mid DOUBLE NOT NULL,
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2)
)
ENGINE='Aria' TRANSACTIONAL=0;	

INSERT INTO KS.SkyPatchDistance_cand
(SkyPatchID_1, SkyPatchID_2, dr_mid)
SELECT
	spdb.SkyPatchID_1, 
	spdb.SkyPatchID_2,
	MIN(spdb.dr_mid) AS dr_mid
FROM
	KS.SkyPatchDistance_batch AS spdb
WHERE spdb.dr_mid < (dr_max + 1.0 / @N)
GROUP BY spdb.SkyPatchID_1, spdb.SkyPatchID_2;

-- Insert the candidates provisionally into SkyPatchDistance; still need to compute dr_min later
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min, IsCrossFace)
SELECT
	spdc.SkyPatchID_1,
	spdc.SkyPatchID_2,
	spdc.dr_mid,
	-1.0 AS dr_min,
	TRUE AS IsCrossFace
FROM
	KS.SkyPatchDistance_cand AS spdc;

-- Clean up the temporary tables
/*
DROP TABLE KS.SkyPatchNeighbors;
DROP TABLE KS.SkyPatchDistance_batch;
DROP TABLE KS.SkyPatchDistance_cand;
*/

END
$$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchDistance_min()
COMMENT "Calculate the minimum distance field accurately by querying over 4 choices for each corner"
BEGIN 

-- Staging table for the distances from one corner to another (16 rows per pair of SkyPatch cells)
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_min_cand (
	SkyPatchID_1 INT NOT NULL,
	SkyPatchID_2 INT NOT NULL,
	cr1 TINYINT NOT NULL,
	cr2 TINYINT NOT NULL,
	dr_mid DOUBLE NOT NULL,
	x1 DOUBLE NOT NULL,
	y1 DOUBLE NOT NULL,
	z1 DOUBLE NOT NULL,
	x2 DOUBLE NOT NULL,
	y2 DOUBLE NOT NULL,
	z2 DOUBLE NOT NULL,
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2, cr1, cr2)
);

-- Populate all 16 candidate distances for each pair of cells
INSERT INTO KS.SkyPatchDistance_min_cand
(SkyPatchID_1, SkyPatchID_2, cr1, cr2, dr_mid, x1, y1, z1, x2, y2, z2)
SELECT
	-- The two SkyPatch cells
	spd.SkyPatchID_1,
	spd.SkyPatchID_2,
 	-- The selected corner
 	cr1._ AS cr1,
 	cr2._ AS cr2,
	-- Midpoint distance
	SQRT(POW(p2.x-p1.x, 2)+POW(p2.y-p1.y, 2)+POW(p2.z-p1.z, 2)) AS dr_mid,
 	-- Selected corner from SkyPatch 1
 	CASE cr1._ WHEN 0 THEN p1.x00 WHEN 1 THEN p1.x01 WHEN 2 THEN p1.x10	WHEN 3 THEN p1.x11 END AS x1,
 	CASE cr1._ WHEN 0 THEN p1.y00 WHEN 1 THEN p1.y01 WHEN 2 THEN p1.y10	WHEN 3 THEN p1.y11 END AS y1,
 	CASE cr1._ WHEN 0 THEN p1.z00 WHEN 1 THEN p1.z01 WHEN 2 THEN p1.z10	WHEN 3 THEN p1.z11 END AS z1,
 	-- Selected corner from SkyPatch 2
 	CASE cr2._ WHEN 0 THEN p2.x00 WHEN 1 THEN p2.x01 WHEN 2 THEN p2.x10	WHEN 3 THEN p2.x11 END AS x2,
 	CASE cr2._ WHEN 0 THEN p2.y00 WHEN 1 THEN p2.y01 WHEN 2 THEN p2.y10	WHEN 3 THEN p2.y11 END AS y2,
 	CASE cr2._ WHEN 0 THEN p2.z00 WHEN 1 THEN p2.z01 WHEN 2 THEN p2.z10	WHEN 3 THEN p2.z11 END AS z2
FROM
	-- Entries on the distance table that we need compute the distances for
	KS.SkyPatchDistance AS spd
	-- The two SkyPatch cells
	INNER JOIN KS.SkyPatch AS p1 ON p1.SkyPatchID = spd.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS p2 ON p2.SkyPatchID = spd.SkyPatchID_2
	-- Counter to choose the corner of grid cell 1
	INNER JOIN KS.Counter AS cr1 ON cr1._ < 4
	-- Counter to choose the corner of grid cell 2
	INNER JOIN KS.Counter AS cr2 ON cr2._ < 4
WHERE
	-- Only need to update distance calculations on pairs that span different faces
	spd.IsCrossFace=True;

-- Staging table for the distances
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_min LIKE KS.SkyPatchDistance;
ALTER TABLE KS.SkyPatchDistance_min DROP COLUMN IsCrossFace;

-- Calculate minimum distance: group by pair of SkyPatch cells and take minimum
INSERT INTO KS.SkyPatchDistance_min
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
SELECT
	b.SkyPatchID_1,
	b.SkyPatchID_2,
	MIN(b.dr_mid) AS dr_mid,
	-- The minimum distance is the minimum over all 16 choices of the corners
	SQRT(MIN(POW(b.x2-b.x1, 2) + POW(b.y2-b.y1, 2) + POW(b.z2-b.z1, 2))) AS dr_min
FROM
	KS.SkyPatchDistance_min_cand AS b
GROUP BY b.SkyPatchID_1, b.SkyPatchID_2;

-- Apply the distances to the main SkyPatchDistance table
UPDATE
	KS.SkyPatchDistance_min AS spdm
	INNER JOIN  KS.SkyPatchDistance AS spd ON
		spd.SkyPatchID_1 = spdm.SkyPatchID_1 AND
		spd.SkyPatchID_2 = spdm.SkyPatchID_2
SET
	spd.dr_mid = spdm.dr_mid,
	spd.dr_min = spdm.dr_min;

-- Clean up the temporary tables
/*
DROP TEMPORARY TABLE KS.SkyPatchDistance_min_cand;
DROP TEMPORARY TABLE KS.SkyPatchDistance_min;
*/
END $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchDistance(
	IN dr_max DOUBLE)
COMMENT "Populate the KS.SkyPatchDistance table"
BEGIN 

-- Empty table first so this procedure can be re-run on demand
TRUNCATE TABLE KS.SkyPatchDistance;

-- Most of the neighbors are on the same major cube face
CALL KS.MakeTable_SkyPatchDistance_face(dr_max);
SELECT 'Added distance for pairs on the same face.' AS msg;

-- Extend to include the paths that span over faces
CALL KS.MakeTable_SkyPatchDistance_extend(dr_max);
SELECT 'Added provisional distance for pairs spanning faces.' AS msg;

-- Calculate the correct midpoint and minimum distance
CALL KS.MakeTable_SkyPatchDistance_min();
SELECT 'Calculated minimum distance.' AS msg;

-- Delete rows that are above the maximum distance
DELETE  
FROM KS.SkyPatchDistance
WHERE dr_min >=dr_max;

END
$$

-- *********************************************************************************
DELIMITER ;

-- ************************************************************************************************
/*
 * Test 
CREATE OR REPLACE TABLE SkyPatchDistance2 LIKE SkyPatchDistance;
INSERT INTO SkyPatchDistance2  SELECT * FROM SkyPatchDistance;

CALL MakeTable_SkyPatchDistance_bf(0.25);

SELECT 
	spd2.SkyPatchID_1,
	spd2.SkyPatchID_2,
	spd2.dr_mid,
	spd.dr_mid AS dr_mid_bf,
	(spd2.dr_mid - spd.dr_mid)*1000 dr_diff
FROM 
	SkyPatchDistance2  AS spd2
	LEFT JOIN KS.SkyPatchDistance AS spd ON
		spd.SkyPatchID_1 = spd2.SkyPatchID_1 AND
		spd.SkyPatchID_2 = spd2.SkyPatchID_2	
ORDER BY dr_diff desc;
*/
