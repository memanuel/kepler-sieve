-- ************************************************************************************************
-- Global variables used at bottom to build the tables
-- Set the variable N
SET @N = 4;
-- SET @N = POW(2,9);

-- Maximum distance for SkyPatchGridDistance
SET @dr_max = 0.25;

-- *********************************************************************************
DELIMITER $$

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
SET @gw = CAST(CEILING(dr_max*n) AS INT);

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

-- Set the bounds on u, v, w
-- UPDATE KS.SkyPatchGrid
-- SET
-- 	uMin = LEAST(u00, u01, u10, u11),
-- 	uMax = GREATEST(u00, u01, u10, u11),
-- 	vMin = LEAST(v00, v01, v10, v11),
-- 	vMax = GREATEST(v00, v01, v10, v11),
-- 	wMin = LEAST(w00, w01, w10, w11),
-- 	wMax = GREATEST(w00, w01, w10, w11);
-- Populate SkyPatchGridDistance by querying over 16 possible pairs of points
-- Compare 4 corners on p1 with 4 corners on p2 and take the minimum

INSERT INTO KS.SkyPatchGridDistance
(i1, j1, i2, j2, dr_mid, dr_min)
WITH t1 AS (
SELECT
	g1.i AS i1,
	g1.j AS j1,
	g2.i AS i2,
	g2.j AS j2,
	-- Distance at midpoint
	SQRT(POW(g2.u-g1.u,2)+POW(g2.v-g1.v,2)+POW(g2.w-g1.w,2)) AS dr_mid,
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
 	CASE cr2._ WHEN 0 THEN g2.w00 WHEN 1 THEN g2.w01 WHEN 2 THEN g2.w10	WHEN 3 THEN g2.w11 END AS w2
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
	INNER JOIN KS.Counter AS cr1 ON cr1._ < 4
	-- Counter to choose the corner of grid cell 2
	INNER JOIN KS.Counter AS cr2 ON cr2._ < 4
), t2 AS(
SELECT
	t1.i1,
	t1.j1,
	t1.i2,
	t1.j2,
	t1.cr1,
	t1.cr2,
	t1.dr_mid,
	SQRT(MIN(POW(u2-u1,2)+POW(v2-v1,2)+POW(w2-w1,2))) AS dr_min
FROM
	t1
GROUP BY t1.i1, t1.j1, t1.i2, t1.j2
)
SELECT
	t2.i1,
	t2.j1,
	t2.i2,
	t2.j2,
	dr_mid,
	t2.dr_min
FROM
	t2
WHERE
	dr_min < dr_max;

END $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatch()
COMMENT "Populate the KS.SkyPatchGrid and KS.SkyPatchGridDistance tables"
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
PROCEDURE KS.MakeTable_SkyPatchDistance_bf(
	IN dr_max DOUBLE)
COMMENT "Populate the KS.SkyPatchDistance table"
BEGIN 

-- Brute force approach - does not scale, but useful for checking with small N
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_max)
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
	t1.dr_mid AS dr_max
FROM
	t1
WHERE 
	t1.dr_mid < dr_max;
	
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

-- Calculate @M and @N
SELECT 
	MAX(spg.i)+1 AS M 
INTO @M 
FROM KS.SkyPatchGrid AS spg;
SET @N = (@M DIV 2);

-- Most of the neighbors are on the same major cube face
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
SELECT
	-- The two SkyPatch cells in this pair
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	-- The midpoint and minimum distance
	gd.dr_mid,
	gd.dr_min
FROM
	-- The starting SkyPatch
	KS.SkyPatch AS p1
	-- All the neighbors of this grid cell
	INNER JOIN KS.SkyPatchGridDistance AS gd ON
		gd.i1=p1.i AND gd.j1=p1.j
	-- The second SkyPatch
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = p1.CubeFaceID AND
		p2.i=gd.i2 AND p2.j=gd.j2;

-- Insert a second batch that are neighbors across different faces
-- dr_min is currently a placeholder.  
-- This is *not* exhaustive, but will pick up all junctions from one face to another
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
WITH t1 AS (
SELECT
	-- The two SkyPatch cells
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	-- The grid cell of each SkyPatch
	p1.i AS i1,
	p1.j AS j1,
	p2.i AS i2,
	p2.j AS j2,
	-- Midpoint distance
	SQRT(POW(p2.x-p1.x, 2)+POW(p2.y-p1.y, 2)+POW(p2.z-p1.z, 2)) AS dr_mid,
	-- Midpoint of SkyPatch 1
	p1.x AS x1,
	p1.y AS y1,
	p1.z AS z1,
	-- Midpoint of SkyPatch 2
	p2.x AS x2,
	p2.y AS y2,
	p2.z AS z2,
	-- The selected perimeter axis (i or j)
	pk._ AS pk,
	-- The perimeter value (0 or 2N-1)
	pv._ AS pv,
	-- Counter of the four permiters in order i0, i1, j0, j1
	pn._ AS pn
FROM
	-- Choose which coordinate we are on the perimeter of
	KS.Counter AS pk
	-- Start with the two values of an index that are on the perimeter
	CROSS JOIN KS.Counter AS pv		
	-- Single counter for the four possibilities in the order i0, i1, j0, j1
	INNER JOIN KS.Counter AS pn ON
		pn._ = 2*(pk._-1) + IF(pv._ > 0, 2, 1)
	-- SkyPatch cells on the perimeter
	INNER JOIN KS.SkyPatch AS p1 ON	
		(p1.i=pv._ AND pk._ = 1) OR
		(p1.j=pv._ AND pk._ = 2)
	-- Neighbors of this cube face
	INNER JOIN KS.CubeFaceNeighbor AS cfn ON cfn.CubeFaceID = p1.CubeFaceID
	-- Relevant neighbor of this cube face
	INNER JOIN KS.CubeFace AS cf2 ON cf2.CubeFaceID = 
		CASE pn._ 
			WHEN 1 THEN cfn.CubeFaceID_i0 
			WHEN 2 THEN cfn.CubeFaceID_i1 
			WHEN 3 THEN cfn.CubeFaceID_j0 
			WHEN 4 THEN cfn.CubeFaceID_j1 
			END
	-- Wrap to the relevant grid cell
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = cf2.CubeFaceID AND
 		p2.i IN (p1.i, p1.j, 0, @M-1) AND
 		p2.j IN (p1.i, p1.j, 0, @M-1)			
WHERE
	-- The coordinate index for the perimeter is 1 (u) or 2 (v)
	(pk._ BETWEEN 1 AND 2) AND
	-- Two possible settings for i or j to lie on the perimeter of a face
	(pv._ = 0 OR pv._ = @M-1)
)
SELECT
	t1.SkyPatchID_1,
	t1.SkyPatchID_2,
	t1.dr_mid AS dr_min,
	t1.dr_mid
FROM
	t1
WHERE
	t1.dr_mid < (1.0 / @N);

-- Create a third batch by searching along both grid faces; save in a temp table
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_batch LIKE KS.SkyPatchDistance;
-- Need to drop PK because this query will generate the same candidate multiple times
ALTER TABLE KS.SkyPatchDistance_batch DROP PRIMARY KEY;
ALTER TABLE KS.SkyPatchDistance_batch DROP dr_mid;
ALTER TABLE KS.SkyPatchDistance_batch DROP dr_min;
-- Upper bound on dr_mid
ALTER TABLE KS.SkyPatchDistance_batch ADD COLUMN dr_mid_ub DOUBLE NOT NULL;

-- Generate candidate paths from one SkyPatch to another within the maximum distance
INSERT INTO KS.SkyPatchDistance_batch
(SkyPatchID_1, SkyPatchID_2, dr_mid_ub)
WITH t1 AS (
SELECT
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	SQRT(POW(p2.x-p1.x,2)+POW(p2.y-p1.y,2)+POW(p2.z-p1.Z,2)) AS dr_anchor
FROM
	-- Start with entries on the distance table
	KS.SkyPatchDistance AS spd
	-- Join the two SkyPatch cells
	INNER JOIN KS.SkyPatch AS p1 ON p1.SkyPatchID = spd.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS p2 ON p2.SkyPatchID = spd.SkyPatchID_2
WHERE
	-- Only take pairs spanning different cube faces
	p1.CubeFaceID <> p2.CubeFaceID
)
SELECT
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	t1.dr_anchor + gd1.dr_mid + gd1.dr_mid AS dr_mid
FROM
	t1
	-- The two anchor SkyPatch cells
	INNER JOIN KS.SkyPatch AS p1a ON p1a.SkyPatchID = t1.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS p2a ON p2a.SkyPatchID = t1.SkyPatchID_2
	-- Search for neighbors on the first face within the maximum distance
	INNER JOIN KS.SkyPatchGridDistance AS gd1 ON
		gd1.i1 = p1a.i AND gd1.j1 = p1a.j AND 
		gd1.dr_mid < @dr_max - t1.dr_anchor
	-- The adjusted SkyPatch 1
	INNER JOIN KS.SkyPatch AS p1 ON
		p1.CubeFaceID = p1a.CubeFaceID AND
		p1.i = gd1.i2 AND p1.j = gd1.j2
	-- Search for neighbors on the second face within the maximum distance
	INNER JOIN KS.SkyPatchGridDistance AS gd2 ON
		gd2.i1 = p2a.i AND gd2.j1 = p2a.j AND 
		gd2.dr_mid < @dr_max - t1.dr_anchor - gd1.dr_mid
	-- The adjusted SkyPatch 2
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = p2a.CubeFaceID AND
		p2.i = gd2.i2 AND p2.j = gd2.j2;

/*
INSERT INTO KS.SkyPatchDistance_batch
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
WITH t1 AS(
SELECT
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	SQRT(POW(p2.x-p1.x,2)+POW(p2.y-p1.y,2)+POW(p2.z-p1.Z,2)) AS dr_mid
FROM
	-- Start with entries on the distance table
	KS.SkyPatchDistance AS spd
	-- Join the two SkyPatch cells
	INNER JOIN KS.SkyPatch AS p1_anchor ON p1_anchor.SkyPatchID = spd.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS p2_anchor ON p2_anchor.SkyPatchID = spd.SkyPatchID_2
	-- Search for additional nearby connections on cube face 1
	INNER JOIN KS.CounterSigned AS di1 ON di1._ BETWEEN -@grid_width AND @grid_width
	INNER JOIN KS.CounterSigned AS dj1 ON dj1._ BETWEEN -@grid_width AND @grid_width
	INNER JOIN KS.SkyPatch AS p1 ON
		p1.CubeFaceID = p1_anchor.CubeFaceID AND
		p1.i = p1.i + di1._ AND p1.j = p1_anchor.j + dj1._
	-- Search for additional nearby connections on cube face 2
	INNER JOIN KS.CounterSigned AS di2 ON di2._ BETWEEN -@grid_width AND @grid_width
	INNER JOIN KS.CounterSigned AS dj2 ON dj2._ BETWEEN -@grid_width AND @grid_width
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = p2_anchor.CubeFaceID AND
		p2.i = p2.i + di2._ AND p2.j = p2_anchor.j + dj2._
WHERE
	-- Only take pairs spanning different cube faces
	p1_anchor.CubeFaceID <> p2_anchor.CubeFaceID
)
SELECT
	t1.SkyPatchID_1,
	t1.SkyPatchID_2,
	t1.dr_mid,
	t1.dr_mid AS dr_min
FROM
	t1;

-- Insert the third batch into the distance table
-- Only insert rows not already present, and no duplicates
INSERT INTO KS.SkyPatchDistance
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
SELECT 
	spdb.SkyPatchID_1,
	spdb.SkyPatchID_2,
	spdb.dr_mid,
	spdb.dr_min
FROM 
	KS.SkyPatchDistance_batch AS spdb
WHERE NOT EXISTS (
	SELECT spd.SkyPatchID_1	
	FROM KS.SkyPatchDistance AS spd
	WHERE 
		spd.SkyPatchID_1 = spdb.SkyPatchID_1 AND 
		spd.SkyPatchID_2 = spdb.SkyPatchID_2
)
GROUP BY spdb.SkyPatchID_1, spdb.SkyPatchID_2;
*/

END
$$

DELIMITER ;

-- *********************************************************************************
-- Build all the tables
-- ************************************************************************************************
CALL KS.MakeTable_SkyPatchGrid(@N, @dr_max);
CALL KS.MakeTable_SkyPatch();
CALL KS.MakeTable_SkyPatchDistance(@dr_max);
