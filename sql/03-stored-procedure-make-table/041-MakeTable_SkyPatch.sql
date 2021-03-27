-- ************************************************************************************************
-- Global variables used at bottom to build the tables
-- Set the variable N
SET @N = 128;
-- SET @N = POW(2,10);

-- Maximum distance for SkyPatchGridDistance
SET @dr_max = GREATEST(1.0 / @N, 1.0/360);

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
	-- The starting SkyPatch
	KS.SkyPatch AS p1
	-- All the neighbors of this grid cell
	INNER JOIN KS.SkyPatchGridDistance AS gd ON
		gd.i1=p1.i AND gd.j1=p1.j
	-- The second SkyPatch
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = p1.CubeFaceID AND
		p2.i=gd.i2 AND p2.j=gd.j2
WHERE
	gd.dr_mid < dr_max;

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
	
-- Staging table for pairs of SkyPatch cells that are connected by edges spanning cube faces
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_edge LIKE KS.SkyPatchDistance;
-- Need to drop PK because this query will generate the same candidate multiple times
ALTER TABLE KS.SkyPatchDistance_edge DROP PRIMARY KEY;
ALTER TABLE KS.SkyPatchDistance_edge DROP dr_min;
ALTER TABLE KS.SkyPatchDistance_edge DROP IsCrossFace;

-- Generate the neighbors connected by an edge
INSERT INTO KS.SkyPatchDistance_edge
(SkyPatchID_1, SkyPatchID_2, dr_mid)
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
	SQRT(POW(p2.x-p1.x, 2)+POW(p2.y-p1.y, 2)+POW(p2.z-p1.z, 2)) AS dr_mid
FROM
	-- Start with one cube face
	KS.CubeFace AS cf
	-- Choose which coordinate we are on the perimeter of
	-- The coordinate index for the perimeter is 1 (u) or 2 (v)
	INNER JOIN KS.Counter AS pk ON (pk._ BETWEEN 1 AND 2)	
	-- Start with the two values of an index that are on the perimeter
	-- Two possible settings for i or j to lie on the perimeter of a face
	INNER JOIN KS.Counter AS pv ON (pv._ = 0 OR pv._ = @M-1)
	-- Single counter for the four possibilities in the order i0, i1, j0, j1
	INNER JOIN KS.Counter AS pn ON pn._ = 2*(pk._-1) + IF(pv._ > 0, 2, 1)
	-- SkyPatch cells matching this CubeFace and on the perimeter
	INNER JOIN KS.SkyPatch AS p1 ON	
		p1.CubeFaceID = cf.CubeFaceID AND
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
)
SELECT
	t1.SkyPatchID_1,
	t1.SkyPatchID_2,
	t1.dr_mid
FROM
	t1
WHERE
	t1.dr_mid < LEAST(1.0 / @N, dr_max);

-- New staging table with just the unique pairs that are candidates to be added
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_cand LIKE KS.SkyPatchDistance_edge;
ALTER TABLE KS.SkyPatchDistance_cand ADD PRIMARY KEY (SkyPatchID_1, SkyPatchID_2);

INSERT INTO KS.SkyPatchDistance_cand
(SkyPatchID_1, SkyPatchID_2, dr_mid)
SELECT
	spdb.SkyPatchID_1, 
	spdb.SkyPatchID_2,
	min(spdb.dr_mid) AS dr_mid
FROM
	KS.SkyPatchDistance_edge AS spdb
WHERE
	spdb.dr_mid < dr_max
GROUP BY spdb.SkyPatchID_1, spdb.SkyPatchID_2;

-- Generate candidate paths from one SkyPatch to another within the maximum distance
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_batch LIKE KS.SkyPatchDistance;
ALTER TABLE KS.SkyPatchDistance_batch DROP PRIMARY KEY;
ALTER TABLE KS.SkyPatchDistance_batch DROP dr_min;
ALTER TABLE KS.SkyPatchDistance_batch DROP IsCrossFace;

INSERT INTO KS.SkyPatchDistance_batch
(SkyPatchID_1, SkyPatchID_2, dr_mid)
WITH t1 AS (
SELECT
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	SQRT(POW(p2.x-p1.x,2)+POW(p2.y-p1.y,2)+POW(p2.z-p1.Z,2)) AS dr_anchor
FROM
	-- Start with the batch of candidates
	KS.SkyPatchDistance_cand AS spd
	-- Join the two SkyPatch cells
	INNER JOIN KS.SkyPatch AS p1 ON p1.SkyPatchID = spd.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS p2 ON p2.SkyPatchID = spd.SkyPatchID_2
)
SELECT
	p1.SkyPatchID AS SkyPatchID_1,
	p2.SkyPatchID AS SkyPatchID_2,
	(t1.dr_anchor + gd1.dr_mid + gd2.dr_mid) AS dr_mid
FROM
	t1
	-- The two anchor SkyPatch cells
	INNER JOIN KS.SkyPatch AS p1a ON p1a.SkyPatchID = t1.SkyPatchID_1
	INNER JOIN KS.SkyPatch AS p2a ON p2a.SkyPatchID = t1.SkyPatchID_2
	-- Search for neighbors on the first face within the maximum distance
	INNER JOIN KS.SkyPatchGridDistance AS gd1 ON
		gd1.i1 = p1a.i AND gd1.j1 = p1a.j AND 
		gd1.dr_mid < dr_max
	-- The adjusted SkyPatch 1
	INNER JOIN KS.SkyPatch AS p1 ON
		p1.CubeFaceID = p1a.CubeFaceID AND
		p1.i = gd1.i2 AND p1.j = gd1.j2
	-- Search for neighbors on the second face within the maximum distance
	INNER JOIN KS.SkyPatchGridDistance AS gd2 ON
		gd2.i1 = p2a.i AND gd2.j1 = p2a.j AND 
		-- gd2.dr_mid < dr_max - gd1.dr_mid
		gd2.dr_mid < dr_max
	-- The adjusted SkyPatch 2
	INNER JOIN KS.SkyPatch AS p2 ON
		p2.CubeFaceID = p2a.CubeFaceID AND
		p2.i = gd2.i2 AND p2.j = gd2.j2;

-- Insert this batch into SkyPatchDistance_cand
INSERT INTO KS.SkyPatchDistance_cand
(SkyPatchID_1, SkyPatchID_2, dr_mid)
SELECT
	spdb.SkyPatchID_1, 
	spdb.SkyPatchID_2,
	spdb.dr_mid
FROM
	KS.SkyPatchDistance_batch AS spdb
WHERE 
	NOT EXISTS(
	SELECT spdc.SkyPatchID_1
	FROM KS.SkyPatchDistance_cand AS spdc
	WHERE 
		spdc.SkyPatchID_1 = spdb.SkyPatchID_1 AND
		spdc.SkyPatchID_2 = spdb.SkyPatchID_2
	)
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

END
$$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_SkyPatchDistance_min()
COMMENT "Calculate the minimum distance field accurately by querying over 4 choices for each corner"
BEGIN 

-- Staging table for the distances from one corner to another (16 rows per pair of SkyPatch cells)
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_batch (
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
INSERT INTO KS.SkyPatchDistance_batch
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
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_dist LIKE KS.SkyPatchDistance;
ALTER TABLE KS.SkyPatchDistance_dist DROP COLUMN IsCrossFace;

-- Calculate minimum distance: group by pair of SkyPatch cells and take minimum
INSERT INTO KS.SkyPatchDistance_dist
(SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
SELECT
	b.SkyPatchID_1,
	b.SkyPatchID_2,
	MIN(b.dr_mid) AS dr_mid,
	-- The minimum distance is the minimum over all 16 choices of the corners
	SQRT(MIN(POW(b.x2-b.x1, 2) + POW(b.y2-b.y1, 2) + POW(b.z2-b.z1, 2))) AS dr_min
FROM
	KS.SkyPatchDistance_batch AS b
GROUP BY b.SkyPatchID_1, b.SkyPatchID_2;

-- Apply the distances to the main SkyPatchDistance table
UPDATE
	KS.SkyPatchDistance_dist AS spdd
	INNER JOIN  KS.SkyPatchDistance AS spd ON
		spd.SkyPatchID_1 = spdd.SkyPatchID_1 AND
		spd.SkyPatchID_2 = spdd.SkyPatchID_2
SET
	spd.dr_mid = spdd.dr_mid,
	spd.dr_min = spdd.dr_min;

-- Clean up the temporary tables
DROP TEMPORARY TABLE KS.SkyPatchDistance_batch;
DROP TEMPORARY TABLE KS.SkyPatchDistance_dist;

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

-- Clean up the temporary tables
DROP TEMPORARY TABLE KS.SkyPatchDistance_edge;
DROP TEMPORARY TABLE KS.SkyPatchDistance_batch;
DROP TEMPORARY TABLE KS.SkyPatchDistance_cand;

END
$$

-- *********************************************************************************
DELIMITER ;

-- *********************************************************************************
-- Build all the tables
CALL KS.MakeTable_SkyPatchGrid(@N, @dr_max);
CALL KS.MakeTable_SkyPatch();
--- CALL KS.MakeTable_SkyPatchDistance_face(@dr_max);
CALL KS.MakeTable_SkyPatchDistance(@dr_max);

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