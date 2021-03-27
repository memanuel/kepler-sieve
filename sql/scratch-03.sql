-- Calculate @M and @N
SELECT 
	MAX(spg.i)+1 AS M 
INTO @M 
FROM KS.SkyPatchGrid AS spg;
SET @N = (@M DIV 2);

-- Staging table for pairs of SkyPatch cells that span cube faces
CREATE OR REPLACE TEMPORARY TABLE KS.SkyPatchDistance_edge LIKE KS.SkyPatchDistance;
-- Need to drop PK because this query will generate the same candidate multiple times
ALTER TABLE KS.SkyPatchDistance_edge DROP PRIMARY KEY;
ALTER TABLE KS.SkyPatchDistance_edge DROP dr_min;
ALTER TABLE KS.SkyPatchDistance_edge DROP IsCrossFace;

-- Insert a second batch that are neighbors across different faces
-- This is *not* exhaustive, but will pick up all junctions from one face to another
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
	t1.dr_mid < LEAST(1.0 / @N, @dr_max);

-- New staging table just the unique pairs that are candidates to be added
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
	spdb.dr_mid < @dr_max
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
		gd1.dr_mid < @dr_max
	-- The adjusted SkyPatch 1
	INNER JOIN KS.SkyPatch AS p1 ON
		p1.CubeFaceID = p1a.CubeFaceID AND
		p1.i = gd1.i2 AND p1.j = gd1.j2
	-- Search for neighbors on the second face within the maximum distance
	INNER JOIN KS.SkyPatchGridDistance AS gd2 ON
		gd2.i1 = p2a.i AND gd2.j1 = p2a.j AND 
		-- gd2.dr_mid < dr_max - gd1.dr_mid
		gd2.dr_mid < @dr_max
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

-- Calculate the correct midpoint and minimum distance
CALL KS.MakeTable_SkyPatchDistance_min();

-- Delete rows that are above the maximum distance
DELETE  
FROM KS.SkyPatchDistance
WHERE dr_min >= @dr_max;

-- Clean up the temporary tables
DROP TEMPORARY TABLE KS.SkyPatchDistance_edge;
DROP TEMPORARY TABLE KS.SkyPatchDistance_batch;
DROP TEMPORARY TABLE KS.SkyPatchDistance_cand;
