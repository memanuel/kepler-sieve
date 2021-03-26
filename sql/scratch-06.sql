-- Get M and calculate N
SELECT
MAX(i)+1 AS M
INTO @M
FROM KS.SkyPatchGrid;
SET @N = @M DIV 2;

# Set maximum distance
SET @dr_max = 1.0/64.0;

-- CALL KS.MakeTable_SkyPatchDistance_face(1.0/64.0);
-- CALL KS.MakeTable_SkyPatchDistance_bf(1.0/64.0);
-- CALL KS.MakeTable_SkyPatchDistance(1.0/64.0);

-- Insert a second batch that are neighbors across different faces
-- dr_min is currently a placeholder.  
-- This is *not* exhaustive, but will pick up all junctions from one face to another
-- INSERT INTO KS.SkyPatchDistance
-- (SkyPatchID_1, SkyPatchID_2, dr_mid, dr_min)
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
	t1.dr_mid,
	-1.0 AS dr_min
FROM
	t1
WHERE
	t1.dr_mid < LEAST(1.0 / @N, dr_max);
