SELECT 
	sp.SkyPatchID,
	sp.x,
	sp.y,
	sp.z
FROM 
	KS.SkyPatch AS sp
ORDER BY SkyPatchID	
LIMIT 100;