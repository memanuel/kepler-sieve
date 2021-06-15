SELECT * FROM KS.SkyPatchGrid LIMIT 100;

SELECT 
	min(spg.i) AS i_min,
	max(spg.i) AS i_max,
	min(spg.j) AS j_min,
	max(spg.j) AS j_max
FROM
	KS.SkyPatchGrid AS spg;
	

SELECT * FROM KS.SkyPatch WHERE SkyPatchID=6382172;