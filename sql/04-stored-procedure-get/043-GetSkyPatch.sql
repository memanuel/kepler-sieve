DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetSkyPatch()
COMMENT "Get the catalog of SkyPatch cells."

BEGIN 

SELECT
	sp.SkyPatchID,
--	sp.CubeFaceID-1 AS f,	
--	sp.i,
--	sp.j,
	sp.x,
	sp.y,
	sp.z
FROM
	KS.SkyPatch AS sp
ORDER BY SkyPatchID;

END
$$

DELIMITER ;
