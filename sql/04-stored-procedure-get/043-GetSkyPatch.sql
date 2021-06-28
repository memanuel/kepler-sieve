DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetSkyPatch()
COMMENT "Get the catalog of SkyPatch cells."

BEGIN 

SELECT
	sp.SkyPatchID,
	sp.CubeFaceID-1 AS f,
	sp.i,
	sp.j,
	sp.x,
	sp.y,
	sp.z
FROM
	KS.SkyPatch AS sp
ORDER BY SkyPatchID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetSkyPatchNeighbors()
COMMENT "Get the catalog of neighbors of each SkyPatch cell."

BEGIN 

SELECT
	spn.SkyPatchID_1,
	spn.SkyPatchID_2,
	spn.dr_mid
FROM
	KS.SkyPatchNeighbor AS spn
ORDER BY SkyPatchID_1, SkyPatchID_2;

END
$$

-- ********************************************************************************
DELIMITER ;
