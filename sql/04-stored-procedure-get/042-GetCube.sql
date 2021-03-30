DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetCubeFace()
COMMENT "Get the catalog of cube faces."

BEGIN 

SELECT
	cf.CubeFaceID,
	cf.CubeFaceCD,
	cf.i,
	cf.j1,
	cf.j2,
	cf.ci,
	cf.alpha,
	cf.beta,
	cf.gamma	
FROM
	KS.CubeFace AS cf
ORDER BY cf.CubeFaceID;

END
$$

DELIMITER ;
