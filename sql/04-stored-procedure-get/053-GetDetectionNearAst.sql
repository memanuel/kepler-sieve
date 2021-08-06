DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionNearAst(
    IN AsteroidID_0 INT,
    IN AsteroidID_1 INT
)
COMMENT "Get detections near asteroids in the given range of asteroids."

BEGIN 

SELECT
	-- Key fields: DetectionID and AsteroidID
	dna.DetectionID,
	dna.AsteroidID,
	-- Time of this detection
	det.mjd AS tObs,
	-- The position of the observatory
	dt.qObs_x,
	dt.qObs_y,
	dt.qObs_z,
	-- The direction of this detection
	det.ux AS uObs_x,
	det.uy AS uObs_y,
	det.uz AS uObs_z
FROM
	-- Start with the detection candidates populated by C++ program detection_near_ast.cpp
	KS.DetectionNearAsteroid AS dna
	-- Join Detection to get the position
	INNER JOIN KS.Detection AS det ON det.DetectionID = dna.DetectionID
	-- Join DetectionTime to get the observatory position
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = det.DetectionTimeID
WHERE
	dna.AsteroidID BETWEEN AsteroidID_0 AND (AsteroidID_1-1)
ORDER BY dna.AsteroidID, dna.DetectionID;

END
$$

-- ********************************************************************************
DELIMITER ;
