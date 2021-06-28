DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidSkyPatch(
    IN n0 INT,
    IN n1 INT
)
COMMENT "Get asteroid sky patch data for asteroids in the specified range."

BEGIN 

SELECT
	asp.AsteroidID,
	asp.Segment,
	asp.SkyPatchID,
	asp.TimeID_0,
	asp.TimeID_1
FROM
	KS.AsteroidSkyPatch AS asp
WHERE
	-- Only selected range of asteroids
	asp.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY asp.AsteroidID, asp.Segment;

END
$$


-- ********************************************************************************
DELIMITER ;
