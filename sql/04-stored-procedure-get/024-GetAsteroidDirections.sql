DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidDirections(
    IN n0 INT,
    IN n1 INT,
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Get the directions from earth to asteroid for a block of asteroids."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;
	
SELECT
    ad.AsteroidID,
	ad.TimeID,
	ad.tObs,
	ad.ux,
	ad.uy,
	ad.uz,
    ad.LightTime,
    -- The time when light left the asteroid
    ad.tObs - ad.LightTime / 1440.0 AS tAst
FROM
	KS.AsteroidDirections AS ad
WHERE 
	ad.TimeID BETWEEN @TimeID_0 AND @TimeID_1 AND
	ad.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY ad.AsteroidID, ad.TimeID;

END
$$

DELIMITER ;
