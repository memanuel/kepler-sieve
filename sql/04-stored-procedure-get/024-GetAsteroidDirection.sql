DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidDirection(
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
	KS.AsteroidDirection AS ad
WHERE 
	ad.TimeID BETWEEN @TimeID_0 AND @TimeID_1 AND
	ad.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY ad.AsteroidID, ad.TimeID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.TestAsteroidDirection(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Test the asteroid directions in the DB by comparing to JPL directions."

BEGIN 

-- Compute TimeID from dates
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;
	
SELECT
	-- Key fields
	adj.AsteroidID,
	adj.TimeID,
	adj.mjd,
	-- The direction and light time according to JPL
	adj.ux_ast AS ux_jpl,
	adj.uy_ast AS uy_jpl,
	adj.uz_ast AS uz_jpl,
	adj.LightTime AS LightTime_jpl,
	-- The direction and light time according to MSE
	ad.ux AS ux_mse,
	ad.uy AS uy_mse,
	ad.uz AS uz_mse,
	ad.LightTime AS LightTime_mse
FROM
	-- Start with JPL asteroid directions
	JPL.AsteroidDirection AS adj
	-- The observatory
	INNER JOIN KS.Observatory AS obs ON obs.ObservatoryID = adj.ObservatoryID
	-- The MSE direction for this asteroid and time
	INNER JOIN KS.AsteroidDirection AS ad ON
		ad.AsteroidID = adj.AsteroidID AND
		ad.TimeID = adj.TimeID
WHERE
	-- Only geocenter 
	obs.ObservatoryShortName = 'Geocenter'
	-- Only selected dates
	AND ad.TimeID BETWEEN @TimeID_0 AND @TimeID_1
ORDER BY adj.AsteroidID, adj.TimeID;

END
$$

DELIMITER ;
