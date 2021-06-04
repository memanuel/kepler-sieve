DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidDirection_Import(
    IN n0 INT,
    IN n1 INT
)
COMMENT "Get asteroid directions imported from Horizons."

BEGIN 

-- Get the asteroid directions from the JPL import
SELECT
	adi.AsteroidID,
	it.TimeID,
	it.mjd,
	adi.RA,
	adi.`DEC`,
	adi.dRAxCosD_dt,
	adi.dDEC_dt,
	adi.Mag,
	adi.Brightness,
	adi.r, 
	adi.rDot,
	adi.delta,
	adi.deltaDot,
	adi.LightTime
FROM
	JPL.AsteroidDirectionImport AS adi
	INNER JOIN KS.IntegrationTime AS it ON it.MJD = (adi.JD - 2400000.5)
WHERE
    adi.AsteroidID BETWEEN n0 AND (n1-1);

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidDirection(
    IN n0 INT,
    IN n1 INT
)
COMMENT "Get asteroid directions imported from Horizons."

BEGIN 

-- Get the asteroid directions from the JPL import
SELECT
	-- Key fields
	ad.AsteroidID,
	ad.TimeID,
	ad.mjd,
	-- Original RA/DEC
	ad.RA,
	ad.`DEC`,
	-- Calculated direction
	ad.ux,
	ad.uy,
	ad.uz,
	-- Range
	ad.r, 
	ad.rDot,
	ad.delta,
	ad.deltaDot,
	-- Misc
	ad.LightTime,
	ad.Mag
FROM
	JPL.AsteroidDirection AS ad
	INNER JOIN JPL.Asteroid
WHERE
    ad.AsteroidID BETWEEN n0 AND (n1-1);

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.AsteroidDirectionTest()
COMMENT "Get asteroid directions and positions to test MSE calculations."

BEGIN 

SELECT
	-- Key fields
	ad.AsteroidID,
	ad.TimeID,
	ad.mjd AS tAst,
	-- The asteroid position
	hv.qx AS qAst_x,
	hv.qy AS qAst_y,
	hv.qz AS qAst_z,
	-- Light time
	ad.LightTime,
	-- Calculated direction
	ad.ux,
	ad.uy,
	ad.uz
FROM
	JPL.AsteroidDirection AS ad
	INNER JOIN JPL.HorizonsBody AS hb ON hb.SmallBodyID = ad.AsteroidID
	INNER JOIN JPL.HorizonsVectors AS hv ON 
		hv.HorizonsBodyID=hb.HorizonsBodyID AND
		hv.TimeID=ad.TimeID
WHERE 
	ad.TimeID BETWEEN 48000*24*60 AND 63000*24*60
ORDER BY ad.AsteroidID, ad.TimeID;

END
$$

-- ********************************************************************************
DELIMITER ;
