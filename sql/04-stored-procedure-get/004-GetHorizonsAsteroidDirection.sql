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
	adi.ObservatoryID,
	it.TimeID,
	it.mjd,
	adi.RA_ast,
	adi.DEC_ast,
	adi.RA_app,
	adi.DEC_app,
	adi.LightTime,
	adi.r,
	adi.rDot,
	adi.delta,
	adi.deltaDot
FROM
	JPL.AsteroidDirectionImport AS adi
	INNER JOIN KS.IntegrationTime AS it ON it.MJD = (adi.JD - 2400000.5)
WHERE
    adi.AsteroidID BETWEEN n0 AND (n1-1)
ORDER BY adi.AsteroidID, adi.ObservatoryID, it.TimeID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.GetAsteroidDirection(
	IN ObservatoryShortName VARCHAR(16),
    IN n0 INT,
    IN n1 INT
)
COMMENT "Get asteroid directions imported from Horizons. ObservatoryShortName one of 'Geocenter' and 'ZTF'"

BEGIN 

-- Get the asteroid directions from JPL, with augmented direction components added by MSE
SELECT
	-- Key fields
	ad.AsteroidID,
	ad.ObservatoryID,
	ad.TimeID,
	-- Departure time of light
	ad.mjd - ad.LightTime / 1440.0 AS tAst,
	-- Quoted astrometric RA/DEC
	ad.RA_ast,
	ad.DEC_ast,
	-- Astrometric direction vector
	ad.ux_ast,
	ad.uy_ast,
	ad.uz_ast,
	-- Quoted apparent RA/DEC
	ad.RA_app,
	ad.DEC_app,
	-- Apparent direction vector
	ad.ux_app,
	ad.uy_app,
	ad.uz_app,
	-- Light time and observation time
	ad.LightTime,
	ad.mjd AS tObs
FROM
	-- The JPL directions
	JPL.AsteroidDirection AS ad
	-- The observatory
	INNER JOIN KS.Observatory AS o ON o.ObservatoryID = ad.ObservatoryID
WHERE
	-- Only the selected asteroids
    ad.AsteroidID BETWEEN n0 AND (n1-1)
    -- Only the selected observatory
    AND o.ObservatoryShortName = ObservatoryShortName
ORDER BY ad.AsteroidID, ad.ObservatoryID, ad.TimeID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.AstrometricDirectionTest()
COMMENT "Get asteroid directions and positions to test MSE calculations; vectors from JPL."

BEGIN 

SELECT
	-- Key fields
	ad.AsteroidID,
	ad.TimeID,
	-- Time light leaves the asteroid
	ad.mjd AS tAst,
	-- The Earth geocenter position
	hve.qx AS qObs_x,
	hve.qy AS qObs_y,
	hve.qz AS qObs_z,
	-- The asteroid position
	hv.qx AS qAst_x,
	hv.qy AS qAst_y,
	hv.qz AS qAst_z,
	-- The asteroid velocity
	hv.vx AS vAst_x,
	hv.vy AS vAst_y,
	hv.vz AS vAst_z,
	-- Light time
	ad.LightTime,
	ad.mjd AS tObs,
	-- Astrometric direction from JPL RA/DEC
	ad.ux_ast AS ux,
	ad.uy_ast AS uy,
	ad.uz_ast AS uz,
	-- Distance
	ad.delta
FROM
	-- The JPL asteroid directions
	JPL.AsteroidDirection AS ad
	-- The observatory
	INNER JOIN KS.Observatory AS o ON o.ObservatoryID=ad.ObservatoryID
	-- Position of Earth
	INNER JOIN JPL.HorizonsBody AS hbe ON hbe.HorizonsBodyName='LB.Earth'
	INNER JOIN JPL.HorizonsVectors AS hve ON 
		hve.HorizonsBodyID=hbe.HorizonsBodyID AND
		hve.TimeID=ad.TimeID
	-- Position of Asteroids
	INNER JOIN JPL.HorizonsBody AS hb ON hb.SmallBodyID = ad.AsteroidID
	INNER JOIN JPL.HorizonsVectors AS hv ON 
		hv.HorizonsBodyID=hb.HorizonsBodyID AND
		hv.TimeID=ad.TimeID
WHERE 
	-- Only test the date range of MSE integrated orbits
	ad.TimeID BETWEEN 48000*24*60 AND 63000*24*60
	-- Only return Earth geocenter directions
	AND o.ObservatoryShortName='Geocenter'
ORDER BY ad.AsteroidID, ad.TimeID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.AstrometricDirectionAndVectorTest()
COMMENT "Get asteroid directions and positions to test MSE calculations; vectors from MSE not JPL."

BEGIN 

SELECT
	-- Key fields
	ad.AsteroidID,
	ad.TimeID,
	-- Time light leaves the asteroid
	ad.mjd AS tAst,
	-- The Earth geocenter position
	sve.qx AS qObs_x,
	sve.qy AS qObs_y,
	sve.qz AS qObs_z,
	-- The asteroid position
	av.qx AS qAst_x,
	av.qy AS qAst_y,
	av.qz AS qAst_z,
	-- The asteroid velocity
	av.vx AS vAst_x,
	av.vy AS vAst_y,
	av.vz AS vAst_z,
	-- Light time
	ad.LightTime,
	ad.mjd AS tObs,
	-- Astrometric direction from JPL RA/DEC
	ad.ux_ast AS ux,
	ad.uy_ast AS uy,
	ad.uz_ast AS uz,
	-- Distance
	ad.delta
FROM
	-- The JPL asteroid directions
	JPL.AsteroidDirection AS ad
	-- The observatory
	INNER JOIN KS.Observatory AS o ON o.ObservatoryID=ad.ObservatoryID
	-- The MSE state vectors for Earth
	INNER JOIN KS.StateVectors_Earth AS sve ON
		sve.TimeID = ad.TimeID
	-- The MSE state vectors for this asteroid
	INNER JOIN KS.AsteroidVectors AS av ON 
		av.AsteroidID = ad.AsteroidID AND
		av.TimeID = ad.TimeID
WHERE
	-- Only test the date range of MSE integrated orbits
	ad.TimeID BETWEEN 48000*24*60 AND 63000*24*60
	-- Only return Earth geocenter directions
	AND o.ObservatoryShortName='Geocenter'
ORDER BY ad.AsteroidID, ad.TimeID;

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
	ad.mjd - ad.LightTime / 1440.0 AS tAst,
	-- The Earth geocenter position
	hve.qx AS qObs_x,
	hve.qy AS qObs_y,
	hve.qz AS qObs_z,
	-- The asteroid position
	hv.qx AS qAst_x,
	hv.qy AS qAst_y,
	hv.qz AS qAst_z,
	-- The asteroid velocity
	hv.vx AS vAst_x,
	hv.vy AS vAst_y,
	hv.vz AS vAst_z,
	-- Light time
	ad.LightTime,
	ad.mjd AS tObs,
	-- Astrometric direction from JPL RA/DEC
	ad.ux_ast,
	ad.uy_ast,
	ad.uz_ast,
	-- Apparent direction from JPL RA/DEC
	ad.ux_app,
	ad.uy_app,
	ad.uz_app
FROM
	JPL.AsteroidDirection AS ad
	-- Position of Earth
	INNER JOIN JPL.HorizonsBody AS hbe ON hbe.HorizonsBodyName='LB.Earth'
	INNER JOIN JPL.HorizonsVectors AS hve ON 
		hve.HorizonsBodyID=hbe.HorizonsBodyID AND
		hve.TimeID=ad.TimeID
	-- Position of Asteroids
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
