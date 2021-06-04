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
-- 	ad.r, 
-- 	ad.rDot,
-- 	ad.delta,
-- 	ad.deltaDot,
--	ad.Mag,	
	-- Misc
	ad.LightTime,
FROM
	JPL.AsteroidDirection AS ad
	INNER JOIN JPL.HorizonsBody AS hb ON hb.SmallBodyID = ad.AsteroidID
	INNER JOIN JPL.HorizonsVectors AS hv ON 
		hv.HorizonsBodyID=hb.HorizonsBodyID AND
		hv.TimeID=ad.TimeID;
