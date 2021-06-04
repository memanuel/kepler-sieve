SELECT * FROM KS.AsteroidDirections2 WHERE AsteroidID=1;

INSERT INTO KS.AsteroidDirections2
(AsteroidID, TimeID, tObs, ux, uy, uz, LightTime)
SELECT
	ad.AsteroidID, 
	ad.TimeID, 
	ad.mjd AS tObs, 
	ad.ux, 
	ad.uy, 
	ad.uz,
	ad.LightTime
FROM
	KS.AsteroidDirections AS ad
WHERE ad.AsteroidID <= 1000;