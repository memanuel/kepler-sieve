SELECT * FROM KS.StateVectors_Planets LIMIT 20;
SELECT * FROM KS.StateVectors_Planets WHERE TimeID = 59000*24*60;
SELECT * FROM KS.StateVectors_Planets WHERE TimeID > 59000*24*60 GROUP BY TimeID LIMIT 100 ;

SELECT * FROM KS.AsteroidVectors WHERE AsteroidID = 1 AND TimeID > 59000*24*60;

SELECT * FROM KS.IntegrationTime;
SELECT * FROM KS.Minutes;


SELECT * FROM KS.HiResTime;

SELECT 
	it1.TimeID + m._ AS TimeID,
	it1.mjd + m.mjd_offset AS mjd,
	it1.jd + m.mjd_offset AS JD
FROM 
	KS.IntegrationTime AS it1
	INNER JOIN KS.IntegrationTime AS it2 ON it2.TimeID = it1.TimeID+5
	CROSS JOIN KS.Minutes AS m ON m._ < 5