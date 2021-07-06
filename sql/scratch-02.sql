SELECT 
	dna.AsteroidID,
	dna.DetectionID,
	dna.s*(180/pi())*3600.0 AS s_sec,
	dna.LightTime
FROM 
	KS.DetectionNearAsteroid AS dna 
WHERE 
	dna.AsteroidID=3;