SELECT 
	dna.AsteroidID,
	dna.DetectionID,
	dna.s*(180/pi())*3600.0 AS s_sec,
	dna1.s*(180/pi())*3600.0 AS s1_sec
FROM 
	KS.DetectionNearAsteroid AS dna
	LEFT JOIN KS.DetectionNearAsteroid_v1 AS dna1 ON
		dna1.DetectionID = dna.DetectionID AND
		dna1.AsteroidID = dna.AsteroidID
WHERE 
	dna.AsteroidID<=10
ORDER BY dna.AsteroidID, dna.DetectionID;