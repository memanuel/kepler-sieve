SELECT
	dna.DetectionID,
	dna.AsteroidID,
	dna.s*3600*60 AS sec,
	dna.LightTime
FROM
	KS.DetectionNearAsteroid AS dna
WHERE dna.AsteroidID < 100 AND dna.s*3600*60 < 50
ORDER BY dna.s;