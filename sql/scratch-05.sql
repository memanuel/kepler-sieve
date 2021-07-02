SELECT
	min(dna.DetectionID) AS MinDetectionID
FROM
	KS.DetectionNearAsteroid_v1 dna;
	
SELECT
	dna.DetectionID,
	dna.AsteroidID
FROM
	KS.DetectionNearAsteroid_v1 dna
WHERE dna.DetectionId=3;