SELECT
	COUNT(dna.DetectionID) AS Hits
FROM
	KS.DetectionNearAsteroidCandidate_v1 AS dna
WHERE
	dna.AsteroidID BETWEEN 0 AND 1000;

SELECT
	dna.AsteroidID,
	COUNT(dna.DetectionID) AS Hits
FROM
	KS.DetectionNearAsteroidCandidate_v1 AS dna
WHERE
	dna.AsteroidID BETWEEN 0 AND 1000
GROUP BY dna.AsteroidID;

SELECT
	dna.AsteroidID,
	dna.DetectionID,
	d.SkyPatchID
FROM
	KS.DetectionNearAsteroidCandidate_v1 AS dna
	INNER JOIN KS.Detection AS d ON d.DetectionId = dna.DetectionID
WHERE
	dna.AsteroidID = 86
ORDER BY DetectionID;