SELECT
	COUNT(dna.DetectionID) AS Hits
FROM
	KS.DetectionNearAsteroid_v1 AS dna
WHERE
	dna.AsteroidID BETWEEN 0 AND 1000;

SELECT
	dna.AsteroidID,
	COUNT(dna.DetectionID) AS Hits
FROM
	KS.DetectionNearAsteroid_v1 AS dna
WHERE
	dna.AsteroidID BETWEEN 0 AND 1000
GROUP BY dna.AsteroidID;

SELECT
	dna.AsteroidID,
	dna.DetectionID,
	d.SkyPatchID,
	d.TimeID
FROM
	KS.DetectionNearAsteroid_v1 AS dna
	INNER JOIN KS.Detection AS d ON d.DetectionID = dna.DetectionID
WHERE
	dna.AsteroidID = 3
ORDER BY DetectionID;

SELECT
	asp.AsteroidID,
	asp.Segment,
	asp.SkyPatchID,
	asp.TimeID_0,
	asp.TimeID_1
FROM
	KS.AsteroidSkyPatch AS asp
WHERE
	asp.AsteroidID=3 AND asp.SkyPatchID=14173326;