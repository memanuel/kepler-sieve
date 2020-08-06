CREATE OR REPLACE 
DEFINER = kepler
VIEW KS.BodyCollectionDisplay
AS
SELECT
	bc.BodyCollectionName,
	b.BodyName,
	bce.BodyNumber,
	mb.M * 1.0E9 AS SolarMass_Em9,
	ROW_NUMBER() OVER (ORDER BY bc.BodyCollectionID, b.SortOrder) AS SortOrder
FROM
	KS.BodyCollectionEntry AS bce
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionID = bce.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bce.BodyID
	LEFT JOIN JPL.HorizonsBody AS hb ON hb.BodyID = b.BodyID
	LEFT JOIN JPL.MassiveBody AS mb ON mb.HorizonsBodyID = hb.HorizonsBodyID;
