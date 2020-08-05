CREATE OR REPLACE 
DEFINER = kepler
VIEW KS.BodyCollectionDisplay
AS
SELECT
	bc.BodyCollectionName,
	b.BodyName,
	bce.BodyNumber,
	ROW_NUMBER() OVER (ORDER BY bc.BodyCollectionID, bce.BodyNumber) AS SortOrder	
FROM
	KS.BodyCollectionEntry AS bce
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionID = bce.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bce.BodyID;
