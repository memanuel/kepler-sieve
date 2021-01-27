SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeID = b.BodyTypeID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'Planets'
WHERE bt.IsLargeBody_JPL;