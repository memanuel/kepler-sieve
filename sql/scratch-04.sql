SET @TimeID=59000*24*60;

SELECT
	ast.AsteroidID,
	dt0.MJD AS epoch
FROM
	-- Start with JPL reference elements
	JPL.AsteroidElement AS elt_jpl
	-- The matching asteroid and body to this quote from JPL
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.AsteroidNumber
	INNER JOIN KS.Body AS b ON b.BodyID = ast.BodyID
	-- The epoch when the elements are quoted
	INNER JOIN KS.DailyTime AS dt0 ON dt0.MJD = elt_jpl.epoch
	-- The epoch as of which we want the results
	INNER JOIN KS.DailyTime AS dt1 ON dt1.TimeID = @TimeID
	-- Try to join AsteroidElement_Ref on the date we want, dt1
	LEFT JOIN KS.AsteroidElement_Ref AS elt_ks ON
		elt_ks.AsteroidID = ast.AsteroidID AND
		elt_ks.TimeID = dt1.TimeID
WHERE
	-- Only take rows that are not already in KS.AsteroidElement_Ref
	elt_ks.AsteroidID IS NULL AND
	-- Only synchronize elements quoted PRIOR to the desired epoch
	dt0.TimeID < dt1.TimeID
ORDER BY ast.AsteroidID;