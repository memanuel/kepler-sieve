SELECT
	ast.AsteroidID,
	ast.AsteroidName,
	elt_jpl.a AS a_jpl,
	elt_mse.a AS a_mse,
	(elt_mse.a - elt_jpl.a) AS delta_a
FROM
	JPL.AsteroidElement AS elt_jpl
	INNER JOIN KS.DailyTime AS dt ON dt.MJD=elt_jpl.epoch
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.AsteroidNumber
	INNER JOIN KS.AsteroidElements AS elt_mse ON
		elt_mse.TimeID = dt.TimeID AND
		elt_mse.AsteroidID = ast.AsteroidID
WHERE elt_jpl.AsteroidNumber<=10;
