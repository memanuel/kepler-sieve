SELECT
	ast.AsteroidID,
	it.TimeID,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node AS Omega,
	elt.omega_peri AS omega,
	elt.M AS M
FROM
	JPL.AsteroidElement AS elt
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt.AsteroidNumber
	INNER JOIN KS.IntegrationTime AS it ON it.MJD = elt.epoch
ORDER BY ast.AsteroidID;