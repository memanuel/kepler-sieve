SELECT
	elt_jpl.AsteroidNumber,
	elt_jpl.AsteroidName,
	elt_jpl.a AS a_jpl,
	elt_ks.a AS a_ks,	
	elt_ks.a - elt_jpl.a AS a_err,
	elt_ks.e - elt_jpl.e AS e_err,
	elt_ks.inc - elt_jpl.inc AS inc_err,
	elt_ks.Omega_node - elt_jpl.Omega_node AS Omega_err,
	elt_ks.omega_peri - elt_jpl.omega_peri AS omega_err,
	elt_ks.M - elt_jpl.M AS M_err,
	elt_ks.f - elt_jpl.f AS f_err
FROM
	JPL.AsteroidElement AS elt_jpl
	LEFT JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.asteroidNumber
	INNER JOIN KS.DailyTime AS dt ON dt.MJD = elt_jpl.epoch
	LEFT JOIN KS.AsteroidElement_Ref AS elt_ks ON
		elt_ks.AsteroidID = ast.AsteroidID AND
		elt_ks.TimeID = dt.TimeID
ORDER BY elt_jpl.AsteroidNumber
LIMIT 100;

WITH t1 AS(
SELECT
	elt_jpl.AsteroidNumber,
	elt_jpl.AsteroidName,
	elt_jpl.a AS a_jpl,
	elt_ks.a AS a_ks,	
	elt_ks.a - elt_jpl.a AS a_err,
	elt_ks.e - elt_jpl.e AS e_err,
	elt_ks.inc - elt_jpl.inc AS inc_err,
	elt_ks.Omega_node - elt_jpl.Omega_node AS Omega_node_err,
	elt_ks.omega_peri - elt_jpl.omega_peri AS omega_peri_err,
	elt_ks.M - elt_jpl.M AS M_err,
	elt_ks.f - elt_jpl.f AS f_err
FROM
	JPL.AsteroidElement AS elt_jpl
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.asteroidNumber
	INNER JOIN KS.DailyTime AS dt ON dt.MJD = elt_jpl.epoch
	INNER JOIN KS.AsteroidElement_Ref AS elt_ks ON
		elt_ks.AsteroidID = ast.AsteroidID AND
		elt_ks.TimeID = dt.TimeID
)
SELECT
	sqrt(AVG(t1.a_err*t1.a_err)) AS a_err
FROM t1;