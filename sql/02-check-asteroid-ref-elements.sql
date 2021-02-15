-- Count number of asteroids that have references pushed forward 
SELECT
	COUNT(elt_ks.AsteroidID) AS ast_count
FROM
	JPL.AsteroidElement AS elt_jpl
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.asteroidNumber
	INNER JOIN KS.DailyTime AS dt ON dt.MJD = 59000
	INNER JOIN KS.AsteroidElement_Ref AS elt_ks ON
		elt_ks.AsteroidID = ast.AsteroidID AND
		elt_ks.TimeID = dt.TimeID
WHERE elt_jpl.epoch < elt_ks.epoch;


-- Compare the elements to see if they are close
SELECT
	elt_jpl.AsteroidNumber,
	elt_jpl.AsteroidName,
	elt_jpl.epoch,
	elt_jpl.a AS a_jpl,
	elt_ks.a AS a_ks,	
	elt_ks.a - elt_jpl.a AS da,
	elt_ks.e - elt_jpl.e AS de,
	elt_ks.inc - elt_jpl.inc AS dinc,
	MOD(elt_ks.Omega_node - elt_jpl.Omega_node + 2.0*PI(), 2.0*PI()) AS dOmega,
	MOD(elt_ks.omega_peri - elt_jpl.omega_peri + 2.0*PI(), 2.0*PI()) AS domega,
	MOD(elt_ks.M - elt_jpl.M + 2.0*PI(), 2.0*PI()) AS dM
FROM
	JPL.AsteroidElement AS elt_jpl
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber = elt_jpl.asteroidNumber
	INNER JOIN KS.DailyTime AS dt ON dt.MJD = 59000
	INNER JOIN KS.AsteroidElement_Ref AS elt_ks ON
		elt_ks.AsteroidID = ast.AsteroidID AND
		elt_ks.TimeID = dt.TimeID
WHERE elt_jpl.epoch < elt_ks.epoch		
ORDER BY elt_jpl.AsteroidNumber
LIMIT 100;

