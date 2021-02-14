SELECT
	ast.AsteroidID,
	dt.TimeID,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node,
	elt.omega_peri,
	elt.f,
	elt.M,
	0.0 AS d,
	0.0 AS v,
	0.0 AS h,
	elt.period,
	elt.mean_motion,
	0.0 AS T_peri,
	0.0 AS pomega
FROM
	JPL.AsteroidElement AS elt
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber=elt.AsteroidNumber
    -- The time as of which the elements were quoted
    INNER JOIN KS.DailyTime AS dt ON dt.MJD = elt.epoch
WHERE elt.epoch = 59000;    