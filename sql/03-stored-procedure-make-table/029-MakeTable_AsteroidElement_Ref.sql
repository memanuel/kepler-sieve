DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_AsteroidElement_Ref()
COMMENT "Populate the AsteroidElement_Ref table from JPL.AsteroidElement"
BEGIN 

REPLACE INTO KS.AsteroidElement_Ref
(AsteroidID, TimeID, epoch, a, e, inc, Omega_node, omega_peri, f, M, eccentric_anomaly, period, mean_motion)
SELECT
	ast.AsteroidID,
	it.TimeID,
	elt.epoch,
	elt.a,
	elt.e,
	elt.inc,
	elt.Omega_node,
	elt.omega_peri,
	elt.f,
	elt.M,
	elt.eccentric_anomaly,
	elt.period,
	elt.mean_motion
FROM
	JPL.AsteroidElement AS elt
	INNER JOIN KS.Asteroid AS ast ON ast.AsteroidNumber=elt.AsteroidNumber
	INNER JOIN KS.IntegrationTime AS it ON it.MJD = elt.epoch;

END
$$

DELIMITER ;
