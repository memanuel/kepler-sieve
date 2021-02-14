DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_AsteroidElement_Ref(
    IN epoch_ INT
)
COMMENT "Populate the AsteroidElement_Ref table from JPL.AsteroidElement"
BEGIN 

REPLACE INTO KS.AsteroidElement_Ref
(AsteroidID, TimeID, epoch, a, e, inc, Omega_node, omega_peri, f, M, d, v, h, period, mean_motion, T_peri, pomega)
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
WHERE
    -- Only match the desired epoch
	elt.epoch = epoch_;

END
$$

DELIMITER ;
