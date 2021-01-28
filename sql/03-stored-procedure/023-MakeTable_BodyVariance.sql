DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_BodyVariance()
COMMENT "Populate the BodyVariance table from HorizonsVectors"
BEGIN 


INSERT IGNORE INTO JPL.BodyVariance 
(BodyID, var_q, var_v, sd_q, sd_v)
WITH t1 as(
    SELECT
        hv.HorizonsBodyID,
        variance(hv.qx) + variance(hv.qy) + variance(hv.qz) AS var_q,
        variance(hv.vx) + variance(hv.vy) + variance(hv.vz) AS var_v
    FROM
        JPL.HorizonsVectors AS hv
        INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = hv.HorizonsBodyID
    GROUP BY hv.HorizonsBodyID
)
SELECT
	hb.BodyID,
	t1.var_q,
	t1.var_v,
	sqrt(t1.var_q) AS sd_q,
	sqrt(t1.var_v) AS sd_v
FROM 
	t1
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = t1.HorizonsBodyID;

END
$$

DELIMITER ;
