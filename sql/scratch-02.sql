DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_StateVectors(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Populate the StateVectors table from StateVectors_Planets"
BEGIN 

-- TimeID for this date range
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

-- Select from StateVectors_Planets into StateVectors
REPLACE INTO KS.StateVectors
(TimeID, BodyID, MJD, qx, qy, qz, vx, vy, vz)
SELECT
    sv.TimeID,
    sv.BodyID,
    sv.MJD,
    sv.qx,
    sv.qy,
    sv.qz,
    sv.vx,
    sv.vy,
    sv.vz
FROM 
    KS.StateVectors_Planets as sv
WHERE
    mjd0 <= sv.TimeID AND sv.TimeID < mjd1;
   
-- Add the Earth-Moon Barycenter
CALL KS.Calc_EarthMoonBarycenter()

END
$$

DELIMITER ;
