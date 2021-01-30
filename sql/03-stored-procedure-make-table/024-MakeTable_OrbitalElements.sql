DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_OrbitalElements(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Populate the OrbitalElements table from OrbitalElements_Planets"
BEGIN 

-- TimeID for this date range
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

-- Select from StateVectors_Planets into StateVectors
REPLACE INTO KS.OrbitalElements
(TimeID, BodyID, MJD, a, e, inc, Omega_node, omega_peri, f, M)
SELECT
    elt.TimeID,
    elt.BodyID,
    elt.MJD,
    elt.a,
    elt.e,
    elt.inc,
    elt.Omega_node,
    elt.omega_peri,
    elt.f,
    elt.M
FROM 
    KS.OrbitalElements_Planets as elt
WHERE
	@TimeID_0 <= elt.TimeID AND elt.TimeID < @TimeID_1;
   
END
$$

DELIMITER ;
