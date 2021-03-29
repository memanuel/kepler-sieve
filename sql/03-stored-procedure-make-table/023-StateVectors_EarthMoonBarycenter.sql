DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.Calc_EarthMoonBarycenter(
    IN mjd0 INT,
    IN mjd1 INT
)
COMMENT "Populate the StateVectors table from StateVectors_Planets"
BEGIN 

-- TimeID for this date range
SET @TimeID_0 = mjd0 * 24 * 60;
SET @TimeID_1 = mjd1 * 24 * 60;

-- Insert EMB into StateVectors_Planets
REPLACE INTO KS.StateVectors
(TimeID, BodyID, mjd, qx, qy, qz, vx, vy, vz)
SELECT
	sv.TimeID,
	sv.BodyID,
	sv.mjd,
	SUM(bw.Weight*sv.qx) AS qx,
	SUM(bw.Weight*sv.qy) AS qy,
	SUM(bw.Weight*sv.qz) AS qz,
	SUM(bw.Weight*sv.vx) AS vx,
	SUM(bw.Weight*sv.vy) AS vy,
	SUM(bw.Weight*sv.vz) AS vz
FROM
	KS.StateVectors_Planets AS sv
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = 'PS3'
	INNER JOIN KS.BarycenterWeight AS bw ON
		bw.BodyCollectionID = bc.BodyCollectionID AND
		bw.BodyID = sv.BodyID
	INNER JOIN KS.Body AS b_emb ON b_emb.BodyName = 'Earth-Moon Barycenter'
WHERE
    mjd0 <= sv.TimeID AND sv.TimeID < mjd1
GROUP BY sv.TimeID;

END
$$

DELIMITER ;
