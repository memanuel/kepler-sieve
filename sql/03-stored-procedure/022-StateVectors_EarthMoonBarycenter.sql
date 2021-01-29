DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.Calc_EarthMoonBarycenter()
COMMENT "Populate the StateVectors table from StateVectors_Planets"
BEGIN 

-- Create temporary table with EMB according to Planets integration
CREATE TEMPORARY TABLE emb_p
AS
SELECT
 	it.TimeID,
 	b_emb.BodyID,
 	it.MJD,
	SUM(bw.Weight * sv.qx) AS qx,
	SUM(bw.Weight * sv.qy) AS qy,
	SUM(bw.Weight * sv.qz) AS qz,
	SUM(bw.Weight * sv.vx) AS vx,
	SUM(bw.Weight * sv.vy) AS vy,
	SUM(bw.Weight * sv.vz) AS vz
FROM 
	KS.Body AS b_emb	
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = 'PS3'
	INNER JOIN KS.BarycenterWeight AS bw ON
		bw.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bw.BodyID
	INNER JOIN KS.IntegrationTime AS it
	INNER JOIN KS.StateVectors_Planets AS sv ON
		sv.TimeID = it.TimeID AND
		sv.BodyID = b.BodyID
WHERE
	b_emb.BodyName = 'Earth-Moon Barycenter'
GROUP BY it.TimeID;

-- Insert EMB into StateVectors_Planets
REPLACE INTO KS.StateVectors_Planets
(TimeID, BodyID, MJD, qx, qy, qz, vx, vy, vz)
SELECT
	emb.TimeID,
	emb.BodyID,
    emb.MJD,
    emb.qx,
    emb.qy,
    emb.qz,
    emb.vx,
    emb.vy,
    emb.vz
FROM emb_p as emb;

-- Create temporary table with EMB according to Planets integration
CREATE TEMPORARY TABLE emb_d
AS
SELECT
 	it.TimeID,
 	b_emb.BodyID,
 	it.MJD,
	SUM(bw.Weight * sv.qx) AS qx,
	SUM(bw.Weight * sv.qy) AS qy,
	SUM(bw.Weight * sv.qz) AS qz,
	SUM(bw.Weight * sv.vx) AS vx,
	SUM(bw.Weight * sv.vy) AS vy,
	SUM(bw.Weight * sv.vz) AS vz
FROM 
	KS.Body AS b_emb	
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = 'PS3'
	INNER JOIN KS.BarycenterWeight AS bw ON
		bw.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bw.BodyID
	INNER JOIN KS.IntegrationTime AS it
	INNER JOIN KS.StateVectors_DE435 AS sv ON
		sv.TimeID = it.TimeID AND
		sv.BodyID = b.BodyID
WHERE
	b_emb.BodyName = 'Earth-Moon Barycenter'
GROUP BY it.TimeID;

-- Insert EMB into StateVectors_DE435
REPLACE INTO KS.StateVectors_DE435
(TimeID, BodyID, MJD, qx, qy, qz, vx, vy, vz)
SELECT
	emb.TimeID,
	emb.BodyID,
    emb.MJD,
    emb.qx,
    emb.qy,
    emb.qz,
    emb.vx,
    emb.vy,
    emb.vz
FROM emb_d as emb;

-- Clean up
DROP TEMPORARY TABLE emb_p
DROP TEMPORARY TABLE emb_d

END
$$

DELIMITER ;
