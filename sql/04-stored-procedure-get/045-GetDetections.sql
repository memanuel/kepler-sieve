DELIMITER $$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetections(
    IN d0 INT,
    IN d1 INT
)
COMMENT "Get asteroid detections in a range of detection IDs."

BEGIN 

SELECT
	-- Integer ID fields
	d.DetectionID,
	d.SkyPatchID,
	d.TimeID,
	-- Data payload
	d.mjd,
	d.ux,
	d.uy,
	d.uz
FROM
	KS.Detection AS d
WHERE
	d.DetectionID BETWEEN d0 AND (d1-1)
ORDER BY d.DetectionID;

END
$$

-- ********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDetectionsObs(
    IN d0 INT,
    IN d1 INT
)
COMMENT "Get asteroid detections in a range of detection IDs; also return the observatory position."

BEGIN 

SELECT
	-- Integer ID fields
	d.DetectionID,
	d.SkyPatchID,
	d.TimeID,
	d.mjd,
	-- Data payload
	d.ux,
	d.uy,
	d.uz,
	-- Location of observatory - from join to DetectionTime
	dt.qObs_x,
	dt.qObs_y,
	dt.qObs_z
FROM
	KS.Detection AS d
	INNER JOIN KS.DetectionTime AS dt ON dt.DetectionTimeID = d.DetectionTimeID
WHERE
	d.DetectionID BETWEEN d0 AND (d1-1)
ORDER BY d.DetectionID;

END
$$

DELIMITER ;
