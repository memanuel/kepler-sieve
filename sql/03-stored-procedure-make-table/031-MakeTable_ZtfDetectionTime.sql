DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.MakeTable_DetectionTime()
COMMENT "Populate the ZTF.DetectionTime table from ZTF.Detection (all distinct times with detections)"
BEGIN 

INSERT INTO ZTF.DetectionTime
(MJD)
SELECT
    det.mjd AS MJD
FROM
    ZTF.Detection AS det
    -- Corresponding time already on ZTF.DetectionTime table
    LEFT JOIN ZTF.DetectionTime AS dt ON dt.MJD = det.mjd
-- Only detection times not already present
WHERE dt.DetectionTimeID IS NULL
-- We want distinct detection times only
GROUP BY det.mjd
ORDER BY det.mjd;

END
$$

DELIMITER ;
