DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.GetDetectionTimes()
COMMENT "Get all the distinct detection times ZTF.Detection table by querying ZTF.DetectionTime."

BEGIN 

SELECT
	dt.DetectionTimeID,
	CAST(FLOOR(dt.mjd*1440) AS INT) AS HiResTimeID,
	dt.mjd,
	ds.DataSourceID,	
	obs.ObservatoryID
FROM
	-- Start with ZTF detection times
	ZTF.DetectionTime AS dt
	-- The DataSource
	INNER JOIN KS.DataSource AS ds ON ds.DataSourceCD = 'ZTF'
	-- The Observatory
	INNER JOIN KS.Observatory AS obs ON obs.ObservatoryShortName = 'ZTF'
ORDER BY dt.DetectionTimeID;

END $$

DELIMITER ;
