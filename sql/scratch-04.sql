SELECT
	CAST(FLOOR(MIN(dt.mjd)) AS INT) AS mjd0,
	CAST(CEILING(Max(dt.mjd)) AS INT) AS mjd1
FROM
	KS.DetectionTime AS dt
WHERE NOT EXISTS (
	SELECT d.DetectionID
	FROM KS.Detection AS d
	WHERE d.DetectionID = dt.DetectionTimeID
    );
    
CALL KS.GetRawDetectionDates();