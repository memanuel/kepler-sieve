DELIMITER $$

-- *********************************************************************************
-- Helper function: count the number of missing rows
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_RawDetection_ZTF_RowCount()
COMMENT "Get the number of rows missing from KS.RawDetection"
BEGIN 

SELECT tbl.TABLE_ROWS 
INTO @rc_tot
FROM INFORMATION_SCHEMA.TABLES AS tbl
WHERE tbl.TABLE_SCHEMA = 'ZTF' AND tbl.TABLE_NAME='Detection';

SELECT tbl.TABLE_ROWS
INTO @rc_exist
FROM INFORMATION_SCHEMA.TABLES AS tbl
WHERE tbl.TABLE_SCHEMA = 'KS' AND tbl.TABLE_NAME='RawDetection';

SELECT
(@rc_tot - @rc_exist) AS RowCount;

END
$$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_RawDetection_ZTF(
	IN sz INT)
COMMENT "Populate the KS.RawDetection table from ZTF.Detection (all distinct times with detections)"
BEGIN 

CREATE OR REPLACE TEMPORARY TABLE KS.RawDetectionInsert LIKE KS.RawDetection;

INSERT INTO KS.RawDetectionInsert
(DataSourceID, ObservatoryID, DetectionTimeID, SourceDetectionID, RA, `DEC`)
SELECT
	ds.DataSourceID,
	obs.ObservatoryID,
	dt.DetectionTimeID,
	det.DetectionID AS SourceDetectionID,
	det.RA,
	det.`DEC`
FROM
	-- Start with all ZTF detections
	ZTF.Detection AS det
	-- Get the data source and observatory
	INNER JOIN KS.DataSource AS ds ON ds.DataSourceCD = 'ZTF'
	INNER JOIN KS.Observatory AS obs ON obs.ObservatoryID = ds.ObservatoryID
	-- Get the DetectionTime by joining on MJD
	INNER JOIN KS.DetectionTime AS dt ON dt.MJD = det.MJD
WHERE
	-- Only take rows that aren't already in KS.RawDetection
	NOT EXISTS (
	SELECT 
	rd.DetectionID
	FROM KS.RawDetection AS rd
	WHERE rd.DataSourceID = ds.DataSourceID AND
	rd.SourceDetectionID = det.DetectionID
	)
ORDER BY det.DetectionID
LIMIT sz;

INSERT INTO KS.RawDetection
(DataSourceID, ObservatoryID, DetectionTimeID, SourceDetectionID, RA, `DEC`)
SELECT
DataSourceID, ObservatoryID, DetectionTimeID, SourceDetectionID, RA, `DEC`
FROM KS.RawDetectionInsert
ORDER BY DetectionID;

DROP TEMPORARY TABLE KS.RawDetectionInsert;

END
$$

DELIMITER ;
