-- truncate TABLE KS.Tracklet;
-- CALL KS.MakeTable_Tracklet(7962056, 8923503);

CALL KS.MakeTable_Tracklet(7962056+0, 7962056+875);

CALL KS.MakeTable_Tracklet(7962056+875, 7962056+1000);

SELECT count(*) AS batch_size_times FROM KS.DetectionTimePair_insert;
SELECT count(*) AS batch_size_rows FROM KS.TrackletBatch;

SELECT * FROM KS.TrackletBatch;

SELECT * FROM KS.Tracklet LIMIT 100;
SELECT max(DetectionTimePairID) AS maxDetectionTimePairID FROM KS.Tracklet;


SELECT * FROM KS.DetectionTimePair_present;
SELECT * FROM KS.DetectionTimePair_insert;

INSERT INTO KS.DetectionTimePair_insert
(DetectionTimePairID)
SELECT
	dtp.DetectionTimePairID
FROM
	KS.DetectionTimePair AS dtp
WHERE	
	-- Only tracklets in the delected range of time pairs
	dtp.DetectionTimePairID BETWEEN DetectionTimePairID_0 AND (DetectionTimePairID_1-1)
	-- Only those not already present
	AND NOT EXISTS (
	SELECT dtpp.DetectionTimePairID
	FROM KS.DetectionTimePair_present AS dtpp
	WHERE dtpp.DetectionTimePairID = dtp.DetectionTimePairID);