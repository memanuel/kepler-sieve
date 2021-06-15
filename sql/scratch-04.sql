SELECT * FROM KS.DetectionTime;

UPDATE KS.DetectionTime 
SET HiResTimeID = CAST(FLOOR(mjd*1440) AS INT);

ALTER TABLE KS.DetectionTime
ADD KEY IDX_DetectionTime_HiResTimeID(HiResTimeID);
