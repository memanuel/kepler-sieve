DELIMITER $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_Tracklet(
	IN DetectionTimeID_0 INT,
	IN DetectionTimeID_1 INT
)
COMMENT "Get the number of rows missing from KS.RawDetection"
BEGIN 

-- Temporary table for inserting into Tracklet
CREATE OR REPLACE TEMPORARY TABLE KS.TrackletBatch LIKE KS.Tracklet;
ALTER TABLE KS.TrackletBatch 
	ADD COLUMN r DOUBLE NOT NULL DEFAULT 0.0,
	ADD COLUMN x0 DOUBLE NOT NULL,
	ADD COLUMN x1 DOUBLE NOT NULL,
	ADD COLUMN y0 DOUBLE NOT NULL,
	ADD COLUMN y1 DOUBLE NOT NULL,
	ADD COLUMN z0 DOUBLE NOT NULL,
	ADD COLUMN z1 DOUBLE NOT NULL;
	
-- Batch of candidate tracklets
INSERT INTO KS.TrackletBatch
(DetectionID_1, DetectionID_2, SkyPatchID, mjd, ux, uy, uz, ux_dot, uy_dot, uz_dot, mag, x0, x1, y0, y1, z0, z1)
SELECT
	-- The two DetectionIDs	
	d1.DetectionID AS DetectionID_1,
	d2.DetectionID AS DetectionID_2,
	-- The SkyPatchID of the mean
	sp.SkyPatchID,
	-- The mean date
	dtp.mjd,
	-- Average direction
	(d1.ux + d2.ux)/2 AS ux,
	(d1.uy + d2.uy)/2 AS uy,
	(d1.uz + d2.uz)/2 AS uz,
	-- Derivative of direction
	(d2.ux - d1.ux)/dtp.dt AS ux_dot,
	(d2.uy - d1.uy)/dtp.dt AS uy_dot,
	(d2.uz - d1.uz)/dtp.dt AS uz_dot,
	-- Average magnitude
	(d1.mag + d2.mag)/2 AS mag,
	-- Bounds on x from SkyPatch
	LEAST(sp.x00, sp.x01, sp.x10, sp.x11) AS x0,
	GREATEST(sp.x00, sp.x01, sp.x10, sp.x11) AS x1,
	LEAST(sp.y00, sp.y01, sp.y10, sp.y11) AS y0,
	GREATEST(sp.y00, sp.y01, sp.y10, sp.y11) AS y1,
	LEAST(sp.z00, sp.z01, sp.z10, sp.z11) AS z0,
	GREATEST(sp.z00, sp.z01, sp.z10, sp.z11) AS z1
FROM
	-- Start with pairs of detection times
	KS.DetectionTimePair AS dtp 
	-- First detection from this pair of times
	INNER JOIN KS.Detection AS d1 ON 
		d1.DetectionTimeID = dtp.DetectionTimeID_1
	INNER JOIN KS.Detection AS d2 ON 
		d2.DetectionTimeID = dtp.DetectionTimeID_2 AND 
		d2.SkyPatchID = d1.SkyPatchID
	-- Possible values of the SkyPatch
	INNER JOIN KS.SkyPatch AS sp ON sp.SkyPatchID IN (d1.SkyPatchID, d2.SkyPatchID)
WHERE
	(DetectionTimeID_0 <= dtp.DetectionTimeID_1) AND (dtp.DetectionTimeID_1 < DetectionTimeID_1);

-- Set the distance r
UPDATE KS.TrackletBatch
SET 
	r = SQRT(POW(ux,2)+POW(uy,2)+POW(uz,2));

-- Insert this batch into the main table
INSERT IGNORE INTO KS.Tracklet
(DetectionID_1, DetectionID_2, SkyPatchID, mjd, ux, uy, uz, ux_dot, uy_dot, uz_dot, mag)
SELECT
	-- The pair of detections
	tb.DetectionID_1,
	tb.DetectionID_2,
	-- The location in the sky
	tb.SkyPatchID,
	-- The mean date
	tb.mjd,
	-- The direction
	tb.ux/r AS ux,
	tb.uy/r AS uy,
	tb.uz/r AS uz,
	-- Time derivative of the location
	tb.ux_dot,
	tb.uy_dot,
	tb.uz_dot,
	-- Magnitude
	tb.mag
FROM
	-- The pair of directions for the batch process
	KS.TrackletBatch AS tb
WHERE
	-- Only when the mean direction matches the SkyPatch
	(tb.ux BETWEEN tb.x0*tb.r AND tb.x1*tb.r) AND
	(tb.uy BETWEEN tb.y0*tb.r AND tb.y1*tb.r) AND
	(tb.uz BETWEEN tb.z0*tb.r AND tb.z1*tb.r);

-- Clean up temp table
DROP TEMPORARY TABLE KS.TrackletBatch;

END $$

DELIMITER ;

