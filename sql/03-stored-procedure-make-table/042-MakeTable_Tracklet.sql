DELIMITER $$

-- *********************************************************************************
CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_Tracklet(
	IN DetectionTimePairID_0 INT,
	IN DetectionTimePairID_1 INT
)
COMMENT "Get the number of rows missing from KS.RawDetection"
BEGIN 

-- Temporary table for detection time pairs already present in the tracklet table
CREATE OR REPLACE TEMPORARY TABLE KS.DetectionTimePair_present (
	DetectionTimePairID INT NOT NULL PRIMARY KEY
);

INSERT INTO KS.DetectionTimePair_present
(DetectionTimePairID)
SELECT
	t.DetectionTimePairID
FROM
	-- Start with tracklets
	KS.Tracklet AS t
WHERE	
	-- Only tracklets in the delected range of time pairs
	t.DetectionTimePairID BETWEEN DetectionTimePairID_0 AND (DetectionTimePairID_1-1)
GROUP BY t.DetectionTimePairID;

-- Temporary table for detection time pairs in the selected range we want to insert
CREATE OR REPLACE TEMPORARY TABLE KS.DetectionTimePair_insert (
	DetectionTimePairID INT NOT NULL PRIMARY KEY
);

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
	
-- Batch of tracklets to insert into KS.Tracklet later
INSERT INTO KS.TrackletBatch
(DetectionTimePairID, DetectionID_1, DetectionID_2, SkyPatchID, 
mjd_bar, mjd1, mjd2,
ux_bar, uy_bar, uz_bar, ux_dot, uy_dot, uz_dot, u_dot,
ux1, uy1, uz1, ux2, uy2, uz2, 
mag_bar, dmag, mag1, mag2,
x0, x1, y0, y1, z0, z1)
SELECT
	-- The time pair
	dtp.DetectionTimePairID,
	-- The two DetectionIDs	
	d1.DetectionID AS DetectionID_1,
	d2.DetectionID AS DetectionID_2,
	-- The SkyPatchID of the mean
	sp.SkyPatchID,
	-- The three dates
	dtp.mjd_bar,
	dtp.mjd1,
	dtp.mjd2,
	-- Average direction
	(d1.ux + d2.ux)/2 AS ux_bar,
	(d1.uy + d2.uy)/2 AS uy_bar,
	(d1.uz + d2.uz)/2 AS uz_bar,
	-- Derivative of direction
	(d2.ux - d1.ux)/dtp.dt AS ux_dot,
	(d2.uy - d1.uy)/dtp.dt AS uy_dot,
	(d2.uz - d1.uz)/dtp.dt AS uz_dot,
	0.0 AS u_dot,
	-- Direction 1
	d1.ux AS ux1,
	d1.uy AS uy1,
	d1.uz AS uz1,
	-- Direction 2
	d2.ux AS ux2,
	d2.uy AS uy2,
	d2.uz AS uz2,	
	-- Magnitude
	(d1.mag + d2.mag)/2 AS mag_bar,
	(d2.mag - d1.mag) AS dmag,
	d1.mag AS mag1,
	d2.mag AS mag2,
	-- Bounds on x, y and z from SkyPatch
	LEAST(sp.x00, sp.x01, sp.x10, sp.x11) AS x0,
	GREATEST(sp.x00, sp.x01, sp.x10, sp.x11) AS x1,
	LEAST(sp.y00, sp.y01, sp.y10, sp.y11) AS y0,
	GREATEST(sp.y00, sp.y01, sp.y10, sp.y11) AS y1,
	LEAST(sp.z00, sp.z01, sp.z10, sp.z11) AS z0,
	GREATEST(sp.z00, sp.z01, sp.z10, sp.z11) AS z1
FROM
	-- Start with the desired pairs of detection times
	KS.DetectionTimePair_insert AS dtpp 
	-- The DetectionTimePair we want to insert
	INNER JOIN KS.DetectionTimePair AS dtp ON dtp.DetectionTimePairID = dtpp.DetectionTimePairID
	-- First detection from this pair of times
	INNER JOIN KS.Detection AS d1 ON 
		d1.DetectionTimeID = dtp.DetectionTimeID_1
	INNER JOIN KS.Detection AS d2 ON 
		d2.DetectionTimeID = dtp.DetectionTimeID_2 AND 
		d2.SkyPatchID = d1.SkyPatchID
	-- Possible values of the SkyPatch
	INNER JOIN KS.SkyPatch AS sp ON sp.SkyPatchID IN (d1.SkyPatchID, d2.SkyPatchID)
WHERE
	-- Only detection time pairs where the DetectionTimePairID is in the selected range
	(dtp.DetectionTimePairID BETWEEN DetectionTimePairID_0 AND (DetectionTimePairID_1-1));

-- Set the distance r and magnitude u of angular velocity
UPDATE KS.TrackletBatch
SET 
	r = SQRT(POW(ux_bar,2)+POW(uy_bar,2)+POW(uz_bar,2)),
	u_dot = SQRT(POW(ux_dot,2)+POW(uy_dot,2)+POW(uz_dot,2));

-- Insert this batch into the main table
INSERT INTO KS.Tracklet
(DetectionTimePairID, DetectionID_1, DetectionID_2, SkyPatchID, 
mjd_bar, mjd1, mjd2,
ux_bar, uy_bar, uz_bar, ux_dot, uy_dot, uz_dot, u_dot,
ux1, uy1, uz1, ux2, uy2, uz2, 
mag_bar, dmag, mag1, mag2)
SELECT
	-- The time pair
	tb.DetectionTimePairID,
	-- The pair of detections
	tb.DetectionID_1,
	tb.DetectionID_2,
	-- The location in the sky
	tb.SkyPatchID,
	-- The three dates
	tb.mjd_bar,
	tb.mjd1,
	tb.mjd2,
	-- The direction
	tb.ux_bar/r AS ux_bar,
	tb.uy_bar/r AS uy_bar,
	tb.uz_bar/r AS uz_bar,
	-- Time derivative of the direction
	tb.ux_dot,
	tb.uy_dot,
	tb.uz_dot,
	tb.u_dot,
	-- Direction of detection 1
	tb.ux1,
	tb.uy1,
	tb.uz1,
	-- Direction of detection 2
	tb.ux2,
	tb.uy2,
	tb.uz2,	
	-- Magnitude
	tb.mag_bar,
	tb.dmag,
	tb.mag1,
	tb.mag2
FROM
	-- The pair of directions for the batch process
	KS.TrackletBatch AS tb
WHERE
	-- Only when the mean direction matches the SkyPatch
	(tb.ux_bar BETWEEN tb.x0*tb.r AND tb.x1*tb.r) AND
	(tb.uy_bar BETWEEN tb.y0*tb.r AND tb.y1*tb.r) AND
	(tb.uz_bar BETWEEN tb.z0*tb.r AND tb.z1*tb.r)
ORDER BY tb.TrackletID;

-- Clean up temp table
DROP TEMPORARY TABLE IF EXISTS KS.DetectionTimePair_present;
DROP TEMPORARY TABLE IF EXISTS KS.DetectionTimePair_insert;
DROP TEMPORARY TABLE IF EXISTS KS.TrackletBatch;

END $$

DELIMITER ;

