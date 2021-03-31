CREATE OR REPLACE TABLE KS.Tracklet(
    -- Primary key
    TrackletID INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY
        COMMENT "Integer ID for the Tracklet table",
    -- The two detections
    DetectionID_1 INT NOT NULL COMMENT "The first detection in this tracklet; foreign key to KS.Detection",
    DetectionID_2 INT NOT NULL COMMENT "The second detection in this tracklet; foreign key to KS.Detection",
    -- The pair of detection times being linked
    DetectionTimePairID INT NOT NULL COMMENT "The pair of detection times linked up; foreign key to KS.DetectionTimePair",
    -- The SkyPatch of this detection
    SkyPatchID INT NOT NULL COMMENT "The SkyPatch of the mean detection; foreign key to SkyPatch table",
    -- Three detection times
    mjd_bar DOUBLE NOT NULL COMMENT "The mean of the two detection times",
    mjd1 double NOT NULL COMMENT "Detection time 1",
    mjd2 double NOT NULL COMMENT "Detection time 2",    	
	-- Mean unit direction u_bar = [ux, uy, uz] in mean ecliptic frame
    ux_bar DOUBLE NOT NULL COMMENT "Mean unit direction (x coordinate) of the two detections",
    uy_bar DOUBLE NOT NULL COMMENT "Mean unit direction (y coordinate) of the two detections",
    uz_bar DOUBLE NOT NULL COMMENT "Mean unit direction (z coordinate) of the two detections",
	-- Time derivative of direction 
    ux_dot DOUBLE NOT NULL COMMENT "Time derivative of unit direction (x coordinate) of the two detections",
    uy_dot DOUBLE NOT NULL COMMENT "Time derivative of unit direction (y coordinate) of the two detections",
    uz_dot DOUBLE NOT NULL COMMENT "Time derivative of unit direction (z coordinate) of the two detections",
    -- Magnitude of angular velocity
    u_dot double NOT NULL COMMENT "Magnitude of time derivative of direction vector u; in radians / day for small angles",
	-- Direction of detection 1
    ux1 DOUBLE NOT NULL COMMENT "Unit direction (x coordinate) of detection 1",
    uy1 DOUBLE NOT NULL COMMENT "Unit direction (y coordinate) of detection 1",
    uz1 DOUBLE NOT NULL COMMENT "Unit direction (z coordinate) of detection 1",
	-- Direction of detection 2
    ux2 DOUBLE NOT NULL COMMENT "Unit direction (x coordinate) of detection 2",
    uy2 DOUBLE NOT NULL COMMENT "Unit direction (y coordinate) of detection 2",
    uz2 DOUBLE NOT NULL COMMENT "Unit direction (z coordinate) of detection 2",
    -- The apparent magnitude
    mag_bar DOUBLE NOT NULL COMMENT "Mean apparent magnitude (brightness) of two detections",
    dmag DOUBLE NOT NULL COMMENT "Change in apparent magnitude (brightness) of two detections",
    mag1 DOUBLE NOT NULL COMMENT "Apparent magnitude (brightness) of detection 1",
    mag2 DOUBLE NOT NULL COMMENT "Apparent magnitude (brightness) of detection 2",
    -- Unique key
    UNIQUE KEY UNQ_Tracklet_DetectionID_12(DetectionID_1, DetectionID_2)
    	COMMENT "The pair of detection IDs uniquely identifies a tracklet",
	-- Foreign keys
	CONSTRAINT FK_Tracklet_DetectionTimePairID
		FOREIGN KEY (DetectionTimePairID) REFERENCES KS.DetectionTimePair(DetectionTimePairID),
    CONSTRAINT FK_Tracklet_DetectionID_1
        FOREIGN KEY (DetectionID_1) REFERENCES KS.Detection(DetectionID),
    CONSTRAINT FK_Tracklet_DetectionID_2
        FOREIGN KEY (DetectionID_2) REFERENCES KS.Detection(DetectionID),
    CONSTRAINT FK_Tracklet_SkyPatch
    	FOREIGN KEY (SkyPatchID) REFERENCES KS.SkyPatch(SkyPatchID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Tracklets: pair of two detections close together in time and in the sky.";
