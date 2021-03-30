CREATE OR REPLACE TABLE KS.Tracklet(
    -- Primary key
    TrackletID INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY
        COMMENT "Integer ID for the Tracklet table",
    -- The two detections
    DetectionID_1 INT NOT NULL
        COMMENT "The first detection in this tracklet; foreign key to KS.Detection",
    DetectionID_2 INT NOT NULL
        COMMENT "The second detection in this tracklet; foreign key to KS.Detection",
    -- The SkyPatch of this detection
    SkyPatchID INT NOT NULL
        COMMENT "The SkyPatch of the mean detection; foreign key to SkyPatch table",
    -- Data payload
    mjd DOUBLE NOT NULL
        COMMENT "The mean of the two detection times",
	-- Mean unit direction u = [ux, uy, uz] in mean ecliptic frame
    ux DOUBLE NOT NULL
        COMMENT "Mean unit direction (x coordinate) of the two detections",
    uy DOUBLE NOT NULL
        COMMENT "Mean unit direction (y coordinate) of the two detections",
    uz DOUBLE NOT NULL
        COMMENT "Mean unit direction (z coordinate) of the two detections",
	-- Time derivative of direction 
    ux_dot DOUBLE NOT NULL
        COMMENT "Time derivative of unit direction (x coordinate) of the two detections",
    uy_dot DOUBLE NOT NULL
        COMMENT "Time derivative of unit direction (y coordinate) of the two detections",
    uz_dot DOUBLE NOT NULL
        COMMENT "Time derivative of unit direction (z coordinate) of the two detections",
    -- The apparent magnitude
    mag DOUBLE NOT NULL
    	COMMENT "Mean apparent magnitude (brightness) of this detection",
    -- Primary key
    UNIQUE KEY (DetectionID_1, DetectionID_2)
    	COMMENT "The pair of detection IDs uniquely identifies a tracklet",
	-- Foreign keys
    CONSTRAINT FK_Tracklet_DetectionID_1
        FOREIGN KEY (DetectionID_1) REFERENCES KS.Detection(DetectionID),
    CONSTRAINT FK_Tracklet_DetectionID_2
        FOREIGN KEY (DetectionID_2) REFERENCES KS.Detection(DetectionID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Tracklets: pair of two detections close together in time and in the sky.";
