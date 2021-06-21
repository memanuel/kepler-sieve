CREATE OR REPLACE TABLE KS.Detection(
    -- Primary key
    DetectionID INT NOT NULL
    	COMMENT "Primary key; also foreign key to the RawDetection table",
    -- Foreign keys
    DetectionTimeID INT NOT NULL
        COMMENT "The time of this detection as a foreign key to the DetectionTime table",
    -- Spatio-temporal index on (SkyPatchID, TimeID)
    SkyPatchID INT NOT NULL
    	COMMENT "Foreign key to SkyPatch table",
    TimeID INT NOT NULL
    	COMMENT "Integer ID for the time in minutes; foreign key to HiResTime table",
	-- Data payload
    mjd DOUBLE NOT NULL
        COMMENT "The Modified Julian Date of this detection in the TDB (barycentric dynamical time) frame",
	-- Unit direction u = [ux, uy, uz] in mean ecliptic frame
    ux DOUBLE NOT NULL
        COMMENT "Unit direction (x coordinate) in the barcycentric mean ecliptic frame",
    uy DOUBLE NOT NULL
        COMMENT "Unit direction (y coordinate) in the barcycentric mean ecliptic frame",
    uz DOUBLE NOT NULL
        COMMENT "Unit direction (z coordinate) in the barcycentric mean ecliptic frame",
    -- The apparent magnitude
    mag DOUBLE NOT NULL
    	COMMENT "Apparent magnitude (brightness) of this detection",
    -- Primary key
    PRIMARY KEY (DetectionID)
    	COMMENT "Primary key optimized for speed; foreign key to RawDetection table",
    INDEX IDX_Detection_SkyPatchID_DetectionTimeID (SkyPatchID, TimeID) 
    	COMMENT "Support matching detections by location and time; SkyPatchID is first b/c it has higher cardinality",
	-- Foreign keys
    CONSTRAINT FK_Detection_DetectionID
        FOREIGN KEY (DetectionID) REFERENCES KS.RawDetection(DetectionID)
    -- Skip FKs to avoid bloating the table
--    CONSTRAINT FK_Detection_DetectionTimeID
--        FOREIGN KEY (DetectionTimeID) REFERENCES KS.DetectionTime(DetectionTimeID)
--    CONSTRAINT FK_Detection_TimeID
--        FOREIGN KEY (TimeID) REFERENCES KS.HiResTime(TimeID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Detections of possible asteroids across multiple data sources; enriched with unit direction and SkyPatch.";
