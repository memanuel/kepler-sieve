CREATE OR REPLACE TABLE KS.Detection(
    -- Primary key and foreign keys
    DetectionTimeID INT NOT NULL
        COMMENT "The time of this detection as a foreign key to the DetectionTime table",
    SkyPatchID INT NOT NULL
    	COMMENT "Foreign key to SkyPatch table",
    k SMALLINT NOT NULL
    	COMMENT "Counter over all the detections sharing this sky patch",
    -- Foreign key fields
    DetectionID INT NOT NULL
    	COMMENT "Foreign key to the RawDetection table",
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
    PRIMARY KEY (DetectionTimeID, SkyPatchID, k)
    	COMMENT "Primary key optimized for spatial searching",
    -- Unique key on DetectionID
    UNIQUE KEY UNQ_Detection_DetectionID(DetectionID)
    	COMMENT "Support joins based on DetectionID",
	-- Foreign keys
    CONSTRAINT FK_Detection_DetectionTimeID
        FOREIGN KEY (DetectionTimeID) REFERENCES KS.DetectionTime(DetectionTimeID),
    CONSTRAINT FK_Detection_DetectionID
        FOREIGN KEY (DetectionID) REFERENCES KS.RawDetection(DetectionID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Detections of possible asteroids across multiple data sources; enriched with unit direction and SkyPatch.";
