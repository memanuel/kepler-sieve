CREATE OR REPLACE TABLE KS.DetectionBatch(
    -- Primary key
    DetectionBatchID INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    -- Foreign key fields
    DetectionTimeID INT NOT NULL
        COMMENT "The time of this detection as a foreign key to the DetectionTime table",
    DetectionID INT NOT NULL
    	COMMENT "Foreign key to the RawDetection table",
	-- Unit direction u = [ux, uy, uz] in mean ecliptic frame
    ux DOUBLE NOT NULL
        COMMENT "Unit direction (x coordinate) in the barcycentric mean ecliptic frame",
    uy DOUBLE NOT NULL
        COMMENT "Unit direction (y coordinate) in the barcycentric mean ecliptic frame",
    uz DOUBLE NOT NULL
        COMMENT "Unit direction (z coordinate) in the barcycentric mean ecliptic frame",
    -- The apparent magnitude
    mag DOUBLE NOT NULL
    	COMMENT "Apparent magnitude (brightness) of this detection"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Staging table for insterting to Detections table; SkyPatchID and k to be calculated on insertion.";
