CREATE OR REPLACE TABLE KS.RawDetection(
    -- Primary key
	DetectionID INT NOT NULL AUTO_INCREMENT PRIMARY KEY
        COMMENT "Integer ID for this detection; unique across multiple sources",
    -- Foreign key fields
    DataSourceID TINYINT NOT NULL
        COMMENT "The data source for this detection; foreign key to DataSource table",
    ObservatoryID TINYINT NOT NULL
        COMMENT "The observatory from which the detections were made; foreign key to Observatory table",
    DetectionTimeID INT NOT NULL
        COMMENT "The time of this detection as a foreign key to the DetectionTime table",
    SourceDetectionID BIGINT NOT NULL
    	COMMENT "Integer ID for this detection according to its data source, e.g. ZTF.Detection.DetectionID for ZTF detections",
    -- Astrometric direction
    -- mjd DOUBLE NOT NULL,
    ra DOUBLE NOT NULL
        COMMENT "Right Ascension in degrees",
    `dec` DOUBLE NOT NULL
        COMMENT "Declination in degrees",
    -- Unique keys
    UNIQUE KEY UNQ_RawDetection_DataSourceID_SourceDetectionID(DataSourceID, SourceDetectionID)
    	COMMENT "Within each source, the SourceDetectionID is unique; this supports fast search by SourceDetectionID",
    -- Foreign keys
    CONSTRAINT FK_RawDetection_DataSourceID
        FOREIGN KEY (DataSourceID) REFERENCES KS.DataSource(DataSourceID),
    CONSTRAINT FK_RawDetection_ObservatoryID
        FOREIGN KEY (ObservatoryID) REFERENCES KS.Observatory(ObservatoryID),
    CONSTRAINT FK_RawDetection_DetectionTimeID
        FOREIGN KEY (DetectionTimeID) REFERENCES KS.DetectionTime(DetectionTimeID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Detections of possible asteroids across multiple data sources, as quoted in RA/DEC by original source.";
