-- ************************************************************************************************
CREATE OR REPLACE TABLE KS.DetectionTimePair(
    -- The pair of detection times
	DetectionTimeID_1 INT NOT NULL
        COMMENT "First detection time in this pair of times; foreign key to DetectionTime table",
	DetectionTimeID_2 INT NOT NULL
        COMMENT "Second detection time in this pair of times",
    -- The data source
    DataSourceID TINYINT NOT NULL
        COMMENT "The data source for this pair of detections; foreign key to DetectionTime table",
    -- The two times
    mjd1 DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the first detection time",
    mjd2 DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the second  detection time",
    -- Primary key
    PRIMARY KEY (DetectionTimeID_1, DetectionTimeID_2), 
    -- Foreign keys
    CONSTRAINT FK_DetectionTimePair_DetectionTimeID_1
        FOREIGN KEY (DetectionTimeID_1) REFERENCES KS.DetectionTime(DetectionTimeID),
    CONSTRAINT FK_DetectionTimePair_DetectionTimeID_2
        FOREIGN KEY (DetectionTimeID_2) REFERENCES KS.DetectionTime(DetectionTimeID),
    CONSTRAINT FK_DetectionTimePair_DataSourceID
        FOREIGN KEY (DataSourceID) REFERENCES KS.DataSource(DataSourceID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Pair of times with available detections close enough to be a tracklet.";