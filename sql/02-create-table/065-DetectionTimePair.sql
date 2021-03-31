-- ************************************************************************************************
CREATE OR REPLACE TABLE KS.DetectionTimePair(
	-- Single integer key on pairs of times
	DetectionTimePairID INT NOT NULL AUTO_INCREMENT PRIMARY KEY
		COMMENT "Pair of detection times is a logical key, but a single integer ID on pairs helps performance building Tracklet table",
    -- The pair of detection times
	DetectionTimeID_1 INT NOT NULL
        COMMENT "First detection time in this pair of times; foreign key to DetectionTime table",
	DetectionTimeID_2 INT NOT NULL
        COMMENT "Second detection time in this pair of times",
    -- The data source
    DataSourceID TINYINT NOT NULL
        COMMENT "The data source for this pair of detections; foreign key to DetectionTime table",
    -- The two times
    mjd1 DOUBLE NOT NULL
        COMMENT "The Modified Julian Date of the first detection time",
    mjd2 DOUBLE NOT NULL
        COMMENT "The Modified Julian Date of the second detection time",
    -- Mean detection time
    mjd_bar DOUBLE NOT NULL
    	COMMENT "The Modified Julian Date of the mean detection time",
    -- Difference in two detection times
    dt DOUBLE NOT NULL
    	COMMENT "Difference in the two detection times as an mjd",
    -- Primary key
    -- Unique key
    UNIQUE KEY UNQ_DetectionTimePair_DetectionID_1_2(DetectionTimeID_1, DetectionTimeID_2)
    	COMMENT "The pair of detection Time IDs is unique", 
    UNIQUE KEY UNQ_DetectionTimePair_mjd_12(mjd1, mjd2)
    	COMMENT "The pair of detection times is unique",
    -- Index
    INDEX IDX_DetectionTimePair_mjd(mjd)
    	COMMENT "Regular index on mean detection time; empirically found that mean detection time was not unique",
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
