CREATE OR REPLACE TABLE KS.DetectionSource(
	DetectionSourceID TINYINT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this detection source",
    DetectionSourceCD VARCHAR(3) NOT NULL UNIQUE
        COMMENT "Short code for this detection source",
	DetectionSourceName VARCHAR(32) NOT NULL UNIQUE
		COMMENT "Name of this source of detections",
    ObservatoryID TINYINT NOT NULL
        COMMENT "Observatory from which these detections were made; foreign key to Observatory table",
	SortOrder TINYINT NOT NULL,
    CONSTRAINT FK_DetectionSource_ObservatoryID
        FOREIGN KEY (ObservatoryID) REFERENCES KS.Observatory(ObservatoryID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Sources of asteroid detections.";

INSERT INTO KS.DetectionSource
(DetectionSourceID, DetectionSourceCD, DetectionSourceName, ObservatoryID, SortOrder)
VALUES
(1, 'ZTF', 'Zwicky Transient Facility', 1, 1);
