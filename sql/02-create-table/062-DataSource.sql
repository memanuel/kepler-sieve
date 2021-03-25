CREATE OR REPLACE TABLE KS.DataSource(
	DataSourceID TINYINT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this data source",
    DataSourceCD VARCHAR(3) NOT NULL UNIQUE
        COMMENT "Short code for this data source",
	DataSourceName VARCHAR(32) NOT NULL UNIQUE
		COMMENT "Name of this data source",
    ObservatoryID TINYINT NOT NULL
        COMMENT "Observatory from which this data (e.g. detections) were obtained; foreign key to Observatory table",
	SortOrder TINYINT NOT NULL,
    CONSTRAINT FK_DataSource_ObservatoryID
        FOREIGN KEY (ObservatoryID) REFERENCES KS.Observatory(ObservatoryID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Data sources for asteroid detections.";

INSERT INTO KS.DataSource
(DataSourceID, DataSourceCD, DataSourceName, ObservatoryID, SortOrder)
VALUES
(1, 'ZTF', 'Zwicky Transient Facility', 1, 1);
