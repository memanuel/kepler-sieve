-- ************************************************************************************************
-- drop tables in reverse dependency order

-- ************************************************************************************************
CREATE OR REPLACE TABLE KS.DetectionTimeSlice(
	DetectionTimeSliceID INT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for a slice of detection times",
    IntegrationTimeID INT NOT NULL
    	COMMENT "The last IntegrationTime prior to or equal to the start of the time slice",
    mjd DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the midpoint of this time slice",
    mjd0 DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the start of this time slice (inclusive)",
    mjd1 DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the end of this time slice (exclusive)",
    -- Keys
    UNIQUE KEY UNQ_mjd0_mjd1 (mjd0, mjd1),
    INDEX IDX_IntegrationTimeID (IntegrationTimeID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Time at which one or more detections were made an observatory.  Cartesian position and velocity includes topos adjustment.";

-- ************************************************************************************************
CREATE OR REPLACE TABLE KS.DetectionTime(
	DetectionTimeID INT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this detection time",
    DetectionTimeSliceID INT NOT NULL
    	COMMENT "The time slice; foreign key to DetectionTimeSlice",
    mjd DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of this detection time in the TDB (barycentric dynamical time) frame",
	CalendarDateTime DATETIME(6) NOT NULL UNIQUE
		COMMENT "The date and time on the Gregorian calendar in the TDB frame of this detection",
    DataSourceID TINYINT NOT NULL
        COMMENT "The data source for this set of detections made at the same time; foreign key to DataSource table",
    ObservatoryID TINYINT NOT NULL
        COMMENT "The observatory from which the detections were made; foreign key to Observatory table",
	-- Position q = [qx, qy, qz]; includes topos adjustment
    qObs_x DOUBLE NOT NULL
        COMMENT "Position of observatory (x coordinate) in AU in the barcycentric mean ecliptic frame; includes topos adjustment",
    qObs_y DOUBLE NOT NULL
        COMMENT "Position of observatory (y coordinate) in AU in the barcycentric mean ecliptic frame; includes topos adjustment",
    qObs_z DOUBLE NOT NULL
        COMMENT "Position of observatory (z coordinate) in AU in the barcycentric mean ecliptic frame; includes topos adjustment",
	-- Velocity v = [vx, vy, vz]; includes topos adjustment
	vObs_x DOUBLE NOT NULL
		COMMENT "Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame; includes topos adjustment",
	vObs_y DOUBLE NOT NULL
		COMMENT "Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame; includes topos adjustment",
	vObs_z DOUBLE NOT NULL
		COMMENT "Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame; includes topos adjustment",       
	-- Position of the sun Sun = [Sun_x, Sun_y, Sun_z]
    qSun_x DOUBLE NOT NULL
        COMMENT "Position of of the sun (x coordinate) in AU in the barcycentric mean ecliptic frame",
    qSun_y DOUBLE NOT NULL
        COMMENT "Position of of the sun (y coordinate) in AU in the barcycentric mean ecliptic frame",
    qSun_z DOUBLE NOT NULL
        COMMENT "Position of of the sun (z coordinate) in AU in the barcycentric mean ecliptic frame",
    -- Foreign keys
    CONSTRAINT FK_DetectionTime_DataSourceID
        FOREIGN KEY (DataSourceID) REFERENCES KS.DataSource(DataSourceID),
    CONSTRAINT FK_DetectionTime_ObservatoryID
        FOREIGN KEY (ObservatoryID) REFERENCES KS.Observatory(ObservatoryID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Time at which one or more detections were made an observatory.  Cartesian position and velocity includes topos adjustment.";


-- ************************************************************************************************
CREATE OR REPLACE TABLE KS.DetectionTimePair(
	DetectionTimeID_1 INT NOT NULL
        COMMENT "First detection time in this pair of times; foreign key to DetectionTime table",
	DetectionTimeID_2 INT NOT NULL
        COMMENT "Second detection time in this pair of times",
    DataSourceID TINYINT NOT NULL
        COMMENT "The data source for this pair of detections; foreign key to DetectionTime table",
    mjd1 DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the first detection time",
    mjd2 DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of the second  detection time",

    -- Foreign keys
    CONSTRAINT FK_DetectionTimePair_DetectionTimeID_1
        FOREIGN KEY (DetecitonTimeID_1) REFERENCES KS.DetectionTime(DetectionTimeID),
    CONSTRAINT FK_DetectionTimePair_DetectionTimeID_2
        FOREIGN KEY (DetecitonTimeID_2) REFERENCES KS.DetectionTime(DetectionTimeID),
    CONSTRAINT FK_DetectionTimePair_DataSourceID
        FOREIGN KEY (DataSourceID) REFERENCES KS.DataSource(DataSourceID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Time at which one or more detections were made an observatory.  Cartesian position and velocity includes topos adjustment.";