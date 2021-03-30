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
