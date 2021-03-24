CREATE OR REPLACE TABLE KS.DetectionTime(
	DetectionTimeID INT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this detection time",
    MJD DOUBLE NOT NULL UNIQUE
        COMMENT "The Modified Julian Date of this detection time in the TDB (barycentric dynamical time) frame",
	CalendarDateTime DATETIME(6) NOT NULL UNIQUE
		COMMENT "The date and time on the Gregorian calendar in the TDB frame of this detection",
    DetectionSourceID TINYINT NOT NULL
        COMMENT "The detection source for this set of detections made at the same time",
    ObservatoryID TINYINT NOT NULL
        COMMENT "The observatory from which the detections were made",
	-- Position q = [qx, qy, qz]; includes topos adjustment
    qx DOUBLE NOT NULL
        COMMENT "Position of observatory (x coordinate) in AU in the barcycentric mean ecliptic frame; includes topos adjustment",
    qy DOUBLE NOT NULL
        COMMENT "Position of observatory (y coordinate) in AU in the barcycentric mean ecliptic frame; includes topos adjustment",
    qz DOUBLE NOT NULL
        COMMENT "Position of observatory (z coordinate) in AU in the barcycentric mean ecliptic frame; includes topos adjustment",
	-- Velocity v = [vx, vy, vz]; includes topos adjustment
	vx DOUBLE NOT NULL
		COMMENT "Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame; includes topos adjustment",
	vy DOUBLE NOT NULL
		COMMENT "Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame; includes topos adjustment",
	vz DOUBLE NOT NULL
		COMMENT "Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame; includes topos adjustment",       
    -- Foreign keys
    CONSTRAINT FK_DetectionSource_DetectionSourceID
        FOREIGN KEY (DetectionSourceID) REFERENCES KS.DetectionSource(DetectionSourceID),
    CONSTRAINT FK_DetectionSource_ObservatoryID
        FOREIGN KEY (ObservatoryID) REFERENCES KS.Observatory(ObservatoryID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Time at which one or more detections were made an observatory.  Cartesian position and velocity includes topos adjustment.";
