-- Create tables for Horizons vectors
CREATE OR REPLACE TABLE JPL.HorizonsVectors(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to JPL.HorizonsTime",
	HorizonsBodyID INT NOT NULL
		COMMENT "The HorizonsBody whose state vectors are described; FK to JPL.HorizonsBody",
	MJD DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from HorizonsTimeID but included for performance.",
	qx DOUBLE NOT NULL
		COMMENT "Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame",
	qy DOUBLE NOT NULL
		COMMENT "Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame",
	qz DOUBLE NOT NULL
		COMMENT "Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame",
	vx DOUBLE NOT NULL
		COMMENT "Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vy DOUBLE NOT NULL
		COMMENT "Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vz DOUBLE NOT NULL
		COMMENT "Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame",
	PRIMARY KEY (TimeID, HorizonsBodyID)
		COMMENT "A state vector is identified by the body and time stamp; use integer TimeStampID for performance.",
	UNIQUE KEY (HorizonsBodyID, TimeID)
		COMMENT "Support searching by HorizonsBodyID first."
)
COMMENT "State vectors (position and velocity) for Solar Systems bodies downloaded from JPL Horizons web server."
