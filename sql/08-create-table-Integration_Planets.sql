-- Create tables for Horizons vectors
CREATE OR REPLACE TABLE KS.Integration_Planets(
	BodyID INT NOT NULL
		COMMENT "The Body whose state vectors are described; FK to JS.Body",
	IntegrationTimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to KS.IntegrationTime",
	MinuteID INT NOT NULL
		COMMENT "MJD as integer number of minutes, e.g. floor(MJD*24*60)",
	MJD DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
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
	PRIMARY KEY (BodyID, IntegrationTimeID)
		COMMENT "A state vector is identified by the body and time stamp; use integer TimeStampID for performance.",
	UNIQUE KEY (MinuteID, BodyID)
		COMMENT "Times are uniquely identified by minute; allow fast search by MinuteID.",
	INDEX MJD_BodyID (MJD, BodyID)
		COMMENT "Index to support filtering directly on the time as an MJD.",
	CONSTRAINT FK_Integration_Planets_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID),
	CONSTRAINT FK_Integration_Planets_IntegrationTimeID
		FOREIGN KEY (IntegrationTimeID) REFERENCES KS.IntegrationTime(IntegrationTimeID)
)
COMMENT "State vectors (position and velocity) for Solar Systems bodies computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000."
