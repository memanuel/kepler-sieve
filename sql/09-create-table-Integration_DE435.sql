-- Create tables for integrated vectors
CREATE OR REPLACE TABLE KS.Integration_DE435(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to KS.IntegrationTime",
	BodyID INT NOT NULL
		COMMENT "The Body whose state vectors are described; FK to JS.Body",
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
	PRIMARY KEY (TimeID, BodyID)
		COMMENT "A state vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY (BodyID, TimeID)
		COMMENT "Allow fast search keyed first by BodyID.",
	CONSTRAINT FK_Integration_DE435_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_Integration_DE435_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
COMMENT "State vectors (position and velocity) for Solar Systems bodies computed in rebound using all the massive bodies from the DE435 integration with initial conditions at MJD 59000."
