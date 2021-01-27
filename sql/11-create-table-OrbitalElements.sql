-- Create tables for calculated orbital elements
CREATE OR REPLACE TABLE KS.OrbitalElements(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these orbital elements; FK to KS.IntegrationTime",
	BodyID INT NOT NULL
		COMMENT "The Body whose orbital elements are described; FK to KS.Body",
	MJD DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
	a DOUBLE NOT NULL COMMENT "The semimajor axis in AU",
	e DOUBLE NOT NULL COMMENT "The eccentricity; dimensionless",
	inc DOUBLE NOT NULL	COMMENT "The inclination in radians",
	Omega_node DOUBLE NOT NULL COMMENT "The longitude of the ascending node in radians",
	omega_peri DOUBLE NOT NULL COMMENT "The argument of periapsis in radians",
	f DOUBLE NOT NULL COMMENT "The true anomaly in radians",
	M DOUBLE NOT NULL COMMENT "The mean anomaly in radians",
	PRIMARY KEY (TimeID, BodyID)
		COMMENT "A orbital elements vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_BodyID_TimeID(BodyID, TimeID)
		COMMENT "Allow fast search keyed first by BodyID.",
	CONSTRAINT FK_OrbitalElements_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_OrbitalElements_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
COMMENT "State vectors (position and velocity) for Solar Systems bodies computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000.";
