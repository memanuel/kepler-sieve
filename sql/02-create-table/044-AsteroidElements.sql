-- Calculated orbital elements for asteroids
CREATE OR REPLACE TABLE KS.AsteroidElements(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these orbital elements; FK to KS.DailyTime",
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose orbital elements are described; FK to KS.Asteroid",
	MJD DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
	a DOUBLE NOT NULL COMMENT "The semimajor axis in AU",
	e DOUBLE NOT NULL COMMENT "The eccentricity; dimensionless",
	inc DOUBLE NOT NULL	COMMENT "The inclination in radians",
	Omega_node DOUBLE NOT NULL COMMENT "The longitude of the ascending node in radians",
	omega_peri DOUBLE NOT NULL COMMENT "The argument of periapsis in radians",
	f DOUBLE NOT NULL COMMENT "The true anomaly in radians",
	M DOUBLE NOT NULL COMMENT "The mean anomaly in radians",
	-- Computed columns
	EA DOUBLE AS (MOD(2.0*ATAN(SQRT((1.0-e)/(1.0+e))*TAN(0.5*f))+2.0*PI(), 2.0*PI())) PERSISTENT
		COMMENT "The eccentric anomaly; derived from the true anomaly",
	PRIMARY KEY (TimeID, AsteroidID)
		COMMENT "A orbital elements vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_AsteroidID_TimeID(AsteroidID, TimeID)
		COMMENT "Allow fast search keyed first by AsteroidID.",
	CONSTRAINT FK_AsteroidElements_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_AsteroidElements_AsteroidID
		FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
COMMENT "State vectors (position and velocity) for Solar Systems bodies computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000.";
