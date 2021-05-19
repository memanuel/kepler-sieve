-- Calculated orbital elements for planets integration
CREATE OR REPLACE TABLE KS.OrbitalElements_Planets(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these orbital elements; FK to KS.IntegrationTime",
	BodyID INT NOT NULL
		COMMENT "The Body whose orbital elements are described; FK to KS.Body",
	mjd DOUBLE NOT NULL
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
	PRIMARY KEY (TimeID, BodyID)
		COMMENT "A orbital elements vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_BodyID_TimeID(BodyID, TimeID)
		COMMENT "Allow fast search keyed first by BodyID.",
	CONSTRAINT FK_OrbitalElements_Planets_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_OrbitalElements_Planets_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Keplerian orbital elements (a, e, i, Omega, omega, f) for Solar Systems bodies computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000.  Includes records for the Sun, Earth, Moon, and planet barycenters.";
