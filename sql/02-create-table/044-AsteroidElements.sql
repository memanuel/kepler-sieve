-- Calculated orbital elements for asteroids
CREATE OR REPLACE TABLE KS.AsteroidElements(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these orbital elements; FK to KS.DailyTime",
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose orbital elements are described; FK to KS.Asteroid",
	mjd DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
	a DOUBLE NOT NULL COMMENT "The semimajor axis in AU",
	e DOUBLE NOT NULL COMMENT "The eccentricity; dimensionless",
	inc DOUBLE NOT NULL	COMMENT "The inclination in radians",
	Omega_node DOUBLE NOT NULL COMMENT "The longitude of the ascending node in radians",
	omega_peri DOUBLE NOT NULL COMMENT "The argument of periapsis in radians",
	f DOUBLE NOT NULL COMMENT "The true anomaly in radians",
	M DOUBLE NOT NULL COMMENT "The mean anomaly in radians",
    WindingNumber INT NOT NULL DEFAULT 0
        COMMENT "The number of complete orbits; so M can be adjusted to be monotonic for splining",
	PRIMARY KEY (TimeID, AsteroidID)
		COMMENT "A orbital elements vector is identified by the body and time stamp; use integer time ID for performance.",
	UNIQUE KEY UNQ_AsteroidID_TimeID(AsteroidID, TimeID)
        COMMENT "Allow fast search keyed first by AsteroidID.",
	CONSTRAINT FK_AsteroidElements_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_AsteroidElements_AsteroidID FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Keplerian orbital elements (a, e, i, Omega, omega, f) for all known asteroids computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000.";
