CREATE OR REPLACE TABLE KS.AsteroidElement_Ref(
	AsteroidID INT NOT NULL 
		COMMENT "Foreign key to Asteroid table",
	TimeID INT NOT NULL
		COMMENT "Foreign key to IntegrationTime table",
	epoch DOUBLE NOT NULL
		COMMENT "The epoch as of which the orbital elements are computed; expressed as an MJD",
	a DOUBLE NOT NULL
		COMMENT "The semimajor axis in AU",
	e DOUBLE NOT NULL
		COMMENT "The eccentricity; dimensionless",
	inc DOUBLE NOT NULL
		COMMENT "The inclination in radians",
	Omega_node DOUBLE NOT NULL
		COMMENT "The longitude of the ascending node in radians",
	omega_peri DOUBLE NOT NULL
		COMMENT "The argument of perhelion omega in radians",
    f DOUBLE NOT NULL
        COMMENT "The true anomaly f in radians",
	M DOUBLE NOT NULL
		COMMENT "The mean anomaly M in radians",
    eccentric_anomaly DOUBLE NOT NULL
        COMMENT "The eccentric anomaly E in radians",
    period DOUBLE NOT NULL
        COMMENT "The orbital period in days",
    mean_motion DOUBLE NOT NULL
        COMMENT "The mean motion in radians per day",
    PRIMARY KEY (AsteroidID, TimeID)
        COMMENT "Reference elements uniquely defined by the asteroid and the epoch",
    UNIQUE KEY (TimeID, AsteroidID)
        COMMENT "Support fast searching by time ID",
    CONSTRAINT FK_AsteroidElement_Ref_AsteroidID
        FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
    CONSTRAINT FK_AsteroidElement_Ref_TimeID
        FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID)
)
COMMENT "Orbital elements of asteroids as of reference dates; used to initialize integrations. Primary for these elements is the Sun, NOT the Solar System Barycenter!";
