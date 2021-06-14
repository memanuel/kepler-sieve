-- State vectors from MSE integration using Rebound
CREATE OR REPLACE TABLE KS.StateVectors_Topos(
    ObservatoryID TINYINT NOT NULL
        COMMENT "The observatory whose vectors are described; FK to KS.Observatory",
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to KS.HiResTime",
	mjd DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
	-- Position of observatory on Earth, dq = [dqx, dqy, dqz]
	qx DOUBLE NOT NULL
		COMMENT "Position of observatory on Earth (x coordinate); offset to geocenter in BME frame; in AU",
	qy DOUBLE NOT NULL
		COMMENT "Position of observatory on Earth (y coordinate); offset to geocenter in BME frame; in AU",
	qz DOUBLE NOT NULL
		COMMENT "Position of observatory on Earth (z coordinate); offset to geocenter in BME frame; in AU",
	-- Velocity of observatory on Earth dv = [dvx, dvy, dvz]
	vx DOUBLE NOT NULL
		COMMENT "Velocity of observatory on Earth (x coordinate); offset to geocenter in BME frame; in AU/day",
	vy DOUBLE NOT NULL
		COMMENT "Velocity of observatory on Earth (y coordinate); offset to geocenter in BME frame; in AU/day",
	vz DOUBLE NOT NULL
		COMMENT "Velocity of observatory on Earth (z coordinate); offset to geocenter in BME frame; in AU/day",
	-- Keys and constraints
	PRIMARY KEY (ObservatoryID, TimeID)
		COMMENT "A topos state vector is identified by the observatory and time stamp; use integer time ID for performance.",
	CONSTRAINT FK_StateVectors_Topos_ObservatoryID
		FOREIGN KEY (ObservatoryID) REFERENCES KS.Observatory(ObservatoryID),
	CONSTRAINT FK_StateVectors_Topos_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.HiResTime(TimeID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "State vectors (position and velocity) for Earth computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000."
