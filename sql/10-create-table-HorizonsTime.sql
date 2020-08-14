CREATE OR REPLACE TABLE JPL.HorizonsTime(
	HorizonsTimeID INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    IntegrationTimeID INT NULL UNIQUE
        COMMENT "Foreign key to KS.IntegrationTime",
	MinuteID INT NOT NULL UNIQUE
		COMMENT "MJD as integer number of minutes, e.g. floor(MJD*24*60)",
	MJD DOUBLE NOT NULL UNIQUE
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame",
	JD DOUBLE NOT NULL UNIQUE
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	CalendarDate DATE NOT NULL
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	CalendarDateTime DATETIME(6) NOT NULL UNIQUE
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	delta_T DOUBLE NOT NULL
		COMMENT "The difference between the TDB and terrestial atomic time frames"
)
	COMMENT "Distinct time stamps at which Horizons data is available.";
