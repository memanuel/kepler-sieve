CREATE OR REPLACE TABLE JPL.HorizonsTime(
	TimeID INT NOT NULL PRIMARY KEY
        COMMENT "MJD as integer number of minutes, e.g. floor(MJD*24*60)",
	MJD DOUBLE NOT NULL UNIQUE
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame",
	JD DOUBLE AS (MJD + 2400000.5)
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	CalendarDate DATE NOT NULL
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	CalendarDateTime DATETIME(6) NOT NULL UNIQUE
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	delta_T DOUBLE NOT NULL
		COMMENT "The difference between the TDB and terrestial atomic time frames"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Distinct time stamps at which Horizons data is available.";
