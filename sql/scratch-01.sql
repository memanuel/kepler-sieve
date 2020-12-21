CREATE OR REPLACE TABLE KS.TestTime(
	IntegrationTimeID INT NOT NULL
		COMMENT "MJD as integer number of minutes, e.g. floor(MJD*24*60)",
	MJD DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame",
	JD DOUBLE NOT NULL
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	CalendarDate DATE NOT NULL
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	CalendarDateTime DATETIME(6) NOT NULL
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	delta_T DOUBLE NOT NULL
		COMMENT "The difference between the TDB and terrestial atomic time frames"
)
	COMMENT "Distinct time stamps at which MSE integrated positions of solar system bodies are available.";