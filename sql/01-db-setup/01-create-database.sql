CREATE OR REPLACE DATABASE KS
	CHARACTER SET = 'latin1'
	COLLATE = 'latin1_general_cs'
	COMMENT = "Main database for Kepler Sieve application.";

CREATE OR REPLACE DATABASE JPL
	CHARACTER SET = 'latin1'
	COLLATE = 'latin1_general_cs'
	COMMENT = "Database for data orginating from the NASA Jet Propulsion Laboratory (JPL)";

CREATE OR REPLACE DATABASE ZTF
	CHARACTER SET = 'latin1'
	COLLATE = 'latin1_general_cs'
	COMMENT = "Database for data orginating ZTF-2 (Zwicky Transient Facility) processed by Alerce data broker";

CREATE OR REPLACE DATABASE temp
	CHARACTER SET = 'latin1'
	COLLATE = 'latin1_general_cs'
	COMMENT = "Temporary database.";