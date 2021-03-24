CREATE OR REPLACE TABLE ZTF.ObjectClass(
	ObjectClassID TINYINT NOT NULL PRIMARY KEY,
	ObjectClassCD VARCHAR(16) NOT NULL UNIQUE
		COMMENT "Short code"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Types of bodies described in ZTF alerce system.";

SELECT * FROM ZTF.ObjectClass;