CREATE OR REPLACE TABLE ZTF.ObjectClass(
	ObjectClassID TINYINT NOT NULL PRIMARY KEY,
	ObjectClassName VARCHAR(16) NOT NULL UNIQUE
		COMMENT "Name or acronym of this object class"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Types of bodies described in ZTF alerce system.";
