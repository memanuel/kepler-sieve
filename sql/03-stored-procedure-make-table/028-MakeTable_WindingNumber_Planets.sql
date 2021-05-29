DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_WindingNumber_Planets()
COMMENT "Populate the WindingNumber table by joining two copies of OrbitalElements_Planets"
BEGIN 

-- Select from StateVectors_Planets into StateVectors
REPLACE INTO KS.WindingNumber_Planets
(TimeID, BodyID, WindingNumber)
SELECT
	elt2.TimeID,
	elt1.BodyID,
	row_number() OVER (PARTITION BY elt1.BodyID ORDER BY elt1.TimeID) AS WindingNumber
FROM
	KS.OrbitalElements_Planets AS elt1
	INNER JOIN KS.OrbitalElements_Planets AS elt2 ON
		elt2.TimeID = elt1.TimeID+5 AND
		elt2.BodyID = elt1.BodyID
WHERE
	elt2.M < elt1.M - 3;

-- Upsample to every available integration time
CREATE OR REPLACE TEMPORARY TABLE KS.WindingNumberBatch
(
	BodyID INT NOT NULL,
	TimeID INT NOT NULL,
	WindingNumber INT NOT NULL,
	PRIMARY KEY (BodyID, TimeID)
);

INSERT INTO KS.WindingNumberBatch
(BodyID, TimeID, WindingNumber)
SELECT
	wn1.BodyID,
	it.TimeID,
	wn1.WindingNumber
FROM
	KS.WindingNumber_Planets AS wn1
	INNER JOIN KS.WindingNumber_Planets AS wn2 ON
		wn2.BodyID = wn1.BodyID AND
		wn2.WindingNumber = wn1.WindingNumber+1
	INNER JOIN KS.IntegrationTime AS it ON 
		wn1.TimeID <= it.TimeID AND it.TimeID < wn2.TimeID;
		
-- Join the winding number back to the main OrbitalElements table
UPDATE
	KS.OrbitalElements_Planets AS elts
	INNER JOIN KS.WindingNumberBatch AS wn ON
		wn.BodyID = elts.BodyID AND
		wn.TimeID = elts.TimeID
SET
	elts.WindingNumber = wn.WindingNumber;

-- Previous join only handles times in between a pair of winding times
-- Do one last update for times on or after the last winding    
CREATE OR REPLACE TEMPORARY TABLE KS.t1 AS
SELECT
	wn.BodyID,
	max(wn.TimeID) AS TimeID,
	max(wn.WindingNumber) AS WindingNumber
FROM
	KS.WindingNumber_Planets AS wn
GROUP BY wn.BodyID;

UPDATE
	KS.OrbitalElements_Planets AS el
	INNER JOIN KS.t1 ON t1.BodyID = el.BodyID 
SET
	el.WindingNumber = t1.WindingNumber
WHERE
	el.TimeID >= t1.TimeID;

DROP TEMPORARY TABLE IF EXISTS KS.t1;
  
END
$$

DELIMITER ;
