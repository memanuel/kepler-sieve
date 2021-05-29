DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.MakeTable_AsteroidWindingNumber()
COMMENT "Populate the WindingNumber table by joining two copies of AsteroidElements"
BEGIN 

REPLACE INTO KS.AsteroidWindingNumber
(TimeID, AsteroidID, WindingNumber)
SELECT
	elt2.TimeID,
	elt1.AsteroidID,
	row_number() OVER (PARTITION BY elt1.AsteroidID ORDER BY elt1.TimeID) AS WindingNumber
FROM
	KS.AsteroidElements AS elt1
	INNER JOIN KS.AsteroidElements AS elt2 ON
		elt2.TimeID = elt1.TimeID + (4*1440) AND
		elt2.AsteroidID = elt1.AsteroidID
WHERE
	elt2.M < elt1.M - 3;

-- Upsample to every available integration time
CREATE OR REPLACE TEMPORARY TABLE KS.WindingNumberBatch
(
	AsteroidID INT NOT NULL,
	TimeID INT NOT NULL,
	WindingNumber INT NOT NULL,
	PRIMARY KEY (AsteroidID, TimeID)
);

INSERT INTO KS.WindingNumberBatch
(AsteroidID, TimeID, WindingNumber)
SELECT
	wn1.AsteroidID,
	dt.TimeID,
	wn1.WindingNumber
FROM
	KS.AsteroidWindingNumber AS wn1
	INNER JOIN KS.AsteroidWindingNumber AS wn2 ON
		wn2.AsteroidID = wn1.AsteroidID AND
		wn2.WindingNumber = wn1.WindingNumber+1
	INNER JOIN KS.DailyTime AS dt ON 
		wn1.TimeID <= dt.TimeID AND dt.TimeID < wn2.TimeID;
		
-- Join the winding number back to the main OrbitalElements table
UPDATE
	KS.AsteroidElements AS elts
	INNER JOIN KS.WindingNumberBatch AS wn ON
		wn.AsteroidID = elts.AsteroidID AND
		wn.TimeID = elts.TimeID
SET
	elts.WindingNumber = wn.WindingNumber;

-- Previous update only handles times in between a pair of winding times
-- Do one last update for times on or after the last winding    
CREATE OR REPLACE TEMPORARY TABLE KS.t1 AS
SELECT
	wn.AsteroidID,
	max(wn.TimeID) AS TimeID,
	max(wn.WindingNumber) AS WindingNumber
FROM
	KS.AsteroidWindingNumber AS wn
GROUP BY wn.AsteroidID;

UPDATE
	KS.AsteroidElements AS el
	INNER JOIN KS.t1 ON t1.AsteroidID = el.AsteroidID 
SET
	el.WindingNumber = t1.WindingNumber
WHERE
	el.TimeID >= t1.TimeID;

DROP TEMPORARY TABLE IF EXISTS KS.t1;

  
END
$$

DELIMITER ;
