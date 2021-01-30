DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetAsteroidRefElementDates(
	IN epoch INT)
COMMENT "Get all dates as of which reference asteroid orbital elements are available, but not as of the desired epoch."

BEGIN 

-- Compute TimeID from epoch
SET @TimeID = epoch * 24 * 60;

-- Select all available orbital elements at this time, sorted by AsteroidID.
SELECT
	elt.TimeID,
	elt.epoch,
	COUNT(elt.TimeID) AS AsteroidCount
FROM
	KS.AsteroidElement_Ref AS elt
WHERE NOT EXISTS (
	SELECT elt2.AsteroidID FROM KS.AsteroidElement_Ref AS elt2 
	WHERE elt2.AsteroidID = elt.AsteroidID AND elt2.TimeID = @TimeID
	)
GROUP BY elt.TimeID
ORDER BY elt.TimeID;	

END
$$

DELIMITER ;
