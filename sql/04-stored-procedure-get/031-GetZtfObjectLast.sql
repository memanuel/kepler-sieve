DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE ZTF.GetObjectLast()
COMMENT "Get the last object code in the ZTF.Object table."

BEGIN 

SELECT
	COALESCE(max(obj.ObjectCD), 'ZTF17') AS LastObjectCD
FROM
	ZTF.`Object` AS obj;

END
$$

DELIMITER ;
