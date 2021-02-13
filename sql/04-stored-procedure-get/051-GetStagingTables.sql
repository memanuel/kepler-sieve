DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetStagingTables(
	IN table_name VARCHAR(64))
COMMENT "Get all temporary staging tables matching the specified parent table name."

BEGIN 

SELECT
	tbl.TABLE_NAME AS table_name_temp,
	-- i0._ AS i0,
	-- i1._ AS i1,
	-- i2._ AS i2,
	-- substring(tbl.TABLE_NAME, i0._, i1._-i0._) AS pid,
	-- substring(tbl.TABLE_NAME, i2._, i2._-i1._) AS chunk,
	pid._ AS pid,
	chunk._ AS chunk
FROM 
	INFORMATION_SCHEMA.TABLES AS tbl
	-- Locate the three delimiter locations for around the pid and chunk number
	INNER JOIN KS.Counter AS i0 ON i0._ = locate('_pid_', tbl.TABLE_NAME) + length('_pid_')
	INNER JOIN KS.Counter AS i1 ON i1._ = locate('_chunk', tbl.TABLE_NAME, i0._)
	INNER JOIN KS.Counter AS i2 ON i2._ = i1._ + length('_chunk_')
	-- Get the pid and chunk number as integers
	INNER JOIN KS.Counter AS pid ON pid._ = CAST(substring(tbl.TABLE_NAME, i0._, i1._-i0._) AS INT)
	INNER JOIN KS.Counter AS chunk ON chunk._ = CAST(substring(tbl.TABLE_NAME, i2._, i2._-i1._) AS INT)
WHERE 
	tbl.TABLE_SCHEMA = 'temp' AND 
	tbl.TABLE_NAME LIKE CONCAT(table_name, '\_pid%\_chunk%')
ORDER BY tbl.table_name;

END
$$

DELIMITER ;
