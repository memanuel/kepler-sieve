DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE KS.GetDataSources()
COMMENT "Get all the data sources."

BEGIN 

SELECT
	ds.DataSourceID,
	ds.DataSourceCD,
    ds.DataSourceName,
    ds.ObservatoryID,
    ds.SortOrder
FROM
	KS.DataSource AS ds
ORDER BY ds.SortOrder;

END
$$

DELIMITER ;
