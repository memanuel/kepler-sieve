-- SHOW processlist;

SELECT
	ps.Id, 
	ps.USER, 
	ps.Host, 
	ps.Command, 
	ps.Time, 
	ps.State, 
	ps.Info	
FROM 
	information_schema.PROCESSLIST AS ps
	WHERE ps.Command <> 'Sleep'	
ORDER BY ps.Id;
