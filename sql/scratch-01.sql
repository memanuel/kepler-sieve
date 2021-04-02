-- SHOW processlist;

SELECT
	ps.Id, ps.USER, ps.Host, ps.Command, ps.Time, ps.State, ps.Info, ps.Progress
FROM 
	information_schema.PROCESSLIST AS ps
	WHERE ps.Command <> 'Sleep'	
ORDER BY ps.Id;

KILL 1272;