SELECT * FROM KS.StateVectors_Planets LIMIT 20;
SELECT * FROM KS.StateVectors_Planets WHERE TimeID = 59000*24*60;
SELECT * FROM KS.StateVectors_Planets WHERE TimeID > 59000*24*60 GROUP BY TimeID LIMIT 100 ;

