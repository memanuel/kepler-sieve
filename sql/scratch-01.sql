SELECT 
COUNT(elt.AsteroidID)
FROM KS.AsteroidElement_Ref AS elt
WHERE elt.TimeID < (59000*24*60)