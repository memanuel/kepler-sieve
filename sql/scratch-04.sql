SELECT
	elts.*
FROM
	KS.AsteroidElements AS elts
WHERE elts.AsteroidID=1;	


CALL KS.GetAsteroidElements(1, 2, 48000, 63000);

CALL KS.GetAsteroidData(1, 11, 48000, 63000);