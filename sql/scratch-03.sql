WITH t1 AS (
SELECT
	ad2.AsteroidID,
	ad2.TimeID,
	ad2.LightTime,
	ad2.ux-ad1.ux AS dx,
	ad2.uy-ad1.uy AS dy,
	ad2.uz-ad1.uz AS dz
FROM
	KS.AsteroidDirections2 AS ad2
	INNER JOIN KS.AsteroidDirections AS ad1 ON 
		ad1.AsteroidID = ad2.AsteroidID AND
		ad1.TimeID = ad2.TimeID
WHERE
	ad2.AsteroidID=1
)
SELECT
	t1.AsteroidID,
	t1.TimeID,
	t1.LightTime,
	SQRT(pow(t1.dx, 2)+pow(t1.dy, 2)+pow(t1.dz, 2))*3600*60 AS du_sec
FROM
	t1;