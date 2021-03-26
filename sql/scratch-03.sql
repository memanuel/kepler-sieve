	CHAR(ASCII('W')+t2.j1) AS alpha,


WITH t1 AS(
SELECT
	(i._+1) AS i,
	((i._+1) MOD 3) + 1 AS j1,
	((i._+2) MOD 3) + 1 AS j2,	
	CAST(2*s._ - 1 AS DOUBLE) AS c
FROM
	KS.Counter AS i
	INNER JOIN KS.Counter AS s ON s._ < 2
WHERE
	i._ < 3
)
SELECT
	row_number() OVER (ORDER BY t1.i*t1.c DESC) AS CubeFaceID,
	CONCAT(CHAR(ASCII('X')+t1.i-1), IF(c>0, '+', '-')) AS CubeFaceCD,
	t1.i,
	t1.j1,
	t1.j2,
	t1.c AS ci,
	-- Names of the three axes
	CHAR(ASCII('W')+t1.j1) AS alpha,
	CHAR(ASCII('W')+t1.j2) AS beta,
	CHAR(ASCII('W')+t1.i) AS gamma
FROM
	t1
ORDER BY (t1.i*t1.c) DESC;