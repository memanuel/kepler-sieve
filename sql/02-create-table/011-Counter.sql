-- Create counter table
CREATE OR REPLACE TABLE KS.Counter(
    _ INT NOT NULL PRIMARY KEY
)
ENGINE='aria' TRANSACTIONAL=1;

-- Populate counter table up to 2^24 (about 4 million)
INSERT INTO KS.Counter ( _ )
VALUES (0), (1), (2), (3);

INSERT IGNORE INTO KS.Counter (_)
SELECT
(power(2,2)*x._ + y._) AS z
FROM
	KS.Counter AS x
	CROSS JOIN KS.Counter AS y;

INSERT IGNORE INTO KS.Counter (_)
SELECT
(power(2,4)*x._ + y._) AS z
FROM
	KS.Counter AS x
	CROSS JOIN KS.Counter AS y;

INSERT IGNORE INTO KS.Counter (_)
SELECT
(power(2,8)*x._ + y._) AS z
FROM
	KS.Counter AS x
	CROSS JOIN KS.Counter AS y;

INSERT IGNORE INTO KS.Counter (_)
SELECT
(power(2,16)*x._ + y._) AS z
FROM
	KS.Counter AS x
	CROSS JOIN KS.Counter AS y
WHERE x._ < power(2,8);