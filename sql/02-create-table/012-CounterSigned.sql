-- Create counter table
CREATE OR REPLACE TABLE KS.CounterSigned(
    _ INT NOT NULL PRIMARY KEY
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Utility table - enables a SQL query to emulate a for loop. Data type is signed integer; range is non-negative integers up to 2^24.";

-- Populate counter table up to +/- 2^20 (about 2 million)
INSERT INTO KS.CounterSigned
( _ )
SELECT
	i._
FROM
	KS.Counter AS i
WHERE i._ < power(2,20);

-- Now insert the negative numbers
INSERT INTO KS.CounterSigned
( _ )
SELECT
	-i._
FROM
	KS.CounterSigned AS i
WHERE i._ > 0;