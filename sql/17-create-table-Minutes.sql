-- Table for number of minutes in a day for integation time
CREATE OR REPLACE TABLE KS.Minutes(
	_ INT NOT NULL PRIMARY KEY,
	mjd_offset double NOT NULL,
	wt0 double NOT NULL,
	wt1 double NOT NULL
);

# Number of minutes in one day
SET @mpd_i = CAST(24*60 AS INT);
SET @mpd_d = CAST(24*60 AS DOUBLE);

INSERT INTO KS.Minutes ( _, mjd_offset, wt0, wt1 )
SELECT 
	minutes._,
	(minutes._ / @mpd_d) AS mjd_offset,
	(@mpd_i - minutes._) / @mpd_d AS wt0,
	minutes._ / @mpd_d AS wt1
FROM KS.Counter AS minutes
WHERE minutes._ < (24*60) AND ((minutes._ % 5) = 0);
