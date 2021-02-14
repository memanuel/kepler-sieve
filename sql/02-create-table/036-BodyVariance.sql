CREATE OR REPLACE TABLE JPL.BodyVariance(
	BodyID INT NOT NULL PRIMARY KEY
		COMMENT "The Body whose state vectors are described; FK to JS.Body",
    var_q DOUBLE NOT NULL
        COMMENT "The variance of the position of this body in AU^2",
    var_v DOUBLE NOT NULL
        COMMENT "The variance of the velocity of this body in (AU/day)^2",       
    sd_q DOUBLE NOT NULL
        COMMENT "The standard deviation of the position of this body in AU",
    sd_v DOUBLE NOT NULL
        COMMENT "The standard deviation of the velocity of this body in AU/day"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Variance and standard deviation of position and velocity of solar system bodies.  Used to compute relative accuracy of an integration compared to Horizons.";
