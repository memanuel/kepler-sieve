CREATE OR REPLACE TABLE KS.PrimaryBody(
    BodyID INT NOT NULL PRIMARY KEY
        COMMENT "Reference to the Body table; the body whose default primary is described.",
    PrimaryBodyID INT NOT NULL
        COMMENT "The default primary body of this body. Almost always the Sun except that the primary of the Moon is Earth",
    CONSTRAINT FK_PrimaryBody_BodyID
        FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID),
    CONSTRAINT FK_PrimaryBody_PrimaryBodyID
        FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Default primary for each body when integrated in the Planets or DE435 collections.  Almost always the Sun except that the primary of the Moon is Earth.";

