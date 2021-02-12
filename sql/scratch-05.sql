ALTER TABLE KS.StateVectors_Planets DROP CONSTRAINT FK_StateVectors_Planets_TimeID;
ALTER TABLE KS.StateVectors_Planets DROP CONSTRAINT FK_StateVectors_Planets_BodyID;
ALTER TABLE KS.StateVectors_Planets engine='aria' transactional=0;
ALTER TABLE KS.StateVectors_Planets ADD CONSTRAINT FK_StateVectors_Planets_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID);
ALTER TABLE KS.StateVectors_Planets ADD CONSTRAINT FK_StateVectors_Planets_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);