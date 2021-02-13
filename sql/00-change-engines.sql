SET SESSION foreign_key_checks=OFF;

-- Failed
ALTER TABLE KS.Asteroid engine='aria' transactional=1;
-- ALTER TABLE KS.AsteroidElement_Ref engine='aria' transactional=1;
-- ALTER TABLE KS.IntegrationTime engine='aria' transactional=1;
ALTER TABLE KS.Body engine='aria' transactional=1;

-- Success
-- ALTER TABLE KS.DailyTime engine='aria' transactional=1;
-- ALTER TABLE KS.Counter engine='aria' transactional=1;
-- ALTER TABLE KS.Minutes engine='aria' transactional=1;

