select
    obj.oid as ObjectID,
    det.candid as CandidateID,
    det.mjd,
    det.ra,
    det.dec,
    det.magpsf as mag_psf,
    det.magap as mag_app,
    det.magnr as nag_nr,
    det.sigmara as sigma_ra,
    det.sigmadec as sigma_dec,
    obj.pclassearly as asteroid_prob,
    has_stamps,
    parent_candid as ParentCandidateID
from 
    detections as det
    inner join objects as obj on obj.oid = det.oid
where
    obj.classearly = 21