select
min(mjd) as mjd_min,
max(mjd) as mjd_max,
min(ra) as ra_min,
max(ra) as ra_max,
min(dec) as dec_min,
max(dec) as dec_max
from detections
where oid = 'ZTF18aayhpyh'
group by oid;