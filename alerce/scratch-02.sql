select
oid,
mjd,
ra,
dec
from detections
where oid = 'ZTF18aayhpyh'
order by mjd;