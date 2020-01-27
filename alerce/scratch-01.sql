select
	oid,
	meanra,
	meandec,
	firstmjd,
	nobs,
	--classearly,
	pclassearly
from objects
where 
	classearly=21 and 
	nobs > 2 and
	pclassearly > 0.95
order by nobs desc	
limit 10;