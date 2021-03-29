select
	-- Grid coordinates
-- 	gd.i1,
-- 	gd.j1,
-- 	gd.i2,
-- 	gd.j2,
	-- Change in grid coordinates
	(gd.i2-gd.i1) as di,
	(gd.j2-gd.j1) as dj,
	-- Tabulated distances
	gd.dr_mid*360 as dr_mid,
	gd.dr_min*360 as dr_min,
	-- Point 1 on the cube
	g1.a as a1,
	g1.b as b1,
	g1.c as c1,
	SQRT(POW(g1.a,2)+POW(g1.b,2)+POW(g1.c,2)) as r1,
	-- Point 2 on cube
	g2.a as a2,
	g2.b as b2,
	g2.c as c2,
	SQRT(POW(g2.a,2)+POW(g2.b,2)+POW(g2.c,2)) as r2,
	-- Change in cube coordinates
	(g2.a-g1.a)*360 as da,
	(g2.b-g1.b)*360 as db,
	(g2.c-g1.c)*360 as dc,
	-- Points on the sphere
	g1.u as u1,
	g1.v as v1,
	g1.w as w1,
	g2.u as u2,
	g2.v as v2,
	g2.w as w2,
	-- Change in sphere coordinates
	(g2.u-g1.u)*360 as du,
	(g2.v-g1.v)*360 as dv,
	(g2.w-g1.w)*360 as dw,
	-- Distance on sphere
	SQRT((pow(g2.u-g1.u,2)+pow(g2.v-g1.v,2)+pow(g2.w-g1.w,2)))*360 as dr_mid2
from 
	KS.SkyPatchGridDistance AS gd
	inner join KS.SkyPatchGrid as g1 on g1.i=gd.i1 and g1.j=gd.j1 
	inner join KS.SkyPatchGrid as g2 on g2.i=gd.i2 and g2.j=gd.j2
where gd.i1=1024 and gd.j1=1024 and i2-i1=3 order by j2;