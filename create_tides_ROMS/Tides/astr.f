	subroutine astr(d1,h,pp,s,p,np,dh,dpp,ds,dp,dnp)
c	this subroutine calculates the following five ephermides
c	of the sun and moon
c	h = mean longitude of the sum
c	pp = mean longitude of the solar perigee
c	s = mean longitude of the moon
c	p = mean longitude of the lunar perigee
c	np = negative of the longitude of the mean ascending node
c	and their rates of change.
c	Units for the ephermides are cycles and for their derivatives
c	are cycles/365 days
c	The formulae for calculating this ephermides were taken from
c	pages 98 and 107 of the Explanatory Supplement to the
c	Astronomical Ephermeris and the American Ephermis and
c	Nautical Almanac (1961)
c
c	implicit real*8(a-h,o-z)
	real*8      np,d1,h,pp,s,p,dh,dpp,ds,dp,dnp,d2,f,f2
	d2=d1*1.d-4
	f=360.d0
	f2=f/365.d0
	h=279.696678d0+.9856473354d0*d1+.00002267d0*d2*d2
	pp=281.220833d0+.0000470684d0*d1+.0000339d0*d2*d2+
     1  .00000007d0*d2**3
	s=270.434164d0+13.1763965268d0*d1-.000085d0*d2*d2+
     1  .000000039d0*d2**3
	p=334.329556d0+.1114040803d0*d1-.0007739d0*d2*d2-
     1  .00000026d0*d2**3
	np=-259.183275d0+.0529539222d0*d1-.0001557d0*d2*d2-
     1  .00000005d0*d2**3
	h=h/f
	pp=pp/f
	s=s/f
	p=p/f
	np=np/f
	h=h-dint(h)
	pp=pp-dint(pp)
	s=s-dint(s)
	p=p-dint(p)
	np=np-dint(np)
	dh=.9856473354d0+2.d-8*.00002267d0*d1
	dpp=.0000470684d0+2.d-8*.0000339d0*d1
     1  +3.d-12*.00000007d0*d1**2
	ds=13.1763965268d0-2.d-8*.000085d0*d1+
     1  3.d-12*.000000039d0*d1**2
	dp=.1114040803d0-2.d-8*.0007739d0*d1-
     1  3.d-12*.00000026d0*d1**2
	dnp=+.0529539222d0-2.d-8*.0001557d0*d1-
     1  3.d-12*.00000005d0*d1**2
	dh=dh/f2
	dpp=dpp/f2
	ds=ds/f2
	dp=dp/f2
	dnp=dnp/f2
	return
	end
