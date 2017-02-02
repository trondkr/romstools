!program  new_astro
subroutine new_astro(dd,astro,ader)
!c	this subroutine calculates the following five ephermides
!c	of the sun and moon
!c	h = mean longitude of the sum
!c	pp = mean longitude of the solar perigee
!c	s = mean longitude of the moon
!c	p = mean longitude of the lunar perigee
!c	np = negative of the longitude of the mean ascending node
!c	and their rates of change.
!c	Units for the ephermides are cycles and for their derivatives
!c	are cycles/365 days
!c	The formulae for calculating this ephermides were taken from
!c	pages 98 and 107 of the Explanatory Supplement to the
!c	Astronomical Ephermeris and the American Ephermis and
!c	Nautical Almanac (1961)

  real :: d1, d,s,h,p,np,pp,dc,tau,dtau
  integer ::sc,hc,pc,npc,ppc
  real,dimension(4) :: args,dargs
  integer :: kd0
  real,dimension(6) :: astro,ader
  real*8  :: dd

!Compute number of days from epoch of 12:00 UT Dec 31, 1899.
  ! (January 0.5 1900 ET)
  call gday(31,12,99,18,kd0,ier)
  if (ier.lt.0) then
     ierr=-1
     return
  end if
  d1=dd-dfloat(kd0)-0.5d0  ! -0.5 skyldes tidpunktet i kommentaren over
  d=d1/10000.0

! Compute astronomical constants at time d1.
  args(1)=d1
  args(2)=d*d
  args(3)=d*d*d


! These are the coefficients of the formulas in the Explan. Suppl.
  s = (270.434164d0 + 13.1763965268d0*args(1) -                   &
       0.000085d0*args(2) + 0.000000039d0*args(3))/360.0d0
  h = (279.696678d0 + 0.9856473354d0*args(1) +                    &
       0.00002267d0*args(2))/360.0d0
  p = (334.329556d0 + 0.1114040803d0*args(1) -                    &
       0.0007739d0*args(2) - 0.00000026d0*args(3))/360.0d0
  np = (-259.183275d0 + 0.0529539222d0*args(1) -                  &
       0.0001557d0*args(2) - 0.00000005d0*args(3))/360.0d0
  pp = (281.220844d0 +  0.0000470684d0*args(1) +                  & 
       0.0000339d0*args(2) + 0.000000070d0*args(3))/360.0d0


  sc=floor(s)
  hc=floor(h)
  pc=floor(p)
  npc=floor(np)
  ppc=floor(pp)
  ! Compute the parameters; we only need the factional part of the cycle.

  astro(2)=s-floor(s)
  astro(3)=h-floor(h)
  astro(4)=p-floor(p)
  astro(5)=np-floor(np)
  astro(6)=pp-floor(pp)

  ! Compute lunar time tau, based on fractional part of solar day.
  ! We add the hour angle to the longitude of the sun and subtract the
  ! longitude of the moon.
  dc=dd-floor(dd)
  tau=dc+astro(3)-astro(2)
  astro(1)=tau
  ! Compute rates of change.
  dargs(1)= 0.0
  dargs(2)= 1.0
  dargs(3)= 2.0E-4*d
  dargs(4)= 3.0E-4*d*d

  ader(2) = (270.434164d0*dargs(1) + 13.1763965268d0*dargs(2) -     &
       0.000085d0*dargs(3) + 0.000000039d0*dargs(4))/360.0d0
  ader(3) = (279.696678d0*dargs(1) + 0.9856473354d0*dargs(2) +      &
       0.00002267d0*dargs(3))/360.0d0
  ader(4) = (334.329556d0*dargs(1) + 0.1114040803d0*dargs(2) -      &
       0.0007739d0*dargs(3) - 0.00000026d0*dargs(4))/360.0d0
  ader(5) = (-259.183275d0*dargs(1) + 0.0529539222d0*dargs(2) -     &
       0.0001557d0*dargs(3) - 0.00000005d0*dargs(4))/360.0d0
  ader(6) = (281.220844d0*dargs(1) +  0.0000470684d0*dargs(2) +     & 
       0.0000339d0*dargs(3) + 0.000000070d0*dargs(4))/360.0d0
  dtau=1.0+ader(3)-ader(1)
  ader(1)=dtau
end subroutine
