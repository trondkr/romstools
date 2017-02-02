subroutine seq_grid(unit,identi, gparam, imax,jmax, lon_seq, lat_seq, angle,topo)

  integer, intent(in), dimension(20) :: identi
  real, intent(in), dimension(6)     :: gparam
  integer, intent(in)                :: imax, jmax, unit
  real, intent(inout),dimension(imax,jmax) :: lon_seq, lat_seq, topo, angle
  !character (len=100), intent(in)    :: seq_grd
  integer :: gridno, igtype, ierr
  real :: xp, yp, dx, ylon
  real :: oy, ox, dy, lon,lat
  real, parameter :: xoff = +1.0
  real, parameter :: yoff = +1.0
  integer :: i, j
  ! --- Mathematical values
  real, parameter :: pi = 3.1415926
  real, parameter :: rad = pi/180.0
  real, parameter :: omega = 2*pi/86400.0
  real, parameter :: Re = 6367.0

! open(12, file=trim(seq_grd),form="unformatted")
!  call getflt(100, unit, 0, identi, imax*jmax, topo_seq, igtype,  gparam, 0, ierr)
 ! print *, imax, jmax

  call getflt(0, unit, 0, identi, imax*jmax, topo, igtype,  gparam, 0, ierr)
  close(unit)
!  print *, "imax",imax,"jmax",jmax


  if (identi(9).ge.1000) then
     gridno = int(identi(9)/1000)
  end if

  if (identi(9).eq.1.or.identi(9).eq.4) then
     ! Polarstereographic grid

     if (identi(17) > 0) then
        ! The standard
        xp = identi(15)*0.01
        yp = identi(16)*0.01
        dx = identi(17)*0.1
        ylon = identi(18)
     else
        ! met.no convention for fine scale grids
        xp = identi(15)
        yp = identi(16)
        dx = -identi(17)*0.1
        ylon = identi(18)
     end if

     ! --- Correct for offset error
     xp = xp + xoff
     yp = yp + yoff
!     print *, "Polarstereographic grid. Grid-parameters: ", xp, yp, dx, ylon

  elseif (identi(9).eq.2.or.gridno.eq.2) then
     ! spherical grid, lat/lon, true north

     if (identi(9).eq.2) then
        ! Only need data from ident
        oy = real(identi(15)/100.0)  ! origos breddegrad
        ox = real(identi(16)/100.0)  ! origos lengdegrad
        dy = real(identi(17)/100.0)  ! dlat
        dx = real(identi(18)/100.0)  ! dlon
     else
        ! Need data from extra ident
        oy = gparam(2)  ! origos breddegrad
        ox = gparam(1)  ! origos lengdegrad
        dy = gparam(4)  ! dlat
        dx = gparam(3)  ! dlon
!        print *, "Geographic grid. Grid-parameters: ", oy, ox, dy, dx
     end if
  else
     print *, 'Input grid on sequential input file is not recognized'
     goto 300

  end if

  if (identi(9).eq.1.or.identi(9).eq.4) then

     do j = 1, jmax
        do i = 1, imax
           call grid2ll(real(i), real(j), lon, lat)
           lon_seq(i,j) = lon
           lat_seq(i,j) = lat
        end do
     end do


     ! Local rotation angle
     do j = 1, jmax
        do i = 1, imax
           angle(i,j) = atan2(xp-i, yp-j)
        end do
     end do


  elseif  (identi(9).eq.2.or.gridno.eq.2) then

     lat=oy
     lon=ox

     do j = 1, jmax
        lat=lat+dy
        lon=ox
        do i = 1, imax
           lon=lon+dx
           lon_seq(i,j) = lon
           lat_seq(i,j) = lat
           angle(i,j) = 0
        end do
     end do
  end if

300 continue
contains

subroutine grid2ll(x,y,lon,lat)
  real, intent(in) :: x, y            ! Grid coordinates
  real, intent(out) :: lon, lat       ! Longitude and latitude
  real, parameter :: deg = 1.0/rad
  real, parameter :: phi0 = 60.0*rad
  real :: r         

  lon = ylon + atan2(x-xp, yp-y)*deg
  r   = dx * sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp))
  lat = 90.0 - 2*atan(r / (Re*(1+sin(phi0))))*deg
end subroutine grid2ll
end subroutine seq_grid


