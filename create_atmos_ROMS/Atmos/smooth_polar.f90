SUBROUTINE smooth_polar(ncep_in, ncep_out, nlonp, nlatp)
IMPLICIT NONE
INTEGER, INTENT(IN)                          :: nlonp, nlatp
REAL, DIMENSION(nlonp,nlatp), INTENT(INOUT)  :: ncep_in, ncep_out
REAL sum
INTEGER i,i_lat,latb1,latb2,latt1,latt2

latb1 = 2
latb2 = 10
latt1 = nlatp-9
latt2 = nlatp-1

ncep_out = ncep_in

DO i_lat=latb1,latb2
  ncep_out(2:nlonp-1,i_lat) = 0.25*ncep_in(1:nlonp-2,i_lat)+0.5*ncep_in(2:nlonp-1,i_lat) &
                               +0.25*ncep_in(3:nlonp,i_lat)
  ncep_out(1,i_lat) = 0.25*ncep_in(nlonp-1,i_lat)+0.5*ncep_in(1,i_lat)+0.25*ncep_in(2,i_lat)
  ncep_out(nlonp,i_lat) = ncep_out(1,i_lat)
ENDDO

DO i_lat=latt1,latt2
  ncep_out(2:nlonp-1,i_lat) = 0.25*ncep_in(1:nlonp-2,i_lat)+0.5*ncep_in(2:nlonp-1,i_lat) &
                               +0.25*ncep_in(3:nlonp,i_lat)
  ncep_out(1,i_lat) = 0.25*ncep_in(nlonp-1,i_lat)+0.5*ncep_in(1,i_lat)+0.25*ncep_in(2,i_lat)
  ncep_out(nlonp,i_lat) = ncep_out(1,i_lat)
ENDDO

ncep_in = ncep_out

DO i_lat=latb1,latb2
  ncep_out(2:nlonp-1,i_lat) = 0.25*ncep_in(1:nlonp-2,i_lat)+0.5*ncep_in(2:nlonp-1,i_lat) &
                               +0.25*ncep_in(3:nlonp,i_lat)
  ncep_out(1,i_lat) = 0.25*ncep_in(nlonp-1,i_lat)+0.5*ncep_in(1,i_lat)+0.25*ncep_in(2,i_lat)
  ncep_out(nlonp,i_lat) = ncep_out(1,i_lat)
ENDDO

DO i_lat=latt1,latt2
  ncep_out(2:nlonp-1,i_lat) = 0.25*ncep_in(1:nlonp-2,i_lat)+0.5*ncep_in(2:nlonp-1,i_lat) &
                               +0.25*ncep_in(3:nlonp,i_lat)
  ncep_out(1,i_lat) = 0.25*ncep_in(nlonp-1,i_lat)+0.5*ncep_in(1,i_lat)+0.25*ncep_in(2,i_lat)
  ncep_out(nlonp,i_lat) = ncep_out(1,i_lat)
ENDDO

! Smooth boundary between zonally smoothed and unsmoothed data
ncep_out(:,latb2) = 0.25*ncep_out(:,latb2-1)+0.5*ncep_out(:,latb2)+0.25*ncep_out(:,latb2+1)

ncep_out(:,latt1) = 0.25*ncep_out(:,latt1-1)+0.5*ncep_out(:,latt1)+0.25*ncep_out(:,latt1+1)

! Smooth around the North Pole

sum = 0.
DO i=2,nlonp-1
   sum = sum + ncep_out(i,nlatp-1)
ENDDO
sum = sum/REAL(nlonp-2)
DO i=1,nlonp
   ncep_out(i,nlatp) = sum
ENDDO

END SUBROUTINE smooth_polar
