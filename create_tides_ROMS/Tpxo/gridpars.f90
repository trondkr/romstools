SUBROUTINE gridpars(icall,ldata,idata,igtype,nx,ny,grid,ierror)

!// NAME:
!//    gridpar

!// PURPOSE:
!//    Conversion between integer*2 field identification and
!//    variables used in programs. Note that some gridtypes also
!//    has extended (geometry) identification behind the field data.

!// SYNOPSIS:
!//    SUBROUTINE gridpar(icall,ldata,idata,igtype,nx,ny,grid,ierror)
!//    integer   icall,ldata,igtype,nx,ny,ierror
!//    integer*2 idata(ldata)
!//    realgrid(6)

!// INPUT:
!//    icall  - +1 : from idata to igtype,nx,ny,grid
!//             -1 : from igtype,nx,ny,grid to idata
!//    ldata  - length of the idata array

!// INPUT/OUTPUT:
!//    idata  - field identification, field data (not touched)
!//             and possibly extra geometry specifications
!//             (max 20 words in this version)
!//    igtype - grid type, 1 = polarstereographic grid (true at 60N)
!//                        2 = geographic
!//                        3 = spherical rotated grid
!//                        4 = polarstereographic grid
!//                        5 = mercator grid (unrotated)
!//                        * = unknown grid type (-32767 to +32 accepted)
!//             (note that 'codes' are used in the field identifcation
!//             when extra geometry identification is stored after
!//             the field data, the reason for the +32 limit)
!//    nx     - no. of gridpoints in x direction (1 - 32767)
!//    ny     - no. of gridpoints in y direction (1 - 32767)
!//    grid   - grid parameters

!//    igtype=1,4, polarstereographic grid:
!//      grid(1) = x position of north pole (xp)
!//      grid(2) = y position of north pole (yp)
!//      grid(3) = no. of grid units between North Pole and Equator
!//      grid(4) = grid rotation angle (degrees), positive east, negative west
!//      grid(5) = projection latitude (degrees), standard is 60 (60 deg. N)
!//      grid(6) = 0. (not used)

!//    igtype=2,3, geographic or spherical rotated grid:
!//      grid(1) = western boundary (degrees) (longitude for x=1)
!//      grid(2) = southern boundary (degrees) (latitude for y=1)
!//      grid(3) = longitude increment (degrees)
!//      grid(4) = latitude  increment (degrees)
!//      grid(5) = longitude position of rotated equator (degrees) 
!//                (0 if geographic grid)
!//      grid(6) = latitude position of rotated equator (degrees)
!//                (0 if geographic grid)

!//    igtype=5, mercator (unrotated) grid:
!//      grid(1) = western boundary (degrees) (longitude for x=1)
!//      grid(2) = southern boundary (degrees) (latitude for y=1)
!//      grid(3) = x (longitude) increment (km)
!//      grid(4) = y (latitude)  increment (km)
!//      grid(5) = reference (construction) latitude (degrees)
!//      grid(6) = 0. (not used)

!//    igtype=*, unknown grid type, only use grid type less than 1 if the
!//              grid parameters have no meaning:
!//      grid(1:6) : unknown grid parameters
