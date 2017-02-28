! Program to generate ROMS forcing files from ERA40 netCDF files produced using cvt_gribfiles.bsh

PROGRAM get_forcing

use netcdf

implicit none

INTEGER i, j, imon, itime, ntimes, ia, ja, IFL
INTEGER nlons, nlonp, nlats ,LonDimID, LatDimID, TimeDimId
INTEGER XiDimID, EtaDimID, Lp, Mp
INTEGER MaskVarId, LonVarId, LatVarId, DiVarId, DjVarId
INTEGER LonRhoVarId, LatRhoVarId, AngleVarId
INTEGER U10VarId, V10VarId, TimeANVarId, TimeFCVarId
INTEGER MSLVarId
INTEGER T2MVarId
INTEGER D2MVarId
INTEGER TCCVarId
INTEGER CPVarId, LSPVarId, SSRDVarId, STRDVarId
INTEGER Iti_AN, ItI_FC

INTEGER XiRhoDimId, EtaRhoDimId, TimeOutDimID
INTEGER ncid_Uwind, UwindVarId, TimeUwindVarId
INTEGER ncid_Vwind, VwindVarId, TimeVwindVarId
INTEGER ncid_Pair, PairVarId, TimePairVarId
INTEGER ncid_Tair, TairVarId, TimeTairVarId
INTEGER ncid_Qair, QairVarId, TimeQairVarId
INTEGER ncid_cloud, cloudVarId, TimecloudVarId
INTEGER ncid_rain, rainVarId, TimerainVarId
INTEGER ncid_swrad, swradVarId, TimeswradVarId
INTEGER ncid_lwrad, lwradVarId, TimelwradVarId
INTEGER currentYear

INTEGER status_era, ncid_era
INTEGER status_roms, ncid_roms
INTEGER status_AN, ncid_AN
INTEGER status_AN_ml, ncid_AN_ml
INTEGER status_FC, ncid_FC

INTEGER ndays, n6hr, n12hr, datetime

INTEGER i1, i2, j1, j2
INTEGER jd, jdref, jdref_era
REAL  Di, Dj, rx, rxm, ry, rym
REAL rimin, rimax, rjmin, rjmax
REAL, DIMENSION(:,:), allocatable :: ipos
REAL, DIMENSION(:,:), allocatable :: jpos
REAL delt

REAL U10_scale_factor, U10_add_offset
REAL V10_scale_factor, V10_add_offset
REAL MSL_scale_factor, MSL_add_offset
REAL T2M_scale_factor, T2M_add_offset
REAL D2M_scale_factor, D2M_add_offset
REAL TCC_scale_factor, TCC_add_offset

REAL CP_scale_factor, CP_add_offset
REAL LSP_scale_factor, LSP_add_offset
REAL SSRD_scale_factor, SSRD_add_offset
REAL STRD_scale_factor, STRD_add_offset

INTEGER nvalue
REAL tx, critx, cor
INTEGER mxs

REAL, ALLOCATABLE, DIMENSION(:) :: lons
REAL, ALLOCATABLE, DIMENSION(:) :: lonp
REAL, ALLOCATABLE, DIMENSION(:) :: lats
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr
REAL, ALLOCATABLE, DIMENSION(:,:) :: scrp
REAL, ALLOCATABLE, DIMENSION(:,:) :: scru
REAL, ALLOCATABLE, DIMENSION(:,:) :: scrv
REAL, ALLOCATABLE, DIMENSION(:,:) :: emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: work
REAL, ALLOCATABLE, DIMENSION(:,:) :: error

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr1_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr2_out

REAL, ALLOCATABLE, DIMENSION(:,:) :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: angle

REAL, ALLOCATABLE, DIMENSION(:,:) :: U10
REAL, ALLOCATABLE, DIMENSION(:,:) :: V10
REAL, ALLOCATABLE, DIMENSION(:,:) :: MSL
REAL, ALLOCATABLE, DIMENSION(:,:) :: T2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: D2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: Q2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: evapres
REAL, ALLOCATABLE, DIMENSION(:,:) :: TCC
REAL, ALLOCATABLE, DIMENSION(:,:) :: TP
REAL, ALLOCATABLE, DIMENSION(:,:) :: CP
REAL, ALLOCATABLE, DIMENSION(:,:) :: LSP
REAL, ALLOCATABLE, DIMENSION(:,:) :: SSRD
REAL, ALLOCATABLE, DIMENSION(:,:) :: STRD

REAL, ALLOCATABLE, DIMENSION(:,:) :: Uwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Vwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Pair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Tair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Qair
REAL, ALLOCATABLE, DIMENSION(:,:) :: cloud
REAL, ALLOCATABLE, DIMENSION(:,:) :: rain
REAL, ALLOCATABLE, DIMENSION(:,:) :: swrad
REAL, ALLOCATABLE, DIMENSION(:,:) :: lwrad
REAL tday

REAL, PARAMETER    :: undef = 2.E+35            ! Undefined land value
REAL pi, DTOR

CHARACTER(len=80) time_dimname
CHARACTER(len=80) maskfile
CHARACTER(len=80) AN_file
CHARACTER(len=80) FC_file
CHARACTER(len=80) gridfile, path_out
CHARACTER(len=80) file_Uwind, file_Vwind
CHARACTER(len=80) file_Pair
CHARACTER(len=80) file_Tair
CHARACTER(len=80) file_Qair
CHARACTER(len=80) file_cloud
CHARACTER(len=80) file_rain
CHARACTER(len=80) file_swrad
CHARACTER(len=80) file_lwrad
CHARACTER(len=80) currentYearString

INTERFACE
  SUBROUTINE spec_hum(td, ta, p, evapres, q, IFL)
  IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(IN) :: td 
    REAL, DIMENSION(:,:), INTENT(IN) :: ta 
    REAL, DIMENSION(:,:), INTENT(IN) :: p(:,:)
    INTEGER, INTENT(IN) :: IFL  
    REAL, DIMENSION(:,:), INTENT(OUT) :: evapres
    REAL, DIMENSION(:,:), INTENT(OUT) ::  q
    REAL, DIMENSION (SIZE(td, 1),SIZE(td, 2)) :: arg
    REAL, DIMENSION (SIZE(td, 1),SIZE(td, 2)) :: f
  END SUBROUTINE spec_hum
END INTERFACE

! Get reference Julian day (1948-01-01 00:00:00Z) for output
jdref = jd(1948,1,1)
! Get reference Julian day (1900-01-01 00:00:00Z) for input
jdref_era = jd(1900,1,1)

pi = ATAN(1.)*4.
DTOR = pi/180.
tx = 0.9*undef
critx = 0.01
cor = 1.6
mxs = 100

!// Read standard input
read(5,'(a)') gridfile  ! Grid file for input grid
read(5,'(a)') AN_file   ! Name of ERA file with AN fields
read(5,'(a)') FC_file   ! Name of ERA file with FC fields (rain and fluxes)
read(5,'(a)') maskfile  ! File with ERA land/sea-mask
read(5,'(a)') path_out  ! Directory where output files will be written
read(5,'(a)') currentYearString
write(*,*) 'Grid file = ',TRIM(gridfile)

read(currentYearString,'(i)') currentYear

! Get info on ROMS grid
status_roms = nf90_open(TRIM(gridfile),nf90_nowrite,ncid_roms)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_dimid(ncid_roms, "xi_rho", XiDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_dimid(ncid_roms, "eta_rho", EtaDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_Inquire_Dimension(ncid_roms,XiDimID,len = Lp)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_Inquire_Dimension(ncid_roms,EtaDimID,len = Mp)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
ALLOCATE(lon_rho(Lp,Mp))
ALLOCATE(lat_rho(Lp,Mp))
ALLOCATE(angle(Lp,Mp))
ALLOCATE(ipos(Lp,Mp))
ALLOCATE(jpos(Lp,Mp))
ALLOCATE(scr1_out(Lp,Mp))
ALLOCATE(scr2_out(Lp,Mp))

write(*,*) 'ROMS Lp, Mp = ',Lp,Mp

! Get lons and lats
status_roms = nf90_inq_varid(ncid_roms, "lon_rho",LonRhoVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "lat_rho",LatRhoVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "angle",AngleVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,LonRhoVarId,lon_rho)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,LatRhoVarId,lat_rho)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,AngleVarId,angle)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_roms = nf90_close(ncid_roms)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! Get dimensions of ERA40 files
status_era = nf90_open(TRIM(maskfile),nf90_nowrite,ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_dimid(ncid_era, "longitude", LonDimID)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_dimid(ncid_era, "latitude", LatDimID)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_Inquire_Dimension(ncid_era,LonDimID,len = nlons)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_Inquire_Dimension(ncid_era,LatDimID,len = nlats)
if(status_era /= nf90_NoErr) call handle_err(status_era)

nlonp = nlons+2

ALLOCATE(lons(nlons))
ALLOCATE(lonp(nlonp))
ALLOCATE(lats(nlats))
ALLOCATE(scr(nlons,nlats))
ALLOCATE(scrp(nlonp,nlats))
ALLOCATE(scru(nlonp,nlats))
ALLOCATE(scrv(nlonp,nlats))
ALLOCATE(emask(nlonp,nlats))
ALLOCATE(work(nlonp,nlats))
ALLOCATE(error(nlonp,nlats))

ALLOCATE(U10(nlonp,nlats))
ALLOCATE(V10(nlonp,nlats))
ALLOCATE(MSL(nlonp,nlats))
ALLOCATE(T2M(nlonp,nlats))
ALLOCATE(D2M(nlonp,nlats))
ALLOCATE(Q2M(nlonp,nlats))
ALLOCATE(evapres(nlonp,nlats))
ALLOCATE(TCC(nlonp,nlats))
ALLOCATE(TP(nlonp,nlats))
ALLOCATE(CP(nlonp,nlats))
ALLOCATE(LSP(nlonp,nlats))
ALLOCATE(SSRD(nlonp,nlats))
ALLOCATE(STRD(nlonp,nlats))

! Get lons and lats
status_era = nf90_inq_varid(ncid_era, "longitude",LonVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_varid(ncid_era, "latitude",LatVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,LonVarId,lons)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,LatVarId,lats)
if(status_era /= nf90_NoErr) call handle_err(status_era)
lonp(2:nlonp-1) = lons
lonp(1) = lons(1)-(lons(2)-lons(1))
lonp(nlonp) = lons(nlons)+(lons(nlons)-lons(nlons-1))

!status_era = nf90_inq_varid(ncid_era, "Di",DiVarId)
!if(status_era /= nf90_NoErr) call handle_err(status_era)
!status_era = nf90_inq_varid(ncid_era, "Dj",DjVarId)
!if(status_era /= nf90_NoErr) call handle_err(status_era)
!status_era = nf90_get_var(ncid_era,DiVarId,Di)
!if(status_era /= nf90_NoErr) call handle_err(status_era)
!status_era = nf90_get_var(ncid_era,DjVarId,Dj)
!if(status_era /= nf90_NoErr) call handle_err(status_era)

Di = 0.7
Dj = 0.7

! Get land-sea mask
status_era = nf90_inq_varid(ncid_era, "lsm",MaskVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,MaskVarId,scr)
if(status_era /= nf90_NoErr) call handle_err(status_era)
emask(2:nlons+1,:) = scr
emask(1,:) = scr(nlons,:)
emask(nlonp,:) = scr(1,:)

write(*,*) 'nlons,nlats = ',nlons,nlats

status_era = nf90_close(ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)


! Get i,j positions of ROMS rho-points on ERA40 grid
WHERE(lon_rho< 0.) lon_rho = lon_rho+360.
DO j=1,Mp
DO i=1,Lp
   ipos(i,j) = (lon_rho(i,j)-lonp(1))/Di + 1.0
   jpos(i,j) = (lats(1)-lat_rho(i,j))/Dj + 1.0
ENDDO
ENDDO
jpos(i,j) = MIN(jpos(i,j),FLOAT(nlats))
! Find subarea in ERA40 grid containing model grid
rimin = 1.e+23
rimax =-1.e+23
rjmin = 1.e+23
rjmax =-1.e+23
DO j=1,Mp
DO i=1,Lp
   rimin = MIN(rimin,ipos(i,j))
   rimax = MAX(rimax,ipos(i,j))
   rjmin = MIN(rjmin,jpos(i,j))
   rjmax = MAX(rjmax,jpos(i,j))
ENDDO
ENDDO
i1 = FLOOR(rimin)
i2 = CEILING(rimax)
j1 = FLOOR(rjmin)
j2 = CEILING(rjmax)

IF(i1<1.OR.i2>nlonp.OR.j1<1.OR.j2>nlats) THEN
   write(*,*) 'i1,i2,j1,j2 = ',i1,i2,j1,j2
   write(*,*) 'nlonp, nlats = ',nlonp,nlats
   STOP 'subdomain outside ERA40 domain'
ENDIF

! Winds
! -------------------------------------------------------------------

! Open output files for winds
! ...
file_Uwind = TRIM(path_out) // 'Uwind.nc'
!
status_roms = nf90_create(TRIM(file_Uwind),nf90_clobber,ncid_Uwind)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Uwind,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Uwind,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Uwind,"wind_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind,nf90_global,"title", &
                  "Xi-component of wind forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Uwind,"Uwind",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), UwindVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind, UwindVarId,"long_Name", &
                 "Xi-component of wind")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind, UwindVarId,"units", &
                 "meter second-1")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind, UwindVarId,"time", &
                 "wind_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Uwind,"wind_time",nf90_float, &
                 TimeOutDimID, TimeUwindVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind, TimeUwindVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Uwind, TimeUwindVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_Uwind)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...

file_Vwind = TRIM(path_out) // 'Vwind.nc'
!
status_roms = nf90_create(TRIM(file_Vwind),nf90_clobber,ncid_Vwind)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Vwind,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Vwind,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Vwind,"wind_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind,nf90_global,"title", &
                  "Eta-component of wind forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Vwind,"Vwind",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), VwindVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind, VwindVarId,"long_Name", &
                 "Eta-component of wind")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind, VwindVarId,"units", &
                 "meter second-1")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind, VwindVarId,"time", &
                 "wind_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Vwind,"wind_time",nf90_float, &
                 TimeOutDimID, TimeVwindVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind, TimeVwindVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Vwind, TimeVwindVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_Vwind)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(Uwind(Lp,Mp))
ALLOCATE(Vwind(Lp,Mp))
! -------------------------------------------------------------------

! Atmospheric pressure
! -------------------------------------------------------------------
file_Pair = TRIM(path_out) // 'Pair.nc'
!
status_roms = nf90_create(TRIM(file_Pair),nf90_clobber,ncid_Pair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Pair,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Pair,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Pair,"pair_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair,nf90_global,"title", &
                  "Atmospheric pressure forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Pair,"Pair",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), PairVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair, PairVarId,"long_Name", &
                 "Atmospheric pressure")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair, PairVarId,"units", &
                 "Pa")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair, PairVarId,"time", &
                 "pair_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Pair,"pair_time",nf90_float, &
                 TimeOutDimID, TimePairVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair, TimePairVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Pair, TimePairVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_Pair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(Pair(Lp,Mp))

! ---------------------------------------------------------------

! Air temperature at 2m
! -------------------------------------------------------------------
file_Tair = TRIM(path_out) // 'Tair.nc'
!
status_roms = nf90_create(TRIM(file_Tair),nf90_clobber,ncid_Tair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Tair,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Tair,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Tair,"Tair_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair,nf90_global,"title", &
                  "Air temperature forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Tair,"Tair",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), TairVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair, TairVarId,"long_Name", &
                 "Air temperature at 2m")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair, TairVarId,"units", &
                 "degree C")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair, TairVarId,"time", &
                 "Tair_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Tair,"Tair_time",nf90_float, &
                 TimeOutDimID, TimeTairVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair, TimeTairVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Tair, TimeTairVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_Tair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(Tair(Lp,Mp))


! ---------------------------------------------------------------

! Specific humidity at 2m
! -------------------------------------------------------------------
file_Qair = TRIM(path_out) // 'Qair.nc'
!
status_roms = nf90_create(TRIM(file_Qair),nf90_clobber,ncid_Qair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Qair,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Qair,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_Qair,"Qair_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair,nf90_global,"title", &
                  "Specific humidity forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Qair,"Qair",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), QairVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair, QairVarId,"long_Name", &
                 "Specific humidity at 2m")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair, QairVarId,"units", &
                 "kg kg-1")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair, QairVarId,"time", &
                 "Qair_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_Qair,"Qair_time",nf90_float, &
                 TimeOutDimID, TimeQairVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair, TimeQairVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_Qair, TimeQairVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_Qair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(Qair(Lp,Mp))

! ---------------------------------------------------------------
! Fraction total cloud cover
! -------------------------------------------------------------------
file_cloud = TRIM(path_out) // 'cloud.nc'
!
status_roms = nf90_create(TRIM(file_cloud),nf90_clobber,ncid_cloud)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_cloud,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_cloud,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_cloud,"cloud_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud,nf90_global,"title", &
                  "Cloud fraction forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_cloud,"cloud",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), cloudVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud, cloudVarId,"long_Name", &
                 "Fraction cloud cover")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud, cloudVarId,"units", &
                 " ")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud, cloudVarId,"time", &
                 "cloud_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_cloud,"cloud_time",nf90_float, &
                 TimeOutDimID, TimecloudVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud, TimecloudVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_cloud, TimecloudVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_cloud)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(cloud(Lp,Mp))

! ---------------------------------------------------------------
! Total precipitation
! -------------------------------------------------------------------
file_rain = TRIM(path_out) // 'rain.nc'
!
status_roms = nf90_create(TRIM(file_rain),nf90_clobber,ncid_rain)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_rain,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_rain,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_rain,"rain_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain,nf90_global,"title", &
                  "Total precipitation forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_rain,"rain",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), rainVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain, rainVarId,"long_Name", &
                 "Total precipiation")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain, rainVarId,"units", &
                 "kg meter-2 second-1")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain, rainVarId,"time", &
                 "rain_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_rain,"rain_time",nf90_float, &
                 TimeOutDimID, TimerainVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain, TimerainVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_rain, TimerainVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_rain)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(rain(Lp,Mp))

! ---------------------------------------------------------------
! Downward shortwave radiation at the surface
! -------------------------------------------------------------------
file_swrad = TRIM(path_out) // 'swrad.nc'
!
status_roms = nf90_create(TRIM(file_swrad),nf90_clobber,ncid_swrad)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_swrad,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_swrad,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_swrad,"swrad_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad,nf90_global,"title", &
                  "Downward shortwave radiation forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_swrad,"swrad",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), swradVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad, swradVarId,"long_Name", &
                 "Downward shortwave radiation")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad, swradVarId,"units", &
                 "Watt meter-2")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad, swradVarId,"time", &
                 "swrad_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_swrad,"swrad_time",nf90_float, &
                 TimeOutDimID, TimeswradVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad, TimeswradVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_swrad, TimeswradVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_swrad)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(swrad(Lp,Mp))

! ---------------------------------------------------------------
! Downward longwave radiation at the surface
! -------------------------------------------------------------------
file_lwrad = TRIM(path_out) // 'lwrad.nc'
!
status_roms = nf90_create(TRIM(file_lwrad),nf90_clobber,ncid_lwrad)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_lwrad,"xi_rho",Lp,XiRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_lwrad,"eta_rho",Mp,EtaRhoDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_dim(ncid_lwrad,"lwrad_time",nf90_unlimited,TimeOutDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad,nf90_global,"title", &
                  "Downward longwave radiation forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad,nf90_global,"type", &
                  "ROMS forcing file")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad,nf90_global,"history", &
                  "Produced with get_forcing.f90")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_lwrad,"lwrad_down",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), lwradVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad, lwradVarId,"long_Name", &
                 "Downward longwave radiation")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad, lwradVarId,"units", &
                 "Watt meter-2")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad, lwradVarId,"time", &
                 "lwrad_time")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_def_var(ncid_lwrad,"lwrad_time",nf90_float, &
                 TimeOutDimID, TimelwradVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad, TimelwradVarId,"long_Name", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_put_att(ncid_lwrad, TimelwradVarId,"units", &
                 "days since 1948-1-1 00:00Z")
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_enddef(ncid_lwrad)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
! ...
 
ALLOCATE(lwrad(Lp,Mp))

! ---------------------------------------------------------------
! ---------------------------------------------------------------

!Do 6-hourly analyzed fields

! Open file of analyzed 6-hour fields

  write(*,'(a)') AN_file
  status_AN = nf90_open(TRIM(AN_file),nf90_nowrite,ncid_AN)
  if(status_AN /= nf90_NoErr) call handle_err(status_AN)

! Get time dimension
  status_AN = nf90_inq_dimid(ncid_AN, "time",TimeDimId)
  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_Inquire_Dimension(ncid_AN,TimeDimId,time_dimname,Iti_AN)
  if(status_AN /= nf90_NoErr) call handle_err(status_AN)

! Get time var id
  status_AN = nf90_inq_varid(ncid_AN, "time",TimeANVarId)
  if(status_AN /= nf90_NoErr) call handle_err(status_AN)

! Get U10 var id

 if ( currentYear <= 1988 .OR. currentYear .EQ. 2013) then
   print*,"Before 1988 ",currentYear
   status_AN = nf90_inq_varid(ncid_AN, "v10u",U10VarId)
 else if ( currentYear > 1988 ) then
  print*,"After 1988 ",currentYear
   status_AN = nf90_inq_varid(ncid_AN, "u10",U10VarId)
 end if

  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, U10VarID, "scale_factor", U10_scale_factor)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, U10VarID, "add_offset", U10_add_offset)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  
! Get V10 var id
  if ( currentYear <= 1988 .OR. currentYear .EQ. 2013) then
   status_AN = nf90_inq_varid(ncid_AN, "v10v",V10VarId)
 else if ( currentYear > 1988 ) then
   status_AN = nf90_inq_varid(ncid_AN, "v10",V10VarId)
 end if

  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, V10VarID, "scale_factor", V10_scale_factor)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, V10VarID, "add_offset", V10_add_offset)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)


! Atmospheric pressure forcing
! _________________________________________________________________

! Get MSL var id
  status_AN = nf90_inq_varid(ncid_AN, "msl",MSLVarId)
  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, MSLVarID, "scale_factor", MSL_scale_factor)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, MSLVarID, "add_offset", MSL_add_offset)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)


! Air temperature at 2m
! _________________________________________________________________

! Get T2M var id
  if ( currentYear <= 1988 .OR. currentYear .EQ. 2013) then
   status_AN = nf90_inq_varid(ncid_AN, "v2t",T2MVarId)
 else if ( currentYear > 1988 ) then
   status_AN = nf90_inq_varid(ncid_AN, "t2m",T2MVarId)
 end if

  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, T2MVarID, "scale_factor", T2M_scale_factor)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, T2MVarID, "add_offset", T2M_add_offset)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)


! Dew point temperature at 2m
! _________________________________________________________________

! Get D2M var id
 if ( currentYear <= 1988 .OR. currentYear .EQ. 2013) then
   status_AN = nf90_inq_varid(ncid_AN, "v2d",D2MVarId)
 else if ( currentYear > 1988 ) then
   status_AN = nf90_inq_varid(ncid_AN, "d2m",D2MVarId)
 end if

  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, D2MVarID, "scale_factor", D2M_scale_factor)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, D2MVarID, "add_offset", D2M_add_offset)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)


! Fraction total cloud cover
! _________________________________________________________________

! Get TCC var id
  status_AN = nf90_inq_varid(ncid_AN, "tcc",TCCVarId)
  if(status_AN /= nf90_NoErr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, TCCVarID, "scale_factor", TCC_scale_factor)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)
  status_AN = nf90_get_att(ncid_AN, TCCVarID, "add_offset", TCC_add_offset)
  if (status_AN /= nf90_noerr) call handle_err(status_AN)

!_____________________________________________________________________
!_____________________________________________________________________

!Do 12-hour integrated forecast fields

! Open file of 12-hour fields

  write(*,'(a)') FC_file
  status_FC = nf90_open(TRIM(FC_file),nf90_nowrite,ncid_FC)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)

! Get time dimension
  status_FC = nf90_inq_dimid(ncid_FC, "time",TimeDimId)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)
  status_FC = nf90_Inquire_Dimension(ncid_FC,TimeDimId,time_dimname,Iti_FC)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)

! Get time var id
  status_FC = nf90_inq_varid(ncid_FC, "time",TimeFCVarId)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)

! Get CP var id
  status_FC = nf90_inq_varid(ncid_FC, "cp",CPVarId)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, CPVarID, "scale_factor", CP_scale_factor)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, CPVarID, "add_offset", CP_add_offset)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)

! Get LSP var id
  status_FC = nf90_inq_varid(ncid_FC, "lsp",LSPVarId)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, LSPVarID, "scale_factor", LSP_scale_factor)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, LSPVarID, "add_offset", LSP_add_offset)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)

! Get SSRD var id
  status_FC = nf90_inq_varid(ncid_FC, "ssrd",SSRDVarId)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, SSRDVarID, "scale_factor", SSRD_scale_factor)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, SSRDVarID, "add_offset", SSRD_add_offset)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)

! Get STRD var id
  status_FC = nf90_inq_varid(ncid_FC, "strd",STRDVarId)
  if(status_FC /= nf90_NoErr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, STRDVarID, "scale_factor", STRD_scale_factor)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)
  status_FC = nf90_get_att(ncid_FC, STRDVarID, "add_offset", STRD_add_offset)
  if (status_FC /= nf90_noerr) call handle_err(status_FC)

! _________________________________________________________________
! _________________________________________________________________

! Number of 6- and 12-hourly records to process

n6hr = Iti_AN
n12hr = Iti_FC

! ..................................
! Start time frame loop

  DO itime = 1,n6hr

! Wind forcing
! _________________________________________________________________

! Read in U10 and V10 values
    status_era = nf90_get_var(ncid_AN,U10VarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    U10(2:nlons+1,:) = scr
    U10(nlonp,:) = scr(1,:)
    U10(1,:) = scr(nlons,:)
    U10 = U10*U10_scale_factor + U10_add_offset

    status_era = nf90_get_var(ncid_AN,V10VarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    V10(2:nlons+1,:) = scr
    V10(nlonp,:) = scr(1,:)
    V10(1,:) = scr(nlons,:)
    V10 = V10*V10_scale_factor + V10_add_offset

! Rotate components to polar stereographic
    DO i=1,nlonp
      scru(i,:) = U10(i,:)*COS(DTOR*lonp(i)) - V10(i,:)*SIN(DTOR*lonp(i))
      scrv(i,:) = V10(i,:)*COS(DTOR*lonp(i)) + U10(i,:)*SIN(DTOR*lonp(i))
    ENDDO

! Read in time value
    status_era = nf90_get_var(ncid_AN,TimeANVarId,datetime, &
                            start = (/itime/))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    tday = REAL(jdref_era) + REAL(datetime)/24. - REAL(jdref)
!    write(*,*) 'tday = ',tday
! ...
! U-comp
! Fill in masked-out values 
    scrp = 0.
    scrp = scru
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work,error,nvalue)

! Smooth zonally
    CALL smooth_polar(scrp, scru, nlonp, nlats)

! V-comp
! Fill in masked-out values 
    scrp = 0.
    scrp = scrv
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
! Smooth zonally
    CALL smooth_polar(scrp, scrv, nlonp, nlats)

    DO i=1,nlonp
      U10(i,:) = scru(i,:)*COS(DTOR*lonp(i)) + scrv(i,:)*SIN(DTOR*lonp(i))
      V10(i,:) = scrv(i,:)*COS(DTOR*lonp(i)) - scru(i,:)*SIN(DTOR*lonp(i))
    ENDDO

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*U10(ia,ja)           +   &
                      rx*rym*U10(ia+1,ja)          +   &
                      rx*ry*U10(ia+1,ja+1)         +   &
                      rxm*ry*U10(ia,ja+1)
   ENDDO
   ENDDO


! Horizontal interpolation
   scr2_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr2_out(i,j) = rxm*rym*V10(ia,ja)           +   &
                      rx*rym*V10(ia+1,ja)          +   &
                      rx*ry*V10(ia+1,ja+1)         +   &
                      rxm*ry*V10(ia,ja+1)
   ENDDO
   ENDDO

! Rotate wind components from E/W and N/S to Xi and Eta directions

   Uwind = scr1_out*cos(angle) + scr2_out*sin(angle)
   Vwind = scr2_out*cos(angle) - scr1_out*sin(angle)

! Blank out Antarctica
  WHERE(ABS(Uwind)>1.e+7)
     Uwind = 0.
  ENDWHERE
  WHERE(ABS(Vwind)>1.e+7)
     Vwind = 0.
  ENDWHERE

! Write out wind components in netCDF file

   status_roms = nf90_put_var(ncid_Uwind,UwindVarId,Uwind,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

   status_roms = nf90_put_var(ncid_Vwind,VwindVarId,Vwind,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_Uwind,TimeUwindVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)
   status_roms = nf90_put_var(ncid_Vwind,TimeVwindVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! _________________________________________________________________
!
! Atmospheric pressure
! _________________________________________________________________

! Read in MSL values
    status_era = nf90_get_var(ncid_AN,MSLVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    MSL(2:nlons+1,:) = scr
    MSL(nlonp,:) = scr(1,:)
    MSL(1,:) = scr(nlons,:)
    MSL = MSL*MSL_scale_factor + MSL_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = MSL
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    MSL = scrp
! Smooth zonally
    CALL smooth_polar(scrp, MSL, nlonp, nlats)

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*MSL(ia,ja)           +   &
                      rx*rym*MSL(ia+1,ja)          +   &
                      rx*ry*MSL(ia+1,ja+1)         +   &
                      rxm*ry*MSL(ia,ja+1)
   ENDDO
   ENDDO

   Pair = scr1_out

! Blank out Antarctica
  WHERE(ABS(Pair)>1.e+7)
     Pair = 1.e+5
  ENDWHERE

! Write out atmospheric pressure in netCDF file

   status_roms = nf90_put_var(ncid_Pair,PairVarId,Pair,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_Pair,TimePairVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! _________________________________________________________________
!
! Atmospheric temperature at 2m
! _________________________________________________________________

! Read in T2M values
    status_era = nf90_get_var(ncid_AN,T2MVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    T2M(2:nlons+1,:) = scr
    T2M(nlonp,:) = scr(1,:)
    T2M(1,:) = scr(nlons,:)
    T2M = T2M*T2M_scale_factor + T2M_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = T2M
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    T2M = scrp
! Smooth zonally
    CALL smooth_polar(scrp, T2M, nlonp, nlats)

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*T2M(ia,ja)           +   &
                      rx*rym*T2M(ia+1,ja)          +   &
                      rx*ry*T2M(ia+1,ja+1)         +   &
                      rxm*ry*T2M(ia,ja+1)
   ENDDO
   ENDDO

! Convert from Kelvin to Celsius
   Tair = scr1_out - 273.16

! Blank out Antarctica
  WHERE(ABS(Tair)>1.e+7)
     Tair = 0.
  ENDWHERE

! Write out air temperature in netCDF file

   status_roms = nf90_put_var(ncid_Tair,TairVarId,Tair,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_Tair,TimeTairVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! _________________________________________________________________
!
! Dew point temperature at 2m -> specific humidity
! _________________________________________________________________

! Read in D2M values
    status_era = nf90_get_var(ncid_AN,D2MVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    D2M(2:nlons+1,:) = scr
    D2M(nlonp,:) = scr(1,:)
    D2M(1,:) = scr(nlons,:)
    D2M = D2M*D2M_scale_factor + D2M_add_offset

! Convert from Kelvin to Celsius
    D2M = D2M - 273.16

! Compute specific humidity on ECMWF points
    scrp = 0.
    IFL = 1 !Input is dew point temperature
    CALL spec_hum(D2M, T2M, MSL, evapres, Q2M, IFL)

! Fill in masked-out values 
    scrp = 0.
    scrp = Q2M
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    Q2M = scrp
! Smooth zonally
    CALL smooth_polar(scrp, Q2M, nlonp, nlats)

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*Q2M(ia,ja)           +   &
                      rx*rym*Q2M(ia+1,ja)          +   &
                      rx*ry*Q2M(ia+1,ja+1)         +   &
                      rxm*ry*Q2M(ia,ja+1)
   ENDDO
   ENDDO


   Qair = scr1_out

! Blank out Antarctica
  WHERE(ABS(Qair)>1.e+7)
     Qair = 1.e-5
  ENDWHERE

! Write out specific humidity in netCDF file

   status_roms = nf90_put_var(ncid_Qair,QairVarId,Qair,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_Qair,TimeQairVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! _________________________________________________________________
!
! Fraction total cloud cover
! _________________________________________________________________

! Read in TCC values
    status_era = nf90_get_var(ncid_AN,TCCVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    TCC(2:nlons+1,:) = scr
    TCC(nlonp,:) = scr(1,:)
    TCC(1,:) = scr(nlons,:)
    TCC = TCC*TCC_scale_factor + TCC_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = TCC
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    TCC = scrp
! Smooth zonally
    CALL smooth_polar(scrp, TCC, nlonp, nlats)

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*TCC(ia,ja)           +   &
                      rx*rym*TCC(ia+1,ja)          +   &
                      rx*ry*TCC(ia+1,ja+1)         +   &
                      rxm*ry*TCC(ia,ja+1)
   ENDDO
   ENDDO

   cloud = scr1_out

! Blank out Antarctica
  WHERE(ABS(cloud)>1.e+7)
     cloud = 0.
  ENDWHERE

! Write out cloud fraction in netCDF file

   status_roms = nf90_put_var(ncid_cloud,cloudVarId,cloud,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_cloud,TimecloudVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


!_________________________________________________________________
!End of 6-hourly data

ENDDO

! End time frame loop
! ..................................
  status_era = nf90_close(ncid_AN)
   if(status_era /= nf90_NoErr) call handle_err(status_era)

! End 6-hourly (AN) data loop
! --------------------------------------------------------------

status_roms = nf90_close(ncid_Uwind)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_Vwind)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_Pair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_Tair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_Qair)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_cloud)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

!--------------------------------------------------------------------

! 12-hour mean data

 delt = 12.*3600.

DO itime=1,n12hr

! Read in time value
    status_era = nf90_get_var(ncid_FC,TimeFCVarId,datetime, &
                            start = (/itime/))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    tday = REAL(jdref_era) + REAL(datetime)/24. - REAL(jdref) - 0.5  !0.25

! _________________________________________________________________
!
! Convective precipitation
! _________________________________________________________________

! Read in CP values
    status_ERA = nf90_get_var(ncid_FC,CPVarId,scr,start=(/1, 1,itime /))
    if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
    CP(2:nlons+1,:) = scr
    CP(nlonp,:) = scr(1,:)
    CP(1,:) = scr(nlons,:)
    CP = CP*CP_scale_factor + CP_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = CP
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    CP = scrp
! Smooth zonally
    CALL smooth_polar(scrp, CP, nlonp, nlats)

! _________________________________________________________________
!
! Large-Scale precipitation
! _________________________________________________________________

! Read in LSP values
    status_ERA = nf90_get_var(ncid_FC,LSPVarId,scr,start=(/1, 1,itime /))
    if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
    LSP(2:nlons+1,:) = scr
    LSP(nlonp,:) = scr(1,:)
    LSP(1,:) = scr(nlons,:)
    LSP = LSP*LSP_scale_factor + LSP_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = LSP
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    LSP = scrp
! Smooth zonally
    CALL smooth_polar(scrp, LSP, nlonp, nlats)


!___________________________________________________________________
!
! Total precipitation
!
!____________________________________________________________________

   TP = CP + LSP


! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*TP(ia,ja)           +   &
                      rx*rym*TP(ia+1,ja)          +   &
                      rx*ry*TP(ia+1,ja+1)         +   &
                      rxm*ry*TP(ia,ja+1)
   ENDDO
   ENDDO

   rain = scr1_out*1000./delt

! Blank out Antarctica
  WHERE(ABS(rain)>1.e+7)
     rain = 0.
  ENDWHERE

! Write out total precipitation rate in netCDF file

   status_roms = nf90_put_var(ncid_rain,rainVarId,rain,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_rain,TimerainVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! _________________________________________________________________
!
! Downward shortwave radiation at the surface
! _________________________________________________________________

! Read in SSRD values
    status_ERA = nf90_get_var(ncid_FC,SSRDVarId,scr,start=(/1, 1,itime /))
    if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
    SSRD(2:nlons+1,:) = scr
    SSRD(nlonp,:) = scr(1,:)
    SSRD(1,:) = scr(nlons,:)
    SSRD = SSRD*SSRD_scale_factor + SSRD_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = SSRD
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    SSRD = scrp
! Smooth zonally
    CALL smooth_polar(scrp, SSRD, nlonp, nlats)

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*SSRD(ia,ja)           +   &
                      rx*rym*SSRD(ia+1,ja)          +   &
                      rx*ry*SSRD(ia+1,ja+1)         +   &
                      rxm*ry*SSRD(ia,ja+1)
   ENDDO
   ENDDO

   swrad = scr1_out/delt

! Blank out Antarctica
  WHERE(ABS(swrad)>1.e+7)
     swrad = 0.
  ENDWHERE

! Write out downward shortwave radiation in netCDF file

   status_roms = nf90_put_var(ncid_swrad,swradVarId,swrad,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF file

   status_roms = nf90_put_var(ncid_swrad,TimeswradVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)

! _________________________________________________________________
!
! Downward longwave radiation at the surface
! _________________________________________________________________

! Read in STRD values
    status_ERA = nf90_get_var(ncid_FC,STRDVarId,scr,start=(/1, 1,itime /))
    if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
    STRD(2:nlons+1,:) = scr
    STRD(nlonp,:) = scr(1,:)
    STRD(1,:) = scr(nlons,:)
    STRD = STRD*STRD_scale_factor + STRD_add_offset

! Fill in masked-out values 
    scrp = 0.
    scrp = STRD
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlonp,nlats,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work, &
              error,nvalue)
    STRD = scrp
! Smooth zonally
    CALL smooth_polar(scrp, STRD, nlonp, nlats)

! Horizontal interpolation
   scr1_out = 0.
   DO j=1,Mp
   DO i=1,Lp
      ia = ifix(ipos(i,j))
      ja = ifix(jpos(i,j))
      rx = (ipos(i,j)-ia)
      ry = (jpos(i,j)-ja)
      rxm = 1.0-rx
      rym = 1.0-ry
      scr1_out(i,j) = rxm*rym*STRD(ia,ja)           +   &
                      rx*rym*STRD(ia+1,ja)          +   &
                      rx*ry*STRD(ia+1,ja+1)         +   &
                      rxm*ry*STRD(ia,ja+1)
   ENDDO
   ENDDO

   lwrad = scr1_out/delt

! Blank out Antarctica
  WHERE(ABS(lwrad)>1.e+7)
     lwrad = 0.
  ENDWHERE

! Write out downward longwave radiation in netCDF file

   status_roms = nf90_put_var(ncid_lwrad,lwradVarId,lwrad,      &
                              start=(/1, 1, itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Write out time variable in netCDF files

   status_roms = nf90_put_var(ncid_lwrad,TimelwradVarId,tday,   &
                              start=(/ itime /))
   if(status_roms /= nf90_NoErr) call handle_err(status_roms)




!__________________________________________________________________

ENDDO

! End time frame loop
! ..................................
  status_era = nf90_close(ncid_FC)
   if(status_era /= nf90_NoErr) call handle_err(status_era)

! End 12-hourly (FC) data loop
! --------------------------------------------------------------


status_roms = nf90_close(ncid_rain)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_lwrad)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_swrad)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)


write(*,'(a)') 'Program get_forcing completed'

END PROGRAM get_forcing
