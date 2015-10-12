PROGRAM roms2agrid
! Program to interpolate all fields from C-grid to A-grid for ROMS output files

USE netcdf

IMPLICIT NONE

INTEGER :: statusi, ncgridi, ncidi
INTEGER :: dim_xi_rhoi, dim_eta_rhoi, dim_xi_ui, dim_eta_ui, dim_xi_vi, dim_eta_vi
INTEGER :: dim_s_rhoi, dim_timei
INTEGER :: ZetaVarIdi, UVarIdi, VVarIdi, SaltVarIdi, TempVarIdi, TimeVarIdi
INTEGER :: LonVarido, LatVarido, DepthVarido, AngleVarido, MaskVarido
INTEGER :: id_s_rhoi, id_theta_si, id_theta_bi, id_tclinei
INTEGER :: id_theta_so, id_theta_bo, id_tclineo
INTEGER :: id_lon_rhoi, id_lat_rhoi, id_anglei, id_hi
INTEGER :: id_rmaski, id_umaski, id_vmaski
INTEGER :: Lpi, Mpi, Li, Mi, Ni, Iti, Kt, Itc
REAL, DIMENSION(:), ALLOCATABLE   :: s_rhoi
REAL, DIMENSION(:,:), ALLOCATABLE :: lon_rhoi
REAL, DIMENSION(:,:), ALLOCATABLE :: lat_rhoi
REAL, DIMENSION(:,:), ALLOCATABLE :: anglei
REAL, DIMENSION(:,:), ALLOCATABLE :: hi
REAL, DIMENSION(:,:), ALLOCATABLE :: rmaski
REAL, DIMENSION(:,:), ALLOCATABLE :: umaski
REAL, DIMENSION(:,:), ALLOCATABLE :: vmaski
REAL :: hmini, Tclinei, theta_si, theta_bi, hci
REAL :: time_in, tday, t_yd, lmon, dt, f1, omf1
CHARACTER(len=100) :: xi_dimnamei, eta_dimnamei, s_dimnamei, time_dimnamei
CHARACTER(len=100) :: time_dimnamet, time_dimnamec, time_dimnamee, time_dimnames
INTEGER :: zetaout, velout, velout_polar, sout, tout

INTEGER :: statuso, ncido
INTEGER :: dim_xi_rhoo, dim_eta_rhoo
INTEGER :: dim_xi_uo, dim_eta_uo
INTEGER :: dim_xi_vo, dim_eta_vo
INTEGER :: dim_s_rhoo, dim_timeo
INTEGER :: UVarido, VVarido, SpeedVarido, DirecVarido, SaltVarido, TempVarido, TimeVarido
INTEGER :: ZetaVarido
INTEGER :: HVarido
INTEGER :: id_lon_rhoo, id_lat_rhoo, id_ho
INTEGER :: Lpo, Mpo, Lo, Mo, No
REAL, DIMENSION(:), ALLOCATABLE   :: s_rhoo
REAL, DIMENSION(:,:), ALLOCATABLE :: lon_rhoo
REAL, DIMENSION(:,:), ALLOCATABLE :: lat_rhoo
REAL, DIMENSION(:,:), ALLOCATABLE :: ho, rmasko
REAL :: Tclineo, theta_so, theta_bo
CHARACTER(len=80) :: xi_dimnameo, eta_dimnameo

INTEGER :: i, j, k, itime, jd, yy, mm, dd, hh
INTEGER :: iout
CHARACTER (len=100) :: ifile, ofile

REAL, DIMENSION(:,:), ALLOCATABLE :: work
REAL, DIMENSION(:,:), ALLOCATABLE :: scru
REAL, DIMENSION(:,:), ALLOCATABLE :: scrv
REAL, DIMENSION(:,:), ALLOCATABLE :: error

INTEGER :: nvalue
REAL, PARAMETER :: undef = 2.E+35            ! Undefined land value
REAL, PARAMETER :: critx = 0.01
REAL, PARAMETER :: cor = 1.6
REAL :: tx
INTEGER, PARAMETER :: mxs = 100

REAL :: rimin, rimax, rjmin, rjmax, rx, ry, rxm, rym, rz1, rz2
INTEGER :: i1, i2, j1, j2
REAL :: pi, deg, cosa, sina, u_east, v_north

REAL, DIMENSION(:,:,:), ALLOCATABLE :: temp_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: salt_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: u_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: v_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: u_rhoi
REAL, DIMENSION(:,:,:), ALLOCATABLE :: v_rhoi
REAL, DIMENSION(:,:),   ALLOCATABLE :: zeta_in

REAL, DIMENSION(:,:,:), ALLOCATABLE :: temp_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: salt_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: u_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: v_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: u_rhoo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: v_rhoo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: speed_rhoo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: direc_rhoo
REAL, DIMENSION(:,:),   ALLOCATABLE :: zeta_out

! -----------------------------------------------------------
!// Read standard input
READ(5,'(A)') ifile
READ(5,*) zetaout, velout, velout_polar, sout, tout
READ(5,'(A)') ofile

PRINT '(2A)', 'Input file with fields on C-grid= ', TRIM(ifile)
PRINT '(2A)', 'Output file with fields on A-grid= ', TRIM(ofile)

!// Set constants related to fill-routine
pi = ACOS(-1.0)
deg = 180./pi
tx = 0.9*undef

! ------------------------------------------------------------
! Get grid parameters from the first input file

statusi = nf90_open(TRIM(ifile),nf90_nowrite,ncgridi)
statusi = nf90_inq_dimid(ncgridi,'xi_rho',dim_xi_rhoi)
statusi = nf90_inq_dimid(ncgridi,'eta_rho',dim_eta_rhoi)
statusi = nf90_inq_dimid(ncgridi,'s_rho',dim_s_rhoi)
statusi = nf90_Inquire_Dimension(ncgridi,dim_xi_rhoi,xi_dimnamei,Lpi)
statusi = nf90_Inquire_Dimension(ncgridi,dim_eta_rhoi,eta_dimnamei,Mpi)
statusi = nf90_Inquire_Dimension(ncgridi,dim_s_rhoi,s_dimnamei,Ni)

statusi = nf90_inq_varid(ncgridi,'theta_s',id_theta_si); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'theta_b',id_theta_bi); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'Tcline',id_Tclinei);   IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_theta_si,theta_si);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_theta_bi,theta_bi);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_Tclinei,Tclinei);      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

Mi = Mpi-1
Li = Lpi-1

ALLOCATE(s_rhoi(Ni))
ALLOCATE(lon_rhoi(Lpi,Mpi))
ALLOCATE(lat_rhoi(Lpi,Mpi))
ALLOCATE(anglei(Lpi,Mpi))
ALLOCATE(hi(Lpi,Mpi))
ALLOCATE(rmaski(Lpi,Mpi))
ALLOCATE(umaski(Li,Mpi))
ALLOCATE(vmaski(Lpi,Mi))
ALLOCATE(work(Lpi,Mpi))
ALLOCATE(scru(Li,Mpi))
ALLOCATE(scrv(Lpi,Mi))
ALLOCATE(error(Lpi,Mpi))

statusi = nf90_inq_varid(ncgridi,'s_rho',id_s_rhoi);     IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'lon_rho',id_lon_rhoi); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'lat_rho',id_lat_rhoi); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'angle',id_anglei);     IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'h',id_hi);             IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'mask_rho',id_rmaski);  IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'mask_u',id_umaski);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'mask_v',id_vmaski);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_s_rhoi,s_rhoi);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_lon_rhoi,lon_rhoi);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_lat_rhoi,lat_rhoi);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_anglei,anglei);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_hi,hi);                IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_rmaski,rmaski);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_umaski,umaski);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_vmaski,vmaski);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

hmini = 1.0e+23
DO j=1,Mpi
DO i=1,Lpi
   hmini = MIN(hmini,hi(i,j))
ENDDO
ENDDO
PRINT '(A,I4,50F7.3)','Ni and s_rho (input-grid): ', Ni, s_rhoi
PRINT '(A,2F10.2)','Hmin and Tcli (input-grid): ', hmini, Tclinei
PRINT '(A,2F10.2)','Theta_s and Theta_b (input-grid): ', theta_si, theta_bi

statusi = nf90_close(ncgridi)

ALLOCATE(temp_in(Lpi,Mpi,Ni))
ALLOCATE(salt_in(Lpi,Mpi,Ni))
ALLOCATE(u_in(Li,Mpi,Ni))
ALLOCATE(v_in(Lpi,Mi,Ni))
ALLOCATE(zeta_in(Lpi,Mpi))

PRINT '(A,3I6)','Dimensions for input grid (Lpi, Mpi, Ni): ', Lpi, Mpi, Ni

! ---------------------------------------------------------------------

! Define info on output file

Lpo = Lpi
Mpo = Mpi
No = Ni
Mo = Mpo-1
Lo = Lpo-1
PRINT '(A,3I6)','Dimensions for output grid (Lpo, Mpo, No): ', Lpo, Mpo, No

ALLOCATE(s_rhoo(No))
ALLOCATE(lon_rhoo(Lpo,Mpo))
ALLOCATE(lat_rhoo(Lpo,Mpo))
ALLOCATE(ho(Lpo,Mpo))
ALLOCATE(rmasko(Lpo,Mpo))

ALLOCATE(temp_out(Lpo,Mpo,No))
ALLOCATE(salt_out(Lpo,Mpo,No))
ALLOCATE(u_out(Lpo,Mpo,No))
ALLOCATE(v_out(Lpo,Mpo,No))
ALLOCATE(u_rhoi(Lpi,Mpi,Ni))
ALLOCATE(v_rhoi(Lpi,Mpi,Ni))
ALLOCATE(u_rhoo(Lpo,Mpo,No))
ALLOCATE(v_rhoo(Lpo,Mpo,No))
ALLOCATE(speed_rhoo(Lpo,Mpo,No))
ALLOCATE(direc_rhoo(Lpo,Mpo,No))
ALLOCATE(zeta_out(Lpo,Mpo))

lon_rhoo = lon_rhoi
lat_rhoo = lat_rhoi

s_rhoo = s_rhoi
ho = hi
rmasko = rmaski

! .......................................................................
! Prepare output netCDF file
statuso = nf90_create(TRIM(ofile),nf90_clobber,ncido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Define dimensions
statuso = nf90_def_dim(ncido,'x',Lpo,dim_xi_rhoo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_dim(ncido,'y',Mpo,dim_eta_rhoo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_dim(ncido,'s_rho',No,dim_s_rhoo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_dim(ncido,'time',nf90_unlimited,dim_timeo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Define variables
statuso = nf90_def_var(ncido,'theta_s',nf90_double,varid=id_theta_so)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'theta_b',nf90_double,varid=id_theta_bo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'Tcline',nf90_double,varid=id_tclineo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_def_var(ncido,'lon',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo/),LonVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'lat',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo/),LatVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'s_rho',nf90_float,(/dim_s_rhoo/),DepthVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'h',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo/),HVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'mask',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo/),MaskVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'angle',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo/),AngleVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

IF (velout == 1) THEN
  statuso = nf90_def_var(ncido,'u',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),UVarido) 
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_def_var(ncido,'v',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),VVarido)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF

IF (velout_polar == 1) THEN
  statuso = nf90_def_var(ncido,'water_spd',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),SpeedVarido) 
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_def_var(ncido,'water_dir',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),DirecVarido) 
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (sout == 1) THEN
  statuso = nf90_def_var(ncido,'salt',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),SaltVarido) 
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (tout == 1) THEN
  statuso = nf90_def_var(ncido,'temp',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),TempVarido) 
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (zetaout == 1) THEN
  statuso = nf90_def_var(ncido,'zeta',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_timeo/),ZetaVarido) 
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF

statuso = nf90_def_var(ncido,'time',nf90_double,dim_timeo,TimeVarido) 
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Include variable attributes

statuso = nf90_put_att(ncido,DepthVarido,'long_name','S-coordinate at RHO-points')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'valid_min','-1.')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'valid_max','0.')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'positive','up')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'standard_name','ocean_s_coordinate_g1')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'field','s_rho, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,id_theta_so,'long_name','S-coordinate surface control parameter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,id_theta_bo,'long_name','S-coordinate bottom control parameter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,id_tclineo,'long_name','S-coordinate surface/bottom layer width')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,id_tclineo,'units','meter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,TimeVarido,'long_name','time since 1948/01/01/00:00')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,TimeVarido,'units','days')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,TimeVarido,'field','time, scalar, series')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,LonVarido,'long_name','longitude')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,LonVarido,'units','degree_east')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,LonVarido,'field','longitude, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,LatVarido,'long_name','latitude')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,LatVarido,'units','degree_north')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,LatVarido,'field','latitude, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,MaskVarido,'long_name','mask on RHO-points')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,MaskVarido,'flag_values','0., 1.')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,MaskVarido,'flag_meanings','land water')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,MaskVarido,'coordinates','lon lat')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,HVarido,'long_name','bathymetry')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,HVarido,'units','meter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,HVarido,'field','bathymetry, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,AngleVarido,'long_name','angle between XI-axis and EAST')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,AngleVarido,'units','radians')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,AngleVarido,'field','angle, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,AngleVarido,'coordinates','lon lat')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

IF (velout == 1) THEN
  statuso = nf90_put_att(ncido,UVarido,'long_name','east-momentum component')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'units','meter second-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'field','east-velocity, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

  statuso = nf90_put_att(ncido,VVarido,'long_name','north-momentum component')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VVarido,'units','meter second-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VVarido,'field','north-velocity, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (velout_polar == 1) THEN
  statuso = nf90_put_att(ncido,SpeedVarido,'long_name','current speed')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SpeedVarido,'units','meter second-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SpeedVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SpeedVarido,'field','current speed, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SpeedVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

  statuso = nf90_put_att(ncido,DirecVarido,'long_name','current direction')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'units','degrees')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'field','current direction, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'reference','clockwise from true north')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'valid_min',0.)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,DirecVarido,'valid_max',360.)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (sout == 1) THEN
  statuso = nf90_put_att(ncido,SaltVarido,'long_name','salinity')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SaltVarido,'units','PSU')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SaltVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SaltVarido,'field','salinity, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SaltVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (tout == 1) THEN
  statuso = nf90_put_att(ncido,TempVarido,'long_name','potential temperature')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,TempVarido,'units','Celsius')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,TempVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,TempVarido,'field','temperature, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,TempVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (zetaout == 1) THEN
  statuso = nf90_put_att(ncido,ZetaVarido,'long_name','free surface elevation')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,ZetaVarido,'units','meter')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,ZetaVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,ZetaVarido,'field','free surface elevation, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,ZetaVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF

statuso = nf90_put_att(ncido,nf90_global,'type','ROMS/TOMS processed result file')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,'title','ROMS/TOMS 3.3')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,'history','Made by program roms2z')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Exit definition mode
statuso = nf90_enddef(ncido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Output static fields to netCDF file
statuso = nf90_put_var(ncido,id_theta_so,theta_si); IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,id_theta_bo,theta_bi); IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,id_tclineo,Tclinei);   IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_var(ncido,LonVarido,lon_rhoo,(/1,1/));   IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,LatVarido,lat_rhoo,(/1,1/));   IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,DepthVarido,s_rhoo,(/1/));     IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,HVarido,ho,(/1,1/));           IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,MaskVarido,rmasko,(/1,1/));    IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,AngleVarido,anglei,(/1,1/));   IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Counter for output time steps
iout = 0

  ! Open input averages file
  statusi = nf90_open(TRIM(ifile),nf90_nowrite,ncidi);                 IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
  statusi = nf90_inq_dimid(ncidi,'ocean_time',dim_timei);              IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
  statusi = nf90_Inquire_Dimension(ncidi,dim_timei,time_dimnamei,Iti); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
  statusi = nf90_inq_varid(ncidi,'ocean_time',TimeVarIdi);             IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

  ! Loop through time steps in ROMS file
  PRINT '(A)','Start looping through time steps in inputfile'

  DO itime = 1,Iti

    ! Read in time
    statusi = nf90_inq_varid(ncidi,'ocean_time',TimeVarIdi)
    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
    statusi = nf90_get_var(ncidi,TimeVarIdi,time_in,start=(/itime/))
    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
    tday = time_in/86400.  ! seconds -> days

    CALL gregorian(tday+jd(1948,1,1),yy,mm,dd,hh)
    PRINT '(A,4I5)','Date on input file is (Y M D H) ', yy, mm, dd, hh

    ! Update counter of time steps in input file
    iout = iout + 1

    IF (zetaout == 1) THEN
 
      PRINT '(A)','Find sea level'

      ! Read in sea surface elevation
      statusi = nf90_inq_varid(ncidi,'zeta',ZetaVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,ZetaVarIdi,zeta_in,start=(/ 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      zeta_out = zeta_in
      WHERE (rmasko(:,:) < 0.01) zeta_out(:,:) = 1.E+37

      ! Output to netCDF file
      statuso = nf90_put_var(ncido,ZetaVarido,zeta_out,(/1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

    ENDIF

    IF (tout == 1) THEN

      PRINT '(A)','Find temperature'

      ! Read in temperature
      statusi = nf90_inq_varid(ncidi,'temp',TempVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,TempVarIdi,temp_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      temp_out = temp_in
      DO k=1,No
        WHERE (rmasko(:,:) < 0.01) temp_out(:,:,k) = 1.E+37
      ENDDO

      ! Output to netCDF file - z-levels
      statuso = nf90_put_var(ncido,TempVarido,temp_out,(/1,1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

    ENDIF

    IF (sout == 1) THEN

      PRINT '(A)','Find salinity'

      ! Read in salinity
      statusi = nf90_inq_varid(ncidi,'salt',SaltVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,SaltVarIdi,salt_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      salt_out = salt_in
      DO k=1,No
        WHERE (rmasko(:,:) < 0.01) salt_out(:,:,k) = 1.E+37
      ENDDO

      ! Output to netCDF file
      statuso = nf90_put_var(ncido,SaltVarido,salt_out,(/1,1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

    ENDIF

    IF (velout == 1 .OR. velout_polar == 1) THEN

      PRINT '(A)','Find 3D velocity'

      ! Read in velocities
      statusi = nf90_inq_varid(ncidi,'u',UVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,UVarIdi,u_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_inq_varid(ncidi,'v',VVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,VVarIdi,v_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      ! Fill in masked-out values before interpolation to A-grid
      i1 = 1
      i2 = Li
      j1 = 1
      j2 = Mpi

      DO k=1,Ni
        scru = 0.
        scru = u_in(:,:,k)
        WHERE (umaski < 0.01) scru = undef
        CALL fill(Li,Mpi,i1,i2,j1,j2,scru,tx,critx,cor,mxs,work,error,nvalue)
        u_in(:,:,k) = scru
      ENDDO

      i1 = 1
      i2 = Lpi
      j1 = 1
      j2 = Mi

      DO k=1,Ni
        scrv = 0.
        scrv = v_in(:,:,k)
        WHERE (vmaski < 0.01) scrv = undef
        CALL fill(Lpi,Mi,i1,i2,j1,j2,scrv,tx,critx,cor,mxs,work,error,nvalue)
        v_in(:,:,k) = scrv
      ENDDO

      ! Average u and v to rho locations
      u_rhoi(2:Li,:,:) = 0.5*(u_in(1:Li-1,:,:)+u_in(2:Li,:,:))
      v_rhoi(:,2:Mi,:) = 0.5*(v_in(:,1:Mi-1,:)+v_in(:,2:Mi,:))
      u_rhoi(1,:,:) = u_in(1,:,:)
      u_rhoi(Lpi,:,:) = u_in(Li,:,:)
      v_rhoi(:,1,:) = v_in(:,1,:)
      v_rhoi(:,Mpi,:) = v_in(:,Mi,:)

      u_rhoo = u_rhoi
      v_rhoo = v_rhoi

      DO k=1,No
        WHERE (rmasko(:,:) < 0.01) u_rhoo(:,:,k) = 1.E+37
        WHERE (rmasko(:,:) < 0.01) v_rhoo(:,:,k) = 1.E+37
      ENDDO

      IF (velout_polar == 1) THEN

        DO k=1,No
        DO j=1,Mpo
        DO i=1,Lpo
          IF (rmasko(i,j) > 0.5) THEN
            cosa = COS(anglei(i,j))
            sina = SIN(anglei(i,j))
            u_east  = u_rhoo(i,j,k) * cosa - v_rhoo(i,j,k) * sina
            v_north = u_rhoo(i,j,k) * sina + v_rhoo(i,j,k) * cosa
            speed_rhoo(i,j,k) = SQRT(u_rhoo(i,j,k)*u_rhoo(i,j,k)+v_rhoo(i,j,k)*v_rhoo(i,j,k))
            IF (ABS(u_east) > 1.E-8) direc_rhoo(i,j,k) = ATAN2(v_north,u_east)*deg  ! ATAN2 elements in [-180,180]
            IF (u_east > 1.E-6) THEN
              direc_rhoo(i,j,k) = 90. - direc_rhoo(i,j,k)
            ELSEIF (u_east < -1.E-8) THEN
              IF (v_north >= 0.) THEN
                direc_rhoo(i,j,k) = 450. - direc_rhoo(i,j,k)
              ELSE
                direc_rhoo(i,j,k) = 90. - direc_rhoo(i,j,k)
              ENDIF
            ELSE
              IF (v_north >= 0.) THEN
                direc_rhoo(i,j,k) = 0.
              ELSE
                direc_rhoo(i,j,k) = 180.
              ENDIF
            ENDIF
          ELSE
            speed_rhoo(i,j,k) = 1.E+37
            direc_rhoo(i,j,k) = 1.E+37
          ENDIF
        ENDDO
        ENDDO
        ENDDO

      ENDIF

      ! Output to netCDF file
      IF (velout == 1) THEN
        statuso = nf90_put_var(ncido,UVarido,u_rhoo,(/1,1,1,iout/))
        IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
        statuso = nf90_put_var(ncido,VVarido,v_rhoo,(/1,1,1,iout/))
        IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
      ENDIF
      IF (velout_polar == 1) THEN
        statuso = nf90_put_var(ncido,SpeedVarido,speed_rhoo,(/1,1,1,iout/))
        IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
        statuso = nf90_put_var(ncido,DirecVarido,direc_rhoo,(/1,1,1,iout/))
        IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
      ENDIF

    ENDIF  ! (velout == 1 .OR. velout_polar == 1)

    ! Write out time
    statuso = nf90_put_var(ncido,TimeVarido,tday,(/iout/))
    IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

  ENDDO  ! timesteps

statuso = nf90_sync(ncido)
statuso = nf90_close(ncido)

END PROGRAM roms2agrid

