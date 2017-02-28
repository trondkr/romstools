PROGRAM roms2z
! Program to obtain lon-lat-lev_z from ROMS output files

USE netcdf

IMPLICIT NONE

INTEGER :: statusi, ncgridi, ncidi
INTEGER :: dim_xi_rhoi, dim_eta_rhoi, dim_xi_ui, dim_eta_ui, dim_xi_vi, dim_eta_vi
INTEGER :: dim_s_rhoi, dim_timei
INTEGER :: ZetaVarIdi, SHfluxVarIdi, SSfluxVarIdi
INTEGER :: UVarIdi, VVarIdi, SaltVarIdi, TempVarIdi, UEVarIdi, VNVarIdi, TimeVarIdi
INTEGER :: Mask3dVarId, LonVarido, LatVarido, DepthVarido, MaskRhoVarido, AngleVarido
INTEGER :: id_theta_si, id_theta_bi, id_tclinei
INTEGER :: id_theta_so, id_theta_bo, id_tclineo
INTEGER :: id_lon_rhoi, id_lat_rhoi, id_anglei, id_hi
INTEGER :: id_rmaski, id_umaski, id_vmaski
INTEGER :: Lpi, Mpi, Li, Mi, Ni, Iti, Kt, Itc
REAL, DIMENSION(:,:), ALLOCATABLE :: lon_rhoi
REAL, DIMENSION(:,:), ALLOCATABLE :: lat_rhoi
REAL, DIMENSION(:,:), ALLOCATABLE :: anglei
REAL, DIMENSION(:,:), ALLOCATABLE :: hi
REAL, DIMENSION(:,:), ALLOCATABLE :: rmaski
REAL, DIMENSION(:,:), ALLOCATABLE :: umaski
REAL, DIMENSION(:,:), ALLOCATABLE :: vmaski
REAL, DIMENSION(:), ALLOCATABLE   :: sc_ri
REAL, DIMENSION(:), ALLOCATABLE   :: Cs_ri
REAL, DIMENSION(:,:,:), ALLOCATABLE :: z_ri
REAL, DIMENSION(:,:,:), ALLOCATABLE :: mask3d
REAL :: hmini, Tclinei, theta_si, theta_bi, hci, cff1i, cff2i
REAL :: cffi, cff_ri, cff1_ri, cff2_ri
REAL :: time_in, tday, tday_prev, t_yd, lmon, dt, f1, omf1
CHARACTER(len=100) :: xi_dimnamei, eta_dimnamei, s_dimnamei, time_dimnamei
CHARACTER(len=100) :: time_dimnamet, time_dimnamec, time_dimnamee, time_dimnames
INTEGER :: zetaout, fluxout, velout, velout_rot, velout_polar, sout, tout

INTEGER :: statuso, ncido
INTEGER :: dim_xi_rhoo, dim_eta_rhoo
INTEGER :: dim_xi_uo, dim_eta_uo
INTEGER :: dim_xi_vo, dim_eta_vo
INTEGER :: dim_s_rhoo, dim_timeo
INTEGER :: UVarido, VVarido, SpeedVarido, DirecVarido, UEVarido, VNVarido, SaltVarido, TempVarido, TimeVarido
INTEGER :: ZetaVarido, SHfluxVarido, SSfluxVarido
INTEGER :: HVarido
INTEGER :: id_lon_rhoo, id_lat_rhoo, id_ho
INTEGER :: Lpo, Mpo, Lo, Mo, No
REAL, DIMENSION(:,:), ALLOCATABLE :: lon_rhoo
REAL, DIMENSION(:,:), ALLOCATABLE :: lat_rhoo
REAL, DIMENSION(:,:), ALLOCATABLE :: ho, rmasko
REAL, DIMENSION(:), ALLOCATABLE   :: sc_ro
REAL, DIMENSION(:), ALLOCATABLE   :: Cs_ro
REAL, DIMENSION(:), ALLOCATABLE :: z_ro
REAL :: hmino, Tclineo, theta_so, theta_bo, hco, cff1o, cff2o
REAL :: cffo, cff_ro, cff1_ro, cff2_ro
CHARACTER(len=80) :: xi_dimnameo, eta_dimnameo

INTEGER :: i, j, k, itime, jd, yy, mm, dd, hh
INTEGER :: iout
CHARACTER (len=300) :: listfile, avgfilei, ofile

REAL, DIMENSION(:,:), ALLOCATABLE :: ipos
REAL, DIMENSION(:,:), ALLOCATABLE :: jpos
REAL, DIMENSION(:,:), ALLOCATABLE :: work
REAL, DIMENSION(:,:), ALLOCATABLE :: scr
REAL, DIMENSION(:,:), ALLOCATABLE :: scr1
REAL, DIMENSION(:,:), ALLOCATABLE :: scru
REAL, DIMENSION(:,:), ALLOCATABLE :: scrv
REAL, DIMENSION(:,:), ALLOCATABLE :: scr1o
REAL, DIMENSION(:,:), ALLOCATABLE :: scr2o
REAL, DIMENSION(:,:), ALLOCATABLE :: scruo
REAL, DIMENSION(:,:), ALLOCATABLE :: scrvo
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
REAL, DIMENSION(:,:),   ALLOCATABLE :: shflux_in
REAL, DIMENSION(:,:),   ALLOCATABLE :: ssflux_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ue_in
REAL, DIMENSION(:,:,:), ALLOCATABLE :: vn_in

REAL, DIMENSION(:,:,:), ALLOCATABLE :: temp_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: salt_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: u_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: v_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: u_rhoo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: v_rhoo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: speed_rhoo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: direc_rhoo
REAL, DIMENSION(:,:),   ALLOCATABLE :: zeta_out
REAL, DIMENSION(:,:),   ALLOCATABLE :: shflux_out
REAL, DIMENSION(:,:),   ALLOCATABLE :: ssflux_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ue_out
REAL, DIMENSION(:,:,:), ALLOCATABLE :: vn_out

REAL, DIMENSION(:,:,:), ALLOCATABLE :: scr3d

REAL, DIMENSION(:), ALLOCATABLE :: zlev
INTEGER :: nlev

INTEGER :: nl, nw

! -----------------------------------------------------------
!// Read standard input
READ(5,'(A)') listfile
READ(5,*) zetaout, fluxout, velout, velout_rot, velout_polar, sout, tout
READ(5,'(A)') ofile
READ(5,*) nlev
ALLOCATE(zlev(nlev))
READ(5,*) zlev

PRINT '(2A)', 'Output file with fields on z-levels= ', TRIM(ofile)

!// Set constants related to fill-routine
pi = ACOS(-1.0)
deg = 180./pi
tx = 0.9*undef

! Read first file in listfile
OPEN(UNIT=55,FILE=listfile,FORM='FORMATTED',ACTION='READ')
READ(55,'(A)',END=999) avgfilei
PRINT '(2A)','First file to process is ', TRIM(avgfilei)
REWIND(55)  ! Start all over within loop over input files

! ------------------------------------------------------------
! Get grid parameters from the first input file

statusi = nf90_open(TRIM(avgfilei),nf90_nowrite,ncgridi)
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

ALLOCATE(lon_rhoi(Lpi,Mpi))
ALLOCATE(lat_rhoi(Lpi,Mpi))
ALLOCATE(anglei(Lpi,Mpi))
ALLOCATE(hi(Lpi,Mpi))
ALLOCATE(rmaski(Lpi,Mpi))
ALLOCATE(umaski(Li,Mpi))
ALLOCATE(vmaski(Lpi,Mi))
ALLOCATE(work(Lpi,Mpi))
ALLOCATE(scr(Lpi,Mpi))
ALLOCATE(scr1(Lpi,Mpi))
ALLOCATE(scru(Li,Mpi))
ALLOCATE(scrv(Lpi,Mi))
ALLOCATE(error(Lpi,Mpi))

statusi = nf90_inq_varid(ncgridi,'lon_rho',id_lon_rhoi); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'lat_rho',id_lat_rhoi); IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'angle',id_anglei);     IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'h',id_hi);             IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'mask_rho',id_rmaski);  IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'mask_u',id_umaski);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_inq_varid(ncgridi,'mask_v',id_vmaski);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_lon_rhoi,lon_rhoi);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_lat_rhoi,lat_rhoi);    IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_anglei,anglei);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_hi,hi);                IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_rmaski,rmaski);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_umaski,umaski);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
statusi = nf90_get_var(ncgridi,id_vmaski,vmaski);        IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

hmini = 1000

hmini=MINVAL(hi)
print*,"test",hmini,Tclinei,theta_si,theta_bi
PRINT '(A,2F10.2)','Hmin and Tcli (input-grid): ', hmini, Tclinei
PRINT '(A,2F10.2)','Theta_s and Theta_b (input-grid): ', theta_si, theta_bi

! Vertical grid specification for input

ALLOCATE(z_ri(Lpi,Mpi,Ni))
ALLOCATE(sc_ri(Ni))
ALLOCATE(Cs_ri(Ni))

hci = MIN(hmini,Tclinei)

IF (theta_si /= 0.0) THEN

  cff1i=1.0/SINH(theta_si)
  cff2i=0.5/TANH(0.5*theta_si)
END IF
cffi=1.0/REAL(Ni)
DO k=1,Ni
  sc_ri(k)=cffi*(REAL(k-Ni)-0.5)
  IF (theta_si /= 0.0) THEN
    Cs_ri(k) = (1.0-theta_bi)*cff1i*SINH(theta_si*sc_ri(k)) +         &
                    theta_bi*(cff2i*TANH(theta_si*(sc_ri(k)+0.5))- 0.5)
  ELSE
    Cs_ri(k)=sc_ri(k)
  END IF
END DO
DO j=1,Mpi
  DO k=1,Ni
    cff_ri=hci*(sc_ri(k)-Cs_ri(k))
    cff1_ri=Cs_ri(k)
    cff2_ri=sc_ri(k)+1.0
    DO i=1,Lpi
      z_ri(i,j,k)=cff_ri+cff1_ri*hi(i,j)
    END DO
  END DO
END DO

PRINT '(A,50F9.2)','z-levels from input grid, z_ri(1,1,:)', z_ri(1,1,:)

statusi = nf90_close(ncgridi)

ALLOCATE(temp_in(Lpi,Mpi,Ni))
ALLOCATE(salt_in(Lpi,Mpi,Ni))
ALLOCATE(u_in(Li,Mpi,Ni))
ALLOCATE(v_in(Lpi,Mi,Ni))
ALLOCATE(zeta_in(Lpi,Mpi))
ALLOCATE(shflux_in(Lpi,Mpi))
ALLOCATE(ssflux_in(Lpi,Mpi))
ALLOCATE(ue_in(Lpi,Mpi,Ni))
ALLOCATE(vn_in(Lpi,Mpi,Ni))

PRINT '(A,3I6)','Dimensions for input grid (Lpi, Mpi, Ni): ', Lpi, Mpi, Ni

! ---------------------------------------------------------------------

! Define info on output file

Lpo = Lpi
Mpo = Mpi
No = nlev
Mo = Mpo-1
Lo = Lpo-1
PRINT '(A,3I6)','Dimensions for output grid (Lpo, Mpo, No): ', Lpo, Mpo, No

ALLOCATE(lon_rhoo(Lpo,Mpo))
ALLOCATE(lat_rhoo(Lpo,Mpo))
ALLOCATE(ho(Lpo,Mpo))
ALLOCATE(rmasko(Lpo,Mpo))
ALLOCATE(sc_ro(No))
ALLOCATE(Cs_ro(No))
ALLOCATE(z_ro(No))
ALLOCATE(scr3d(Lpi,Mpi,No))
ALLOCATE(scr1o(Lpo,Mpo))
ALLOCATE(scr2o(Lpo,Mpo))
ALLOCATE(scruo(Lo,Mpo))
ALLOCATE(scrvo(Lpo,Mo))

ALLOCATE(mask3d(Lpo,Mpo,No))
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
ALLOCATE(shflux_out(Lpo,Mpo))
ALLOCATE(ssflux_out(Lpo,Mpo))
ALLOCATE(ue_out(Lpo,Mpo,No))
ALLOCATE(vn_out(Lpo,Mpo,No))

ALLOCATE(ipos(Lpo,Mpo))
ALLOCATE(jpos(Lpo,Mpo))

lon_rhoo = lon_rhoi
lat_rhoo = lat_rhoi

! Vertical grid specification for output

DO k=1,No
  z_ro(k) = -zlev(k)
END DO
PRINT '(A,50F9.2)','z_ro: ',z_ro

ho = hi
rmasko = rmaski

! Create 3D mask
DO k=1,No
  scr1o(:,:) = 0.
  WHERE (rmasko(:,:) > 0.5)  scr1o(:,:) = 1.
  WHERE (-ho(:,:) > z_ro(k)) scr1o(:,:) = 0.
  mask3d(:,:,k) = scr1o(:,:)
ENDDO

nl = 0; nw = 0
DO j=1,Mpo
DO i=1,Lpo
  IF (rmasko(i,j) < 0.5) THEN
    nl = nl + 1
  ELSE
    nw = nw + 1
  ENDIF
ENDDO
ENDDO
PRINT '(A,4I8)','Mpo*Lpo, n-land and n-wet and sum (from rmasko): ', Mpo*Lpo, nl, nw, nl+nw

DO k=1,no
  nl = 0; nw = 0
  DO j=1,Mpo
  DO i=1,Lpo
    IF (mask3d(i,j,k) < 0.5) THEN
      nl = nl + 1
    ELSE
      nw = nw + 1
    ENDIF
  ENDDO
  ENDDO
  PRINT '(A,5I8)','k, Mpo*Lpo, n-land and n-wet and sum from mask3d: ', k, Mpo*Lpo, nl, nw, nl+nw
ENDDO

! .......................................................................
! Prepare output netCDF file
statuso = nf90_create(TRIM(ofile),nf90_clobber,ncido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Define dimensions
statuso = nf90_def_dim(ncido,'x',Lpo,dim_xi_rhoo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_dim(ncido,'y',Mpo,dim_eta_rhoo)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_dim(ncido,'depth',No,dim_s_rhoo)
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
statuso = nf90_def_var(ncido,'depth',nf90_float,(/dim_s_rhoo/),DepthVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'h',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo/),HVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_def_var(ncido,'mask3d',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo/),Mask3dVarId)
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
IF (fluxout == 1) THEN
  statuso = nf90_def_var(ncido,'shflux',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_timeo/),SHfluxVarido)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_def_var(ncido,'ssflux',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_timeo/),SSfluxVarido)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (velout_rot == 1) THEN
  statuso = nf90_def_var(ncido,'u_eastward',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),UEVarido)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_def_var(ncido,'v_northward',nf90_float,(/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),VNVarido)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF

statuso = nf90_def_var(ncido,'time',nf90_double,dim_timeo,TimeVarido)
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Include variable attributes

statuso = nf90_put_att(ncido,id_theta_so,'long_name','S-coordinate surface control parameter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,id_theta_bo,'long_name','S-coordinate bottom control parameter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,id_tclineo,'long_name','S-coordinate surface/bottom layer width')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,id_tclineo,'units','meter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,TimeVarido,'long_name','averaged time since initialization')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,TimeVarido,'units','days since 1948-01-01 00:00:00')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,TimeVarido,'calendar','gregorian')
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

statuso = nf90_put_att(ncido,DepthVarido,'long_name','depth of level')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'units','m')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'positive','down')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,DepthVarido,'field','depth, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,HVarido,'long_name','bathymetry')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,HVarido,'units','meter')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,HVarido,'coordinates','lon lat')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,HVarido,'field','bathymetry, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,Mask3dVarId,'long_name','3D mask')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,Mask3dVarId,'units','none')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,Mask3dVarId,'flag_values','0., 1.')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,Mask3dVarId,'flag_meanings','land water')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,Mask3dVarId,'coordinates','lon lat')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,Mask3dVarId,'field','mask, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

statuso = nf90_put_att(ncido,AngleVarido,'long_name','angle between XI-axis and EAST')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,AngleVarido,'units','radians')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,AngleVarido,'coordinates','lon lat')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_att(ncido,AngleVarido,'field','angle, scalar')
IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

IF (velout == 1) THEN
  statuso = nf90_put_att(ncido,UVarido,'long_name','east-momentum component')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'units','meter second-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UVarido,'coordinates','lon lat')
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
  statuso = nf90_put_att(ncido,VVarido,'coordinates','lon lat')
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
  statuso = nf90_put_att(ncido,SpeedVarido,'coordinates','lon lat')
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
  statuso = nf90_put_att(ncido,DirecVarido,'coordinates','lon lat')
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
  statuso = nf90_put_att(ncido,SaltVarido,'coordinates','lon lat')
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
  statuso = nf90_put_att(ncido,TempVarido,'coordinates','lon lat')
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
  statuso = nf90_put_att(ncido,ZetaVarido,'coordinates','lon lat')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,ZetaVarido,'field','free surface elevation, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,ZetaVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (fluxout == 1) THEN
  statuso = nf90_put_att(ncido,SHfluxVarido,'long_name','surface net heat flux')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'units','W m-2')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'negative_value','upward flux, cooling')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'positive_value','downward flux, heating')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'coordinates','lon lat')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'field','surface net heat flux, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SHfluxVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'long_name','surface net salt flux, (E-P)*SALT')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'units','m s-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'negative_value','upward flux, freshening (net precipitation)')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'positive_value','downward flux, salting (net evaporation)')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'coordinates','lon lat')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'field','surface net salt flux, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,SSfluxVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
ENDIF
IF (velout_rot == 1) THEN
  statuso = nf90_put_att(ncido,UEVarido,'long_name','eastward momentum component')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UEVarido,'units','m s-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UEVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UEVarido,'coordinates','lon lat')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UEVarido,'field','u_eastward, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,UEVarido,'_FillValue',1.E+37)
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VNVarido,'long_name','northward momentum component')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VNVarido,'units','m s-1')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VNVarido,'time','time')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VNVarido,'coordinates','lon lat')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VNVarido,'field','v_northward, scalar, series')
  IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
  statuso = nf90_put_att(ncido,VNVarido,'_FillValue',1.E+37)
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
statuso = nf90_put_var(ncido,DepthVarido,zlev,(/1/));       IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,HVarido,ho,(/1,1/));           IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,Mask3dVarId,mask3d,(/1,1,1/)); IF (statuso /= nf90_NoErr) CALL handle_err(statuso)
statuso = nf90_put_var(ncido,AngleVarido,anglei,(/1,1/));   IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

! Counter for output time steps
iout = 0

tday = 0.

DO WHILE (.TRUE.)

  ! Read filename in listfile
  READ(55,'(A)',END=999) avgfilei
  PRINT '(2A)', 'Read fields from ', TRIM(avgfilei)

  ! Open input averages file
  statusi = nf90_open(TRIM(avgfilei),nf90_nowrite,ncidi);              IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
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

    tday = time_in/86400.  ! seconds -> days for new time stamp

    CALL gregorian(tday+jd(1948,1,1),yy,mm,dd,hh)
    PRINT '(A,4I5)','Date on input file is (Y M D H) ', yy, mm, dd, hh

    IF (tday <= tday_prev) THEN
      PRINT *,'Repetition of time stamp: ', time_in, tday
    ELSE  ! Continue

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

    IF (fluxout == 1) THEN

      PRINT '(A)','Find surface fluxes'

      ! Read in surface net heat flux
      statusi = nf90_inq_varid(ncidi,'shflux',SHfluxVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,SHfluxVarIdi,shflux_in,start=(/ 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      shflux_out = shflux_in
      WHERE (rmasko(:,:) < 0.01) shflux_out(:,:) = 1.E+37

      ! Output to netCDF file
      statuso = nf90_put_var(ncido,SHfluxVarido,shflux_out,(/1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

      ! Read in surface net salt flux
      statusi = nf90_inq_varid(ncidi,'ssflux',SSfluxVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,SSfluxVarIdi,ssflux_in,start=(/ 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      ssflux_out = ssflux_in
      WHERE (rmasko(:,:) < 0.01) ssflux_out(:,:) = 1.E+37

      ! Output to netCDF file
      statuso = nf90_put_var(ncido,SSfluxVarido,ssflux_out,(/1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

    ENDIF

    IF (tout == 1) THEN

      PRINT '(A)','Find temperature'

      ! Read in temperature
      statusi = nf90_inq_varid(ncidi,'temp',TempVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,TempVarIdi,temp_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      ! Vertical interpolation
      scr3d = 0.
      DO j=1,Mpi
      DO i=1,Lpi
        DO k=1,No
          IF (z_ro(k).LT.z_ri(i,j,1)) THEN
            scr3d(i,j,k) = temp_in(i,j,1)
          ELSEIF (z_ro(k).GT.z_ri(i,j,Ni)) THEN
            scr3d(i,j,k) = temp_in(i,j,Ni)
          ELSE
            DO kT=1,Ni
              IF (z_ro(k).LT.z_ri(i,j,kT+1) .AND. z_ro(k).GE.z_ri(i,j,kT)) THEN
                rz2 = (z_ro(k)-z_ri(i,j,kT))/(z_ri(i,j,kT+1)-z_ri(i,j,kT))
                rz1 = 1.0 - rz2
                scr3d(i,j,k) = rz1*temp_in(i,j,kT) + rz2*temp_in(i,j,kT+1)
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      temp_out = scr3d
      DO k=1,No
        WHERE (mask3d(:,:,k) < 0.5) temp_out(:,:,k) = 1.E+37
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

      ! Vertical interpolation
      scr3d = 0.
      DO j=1,Mpi
      DO i=1,Lpi
        DO k=1,No
          IF (z_ro(k).LT.z_ri(i,j,1)) THEN
            scr3d(i,j,k) = salt_in(i,j,1)
          ELSEIF (z_ro(k).GT.z_ri(i,j,Ni)) THEN
            scr3d(i,j,k) = salt_in(i,j,Ni)
          ELSE
            DO kT=1,Ni
              IF (z_ro(k).LT.z_ri(i,j,kT+1) .AND. z_ro(k).GE.z_ri(i,j,kT)) THEN
                rz2 = (z_ro(k)-z_ri(i,j,kT))/(z_ri(i,j,kT+1)-z_ri(i,j,kT))
                rz1 = 1.0 - rz2
                scr3d(i,j,k) = rz1*salt_in(i,j,kT) + rz2*salt_in(i,j,kT+1)
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      salt_out = scr3d
      DO k=1,No
        WHERE (mask3d(:,:,k) < 0.01) salt_out(:,:,k) = 1.E+37
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

      ! Vertical interpolation of u_rho
      scr3d = 0.
      DO j=1,Mpi
      DO i=1,Lpi
        DO k=1,No
          IF (z_ro(k).LT.z_ri(i,j,1)) THEN
            scr3d(i,j,k) = u_rhoi(i,j,1)
          ELSEIF (z_ro(k).GT.z_ri(i,j,Ni)) THEN
            scr3d(i,j,k) = u_rhoi(i,j,Ni)
          ELSE
            DO kT=1,Ni
              IF (z_ro(k).LT.z_ri(i,j,kT+1) .AND. z_ro(k).GE.z_ri(i,j,kT)) THEN
                rz2 = (z_ro(k)-z_ri(i,j,kT))/(z_ri(i,j,kT+1)-z_ri(i,j,kT))
                rz1 = 1.0 - rz2
                scr3d(i,j,k) = rz1*u_rhoi(i,j,kT) + rz2*u_rhoi(i,j,kT+1)
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      u_rhoo = scr3d

      ! Vertical interpolation of v_rho
      scr3d = 0.
      DO j=1,Mpi
      DO i=1,Lpi
        DO k=1,No
          IF (z_ro(k).LT.z_ri(i,j,1)) THEN
            scr3d(i,j,k) = v_rhoi(i,j,1)
          ELSEIF (z_ro(k).GT.z_ri(i,j,Ni)) THEN
            scr3d(i,j,k) =v_rhoi(i,j,Ni)
          ELSE
            DO kT=1,Ni
              IF (z_ro(k).LT.z_ri(i,j,kT+1) .AND. z_ro(k).GE.z_ri(i,j,kT)) THEN
                rz2 = (z_ro(k)-z_ri(i,j,kT))/(z_ri(i,j,kT+1)-z_ri(i,j,kT))
                rz1 = 1.0 - rz2
                scr3d(i,j,k) = rz1*v_rhoi(i,j,kT) + rz2*v_rhoi(i,j,kT+1)
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      v_rhoo = scr3d

      DO k=1,No
        WHERE(mask3d(:,:,k) < 0.01) u_rhoo(:,:,k) = 1.E+37
        WHERE(mask3d(:,:,k) < 0.01) v_rhoo(:,:,k) = 1.E+37
      ENDDO

      IF (velout_polar == 1) THEN

        DO k=1,No
        DO j=1,Mpo
        DO i=1,Lpo
          IF (mask3d(i,j,k) > 0.5) THEN
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

    IF (velout_rot == 1) THEN

      PRINT '(A)','Find eastward and northward current component'

      ! Read in eastward component
      statusi = nf90_inq_varid(ncidi,'u_eastward',UEVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,UEVarIdi,ue_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      ! Vertical interpolation
      scr3d = 0.
      DO j=1,Mpi
      DO i=1,Lpi
        DO k=1,No
          IF (z_ro(k).LT.z_ri(i,j,1)) THEN
            scr3d(i,j,k) = ue_in(i,j,1)
          ELSEIF (z_ro(k).GT.z_ri(i,j,Ni)) THEN
            scr3d(i,j,k) = ue_in(i,j,Ni)
          ELSE
            DO kT=1,Ni
              IF (z_ro(k).LT.z_ri(i,j,kT+1) .AND. z_ro(k).GE.z_ri(i,j,kT)) THEN
                rz2 = (z_ro(k)-z_ri(i,j,kT))/(z_ri(i,j,kT+1)-z_ri(i,j,kT))
                rz1 = 1.0 - rz2
                scr3d(i,j,k) = rz1*ue_in(i,j,kT) + rz2*ue_in(i,j,kT+1)
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      ue_out = scr3d
      DO k=1,No
        WHERE (mask3d(:,:,k) < 0.5) ue_out(:,:,k) = 1.E+37
      ENDDO

      ! Output to netCDF file - z-levels
      statuso = nf90_put_var(ncido,UEVarido,ue_out,(/1,1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

      ! Read in northward component
      statusi = nf90_inq_varid(ncidi,'v_northward',VNVarIdi)
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)
      statusi = nf90_get_var(ncidi,VNVarIdi,vn_in,start=(/ 1, 1, 1, itime/))
      IF (statusi /= nf90_NoErr) CALL handle_err(statusi)

      ! Vertical interpolation
      scr3d = 0.
      DO j=1,Mpi
      DO i=1,Lpi
        DO k=1,No
          IF (z_ro(k).LT.z_ri(i,j,1)) THEN
            scr3d(i,j,k) = vn_in(i,j,1)
          ELSEIF (z_ro(k).GT.z_ri(i,j,Ni)) THEN
            scr3d(i,j,k) = vn_in(i,j,Ni)
          ELSE
            DO kT=1,Ni
              IF (z_ro(k).LT.z_ri(i,j,kT+1) .AND. z_ro(k).GE.z_ri(i,j,kT)) THEN
                rz2 = (z_ro(k)-z_ri(i,j,kT))/(z_ri(i,j,kT+1)-z_ri(i,j,kT))
                rz1 = 1.0 - rz2
                scr3d(i,j,k) = rz1*vn_in(i,j,kT) + rz2*vn_in(i,j,kT+1)
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ENDDO

      vn_out = scr3d
      DO k=1,No
        WHERE (mask3d(:,:,k) < 0.5) vn_out(:,:,k) = 1.E+37
      ENDDO

      ! Output to netCDF file - z-levels
      statuso = nf90_put_var(ncido,VNVarido,vn_out,(/1,1,1,iout/))
      IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

    ENDIF

    ! Write out time
    statuso = nf90_put_var(ncido,TimeVarido,tday,(/iout/))
    IF (statuso /= nf90_NoErr) CALL handle_err(statuso)

    tday_prev = tday       ! store last time stamp

  ENDIF  ! (tday <= tday_prev)

  ENDDO  ! timesteps

ENDDO  ! WHILE reading listfile

999 CONTINUE

statuso = nf90_sync(ncido)
statuso = nf90_close(ncido)

END PROGRAM roms2z

