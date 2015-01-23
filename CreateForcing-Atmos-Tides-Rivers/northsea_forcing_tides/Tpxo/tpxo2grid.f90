program tpxo2grid

  ! Interpolate TPXO data from netcdf-files to a given grid.
  ! It is assumed that TPXO fields are given on a spherical grid.
  ! The program can take both sequential feltfiles and ROMS-grid file
  ! as input.
  ! Written by Ann Kristin Sperrevik, met.no, 2010

  implicit none
  include "netcdf.inc"

  ! TPXO files
  character (len=200)       :: tpxo_grd
  character (len=200)       :: tpxo_h
  character (len=200)       :: tpxo_u

  ! TPXO grid
  integer                                       :: nx, ny, nc, nct       ! dimensions
  real, dimension(:,:), allocatable             :: lat_h, lat_u, lat_v
  real, dimension(:,:), allocatable             :: lon_h, lon_u, lon_v
  real, dimension(:,:), allocatable             :: topo
  character (len=4), dimension(:), allocatable  :: con
  integer, dimension(:,:), allocatable          :: topo_h, topo_u, topo_v

  character (len=150)        :: grid, outfile, datadir

  ! Output grid
  integer                              :: imax, jmax
  real, dimension(:,:), allocatable    :: lat_grid, lon_grid, angle, lon_rho, lat_rho
  integer, dimension(:,:), allocatable :: topo_grid

  ! amplitude and phase variables
  real, dimension(:,:,:), allocatable  :: h_amp_tpxo, h_amp_grid
  real, dimension(:,:,:), allocatable  :: h_pha_tpxo, h_pha_grid
  real, dimension(:,:,:), allocatable  :: u_amp_tpxo, u_amp_grid
  real, dimension(:,:,:), allocatable  :: u_pha_tpxo, u_pha_grid
  real, dimension(:,:,:), allocatable  :: v_amp_tpxo, v_amp_grid
  real, dimension(:,:,:), allocatable  :: v_pha_tpxo, v_pha_grid

  ! variables for sequential file reading
  integer(kind=2), dimension(20)    ::  ident
  integer,dimension(20)             :: identi
  integer                           ::  igtype, ierr
  real, dimension(6)                ::  gparam
  real, parameter                   :: xoff = +1.0
  real, parameter                   :: yoff = +1.0

  ! Variables needed for interpolation
  real, dimension(:,:), allocatable      :: r1_h, r2_h, r3_h, r4_h
  real, dimension(:,:), allocatable      :: r1_u, r2_u, r3_u, r4_u
  real, dimension(:,:), allocatable      :: r1_v, r2_v, r3_v, r4_v
  integer, dimension(:,:), allocatable   :: i1_h, j1_h
  integer, dimension(:,:), allocatable   :: i1_u, j1_u
  integer, dimension(:,:), allocatable   :: i1_v, j1_v

  ! Variables for phase calculations
  integer                                :: count, npos, nneg, nvalue
  real                                   :: abs_sum_pha
  integer, dimension(4)                  :: iindex, jindex
  real,dimension(4)                      :: tmp

  ! Variables for fill-routines
  real,dimension(:,:), allocatable       :: rmask, work, error

  ! misc. netCDF
  integer                                :: ncid_grd, ncid_u, ncid_h, ncid
  integer                                :: imax_id, jmax_id, nc_id, nct_id, xi_id, eta_id
  integer                                :: status, dimid, varid
  integer, dimension(:), allocatable     :: dimensions
  real, parameter                        :: undef=1e36

  ! misc
  integer             :: i, j, n
  real                :: dum
  character (len=5)   :: type
  character (len=100) :: str, str2

  ! Inquire user for input gridfile

  CALL GETARG(1,grid)
  CALL GETARG(2,type)
  CALL GETARG(3,outfile)
  CALL GETARG(4,datadir)

  write (*,'(2A)') "Directory with TPXO datafiles is located: ", TRIM(datadir)
  !!!tpxo_grd = TRIM(datadir) // "/grid_tpxo7.2.nc"
  !!!tpxo_h   = TRIM(datadir) // "/h_tpxo7.2.nc"
  !!!tpxo_u   = TRIM(datadir) // "/u_tpxo7.2.nc"
  tpxo_grd = TRIM(datadir) // "/gridAO_atlas.nc"
  tpxo_h   = TRIM(datadir) // "/hf.AO_2011atlas.nc"
  tpxo_u   = TRIM(datadir) // "/uv.AO_2011atlas.nc"
  write (*,'(4A)') "The following 3 datafiles will be read: ", TRIM(tpxo_grd), TRIM(tpxo_h), TRIM(tpxo_u)

  write (*,'(/,a14,a,a9,a)') "Input file is ", trim(grid), " of type ", trim(type)
  if (verify(trim(type),"seq").eq.0) then  ! input topography is sequential feltfile

     !********************
     ! Get Met.no grid   *
     !********************

     open(12, file=trim(grid),form="unformatted")
     call getflt(100, 12, 0, identi, imax*jmax, topo, igtype,  gparam, 0, ierr)
     !  close(12)

     ! - Check that the 1st field is a depth matrix
     if (identi(6) /= 351) then
        print *, "***ERROR***"
        print *, "First field in ", trim(grid), " is not depth"
        print *, "identi(6) = ", identi(6)
        print *, "tpxo2grid: terminating"
        stop
     end if

     imax = identi(10)
     jmax = identi(11)
     allocate(lon_grid(imax,jmax),lat_grid(imax,jmax))
     allocate(angle(imax,jmax), topo_grid(imax,jmax), topo(imax,jmax ))

     call seq_grid(12,identi, gparam, imax,jmax, lon_grid, lat_grid, angle,topo)
     close(12)
     where (topo .gt. -30000)
        topo=1
     elsewhere
        topo=0
     end where
     topo_grid=ifix(topo)
     deallocate(angle)

     allocate(lon_rho(0:imax+1, 0:jmax+1), lat_rho(0:imax+1, 0:jmax+1))
     allocate(angle(0:imax+1, 0:jmax+1))

     open(12, file=trim(grid),form="unformatted")
     call seq_grid2roms(12,identi, gparam, imax,jmax, lon_rho, lat_rho, angle,topo)
     deallocate(topo)

  else if (verify(trim(type),"nc").eq.0) then ! Input topography is ROMS grid file (netcdf)

     status=nf_open(trim(grid), NF_NOWRITE, ncid); call check_err(status)
     status=nf_inq_dimid(ncid,'xi_rho', dimid); call check_err(status)
     status = nf_inq_dimlen(ncid, dimid, imax); call check_err(status)
     status=nf_inq_dimid(ncid,'eta_rho', dimid); call check_err(status)
     status = nf_inq_dimlen(ncid, dimid, jmax); call check_err(status)

     allocate(lon_rho(imax,jmax),lat_rho(imax,jmax))
     allocate(angle(imax,jmax), topo(imax,jmax))

     status=nf_inq_varid(ncid,"lon_rho",varid); call check_err(status)
     status=nf_get_var_real(ncid,varid,lon_rho); call check_err(status)
     status=nf_inq_varid(ncid,"lat_rho",varid); call check_err(status)
     status=nf_get_var_real(ncid,varid,lat_rho); call check_err(status)
     status=nf_inq_varid(ncid,"angle",varid); call check_err(status)
     status=nf_get_var_real(ncid,varid,angle); call check_err(status)
     status=nf_inq_varid(ncid,"mask_rho",varid); call check_err(status)
     status=nf_get_var_real(ncid,varid,topo); call check_err(status)
     status=nf_close(ncid); call check_err(status)

     allocate(lon_grid(imax-2,jmax-2),lat_grid(imax-2,jmax-2))
     allocate(topo_grid(imax-2,jmax-2))
     lon_grid(:,:)=lon_rho(2:imax-1,2:jmax-1)
     lat_grid(:,:)=lat_rho(2:imax-1,2:jmax-1)
     topo_grid(:,:)=aint(topo(2:imax-1,2:jmax-1))
     deallocate(topo)
     imax=imax-2; jmax=jmax-2


  else
   print *, "***ERROR***"
   print *, "File format of input file must be either "
   print *, "sequential feltfile (seq) or netCDF (nc) "
   print *, "tpxo2grid: terminating"
   stop
  end if

  write (*, '(/,a29,a,a4)' ) "Grid dimensions on inputfile ", trim(grid), " are"
  write (*, '(a7,i5, x, a11, i5,/)' ) "imax = ", imax, "and jmax = ", jmax

  !*****************
  !* Get TPXO grid *
  !*****************
  status=nf_open(tpxo_grd, NF_NOWRITE, ncid_grd); call check_err(status)
  status=nf_inq_dimid(ncid_grd,'nx', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid_grd, dimid, nx); call check_err(status)
  status=nf_inq_dimid(ncid_grd,'ny', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid_grd, dimid, ny); call check_err(status)
  status=nf_close(ncid_grd); call check_err(status)

  status=nf_open(tpxo_h, NF_NOWRITE, ncid_h); call check_err(status)
  status=nf_inq_dimid(ncid_h,'nc', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid_h, dimid, nc); call check_err(status)
  status=nf_inq_dimid(ncid_h,'nct', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid_h, dimid, nct); call check_err(status)

  allocate(con(nc))
  status=nf_inq_varid(ncid_h,"con",varid); call check_err(status)
  status=nf_get_var_text(ncid_h,varid,con); call check_err(status)
  status=nf_close(ncid_h); call check_err(status)

  write (*, '(a34,a)' ) "Grid dimensions on tpxo inputfile ", trim(tpxo_grd)
  write (*, '(a5,i5, x, a9, i5,/)' ) "nx = ", nx, "and ny = ", ny


  write (*, '(a)') "Tidal constituents on tpxo inputfile:"
  print *, con
  write (*,*)


  ! allocate tpxo grid arrays
  allocate(lon_h(ny,nx), lat_h(ny,nx))
  allocate(lon_u(ny,nx), lat_u(ny,nx))
  allocate(lon_v(ny,nx), lat_v(ny,nx))
  allocate(topo_h(ny,nx), topo_u(ny,nx), topo_v(ny,nx), topo(ny,nx))

  ! Read TPXO grid and bathymetry
  call read_tpxo(lon_h,dum,"lon_z",nx,ny,nc,2,tpxo_grd,tpxo_grd,'z')
  call read_tpxo(lat_h,dum,"lat_z",nx,ny,nc,2,tpxo_grd,tpxo_grd,'z')
  call read_tpxo(topo,dum,"hz   ",nx,ny,nc,2,tpxo_grd,tpxo_grd,'z')

  where (topo .ne. 0.)
     topo=1
  end where
  topo_h=ifix(topo)

  call read_tpxo(lon_u,dum,"lon_u",nx,ny,nc,2,tpxo_grd,tpxo_grd,'u')
  call read_tpxo(lat_u,dum,"lat_u",nx,ny,nc,2,tpxo_grd,tpxo_grd,'u')
  call read_tpxo(topo,dum,"hu   ",nx,ny,nc,2,tpxo_grd,tpxo_grd,'u')

  where (topo .ne. 0.)
     topo=1
  end where
  topo_u=ifix(topo)


  call read_tpxo(lon_v,dum,"lon_v",nx,ny,nc,2,tpxo_grd,tpxo_grd,'u')
  call read_tpxo(lat_v,dum,"lat_v",nx,ny,nc,2,tpxo_grd,tpxo_grd,'u')
  call read_tpxo(topo,dum,"hv   ",nx,ny,nc,2,tpxo_grd,tpxo_grd,'u')
  where (topo .ne. 0.)
     topo=1
  end where
  topo_v=ifix(topo)


  !*******************************
  !*  Set interpolation weights  *
  !*******************************
  ! Note that this part of the code is hard coded for a spherical TPXO grid, with fixed dx, dy
  allocate(i1_h(imax,jmax), j1_h(imax,jmax))
  allocate(r1_h(imax,jmax), r2_h(imax,jmax))
  allocate(r3_h(imax,jmax), r4_h(imax,jmax))
  allocate(i1_u(imax,jmax), j1_u(imax,jmax))
  allocate(r1_u(imax,jmax), r2_u(imax,jmax))
  allocate(r3_u(imax,jmax), r4_u(imax,jmax))
  allocate(i1_v(imax,jmax), j1_v(imax,jmax))
  allocate(r1_v(imax,jmax), r2_v(imax,jmax))
  allocate(r3_v(imax,jmax), r4_v(imax,jmax))

  write (*, '(/,a,/)') "Calculate interpolation weights"

  call interpweights(lon_h, lat_h,nx,ny, lon_grid,lat_grid,imax, jmax, r1_h, r2_h, r3_h, r4_h, i1_h, j1_h,topo_h)
  call interpweights(lon_u, lat_u,nx,ny, lon_grid,lat_grid,imax, jmax, r1_u, r2_u, r3_u, r4_u, i1_u, j1_u,topo_u)
  call interpweights(lon_v, lat_v,nx,ny, lon_grid,lat_grid,imax, jmax, r1_v, r2_v, r3_v, r4_v, i1_v, j1_v,topo_v)


  !****************************************************
  !*  Get Amplitude and phase for surface elevation   *
  !****************************************************

  allocate(h_amp_tpxo(ny,nx,nc), h_amp_grid(imax,jmax,nc))
  allocate(h_pha_tpxo(ny,nx,nc), h_pha_grid(imax,jmax,nc))
  allocate(u_amp_tpxo(ny,nx,nc), u_amp_grid(imax,jmax,nc))
  allocate(u_pha_tpxo(ny,nx,nc), u_pha_grid(imax,jmax,nc))
  allocate(v_amp_tpxo(ny,nx,nc), v_amp_grid(imax,jmax,nc))
  allocate(v_pha_tpxo(ny,nx,nc), v_pha_grid(imax,jmax,nc))


  call read_tpxo(dum,h_amp_tpxo,"ha   ",nx,ny,nc,3,tpxo_h,tpxo_grd,'z')
  call read_tpxo(dum,h_pha_tpxo,"hp   ",nx,ny,nc,3,tpxo_h,tpxo_grd,'z')

  call read_tpxo(dum,u_amp_tpxo,"ua   ",nx,ny,nc,3,tpxo_u,tpxo_grd,'u')
  u_amp_tpxo=u_amp_tpxo/100. ! Convert from cm/s to m/s
  call read_tpxo(dum,u_pha_tpxo,"up   ",nx,ny,nc,3,tpxo_u,tpxo_grd,'u')

  call read_tpxo(dum,v_amp_tpxo,"va   ",nx,ny,nc,3,tpxo_u,tpxo_grd,'v')
  v_amp_tpxo=v_amp_tpxo/100. ! Convert from cm/s to m/s
  call read_tpxo(dum,v_pha_tpxo,"vp   ",nx,ny,nc,3,tpxo_u,tpxo_grd,'v')

  allocate(work(ny,nx),rmask(ny,nx),error(ny,nx))

  write (*, '(/,a )') "Interpolate to output grid"

  ! Interpolate to output grid:
  do n= 1, nc
     ! First fill in landpoints in input fields...
     work=h_amp_tpxo(:,:,n)
     where (topo_h .ne. 1)
        work=undef
     end where
     call creep_fill(work,real(topo_h),1,ny,1,nx,1000,1.0)
     call fill(ny,nx,1,ny,1,nx,work,1.0e35,0.01,1.6,1,rmask,error,nvalue)
     h_amp_tpxo(:,:,n)=work
     work=u_amp_tpxo(:,:,n)
     where (topo_u .ne. 1)
        work=undef
     end where
     call creep_fill(work,real(topo_u),1,ny,1,nx,1000,1.0)
     call fill(ny,nx,1,ny,1,nx,work,1.0e35,0.01,1.6,1,rmask,error,nvalue)
     u_amp_tpxo(:,:,n)=work

     work=v_amp_tpxo(:,:,n)
     where (topo_v .ne. 1)
        work=undef
     end where
     call creep_fill(work,real(topo_v),1,ny,1,nx,1000,1.0)
     call fill(ny,nx,1,ny,1,nx,work,1.0e35,0.01,1.6,1,rmask,error,nvalue)
     v_amp_tpxo(:,:,n)=work

     work=h_pha_tpxo(:,:,n)
     where (topo_h .ne. 1)
        work=undef
     end where
     call creep_fill(work,real(topo_h),1,ny,1,nx,1000,1.0)
     call fill(ny,nx,1,ny,1,nx,work,1.0e35,0.01,1.6,1,rmask,error,nvalue)
     where(work.lt.0.)
        work=work+360.0
     end where
     where(work.ge.360.)
        work=work-360.0
     end where

     where (work.gt.180.)
        work=work-360.
     end where

     h_pha_tpxo(:,:,n)=work
     work=u_pha_tpxo(:,:,n)
     where (topo_u .ne. 1)
        work=undef
     end where
     call creep_fill(work,real(topo_u),1,ny,1,nx,1000,1.0)
     call fill(ny,nx,1,ny,1,nx,work,1.0e35,0.01,1.6,1,rmask,error,nvalue)
     where(work.lt.0.)
        work=work+360.0
     end where
     where(work.ge.360.)
        work=work-360.0
     end where
     where (work.gt.180.)
        work=work-360.
     end where

     u_pha_tpxo(:,:,n)=work

     work=v_pha_tpxo(:,:,n)
     where (topo_v .ne. 1)
        work=undef
     end where
     call creep_fill(work,real(topo_v),1,ny,1,nx,1000,1.0)
     call fill(ny,nx,1,ny,1,nx,work,1.0e35,0.01,1.6,1,rmask,error,nvalue)
     where(work.lt.0.)
        work=work+360.0
     end where
     where(work.ge.360.)
        work=work-360.0
     end where
     where (work.gt.180.)
        work=work-360.
     end where
     v_pha_tpxo(:,:,n)=work
     ! Interpolate tidal data to output grid


     do  j = 1, jmax
        do i = 1, imax

           h_amp_grid(i,j,n) = r1_h(i,j) * h_amp_tpxo(j1_h(i,j)    , i1_h(i,j)    , n) + &
                r2_h(i,j) * h_amp_tpxo(j1_h(i,j)    , i1_h(i,j) + 1, n) + &
                r3_h(i,j) * h_amp_tpxo(j1_h(i,j) + 1, i1_h(i,j) + 1, n) + &
                r4_h(i,j) * h_amp_tpxo(j1_h(i,j) + 1, i1_h(i,j)    , n)

           npos=0; nneg=0; abs_sum_pha=0.0
           iindex(1)=i1_h(i,j); iindex(2)=i1_h(i,j)+1; iindex(3)=iindex(2); iindex(4)=iindex(1)
           jindex(1)=j1_h(i,j); jindex(2)=j1_h(i,j); jindex(3)=j1_h(i,j)+1; jindex(4)=jindex(3)
           tmp=0.0
           do count = 1,4
              if (h_pha_tpxo(jindex(count),iindex(count),n).ge.0.0)then
                 npos=npos+1;
                 tmp(count)=0.0
              else
                 nneg=nneg+1
                 tmp(count)=360.0
              end if
              abs_sum_pha=abs_sum_pha+abs(h_pha_tpxo(jindex(count),iindex(count),n))
           end do
           if (npos.eq.0 .or. nneg.eq.0 .or. (abs_sum_pha/count).lt.90.0) then

              h_pha_grid(i,j,n) = r1_h(i,j) * h_pha_tpxo(j1_h(i,j)    , i1_h(i,j)    , n) + &
                   r2_h(i,j) * h_pha_tpxo(j1_h(i,j)    , i1_h(i,j) + 1, n) + &
                   r3_h(i,j) * h_pha_tpxo(j1_h(i,j) + 1, i1_h(i,j) + 1, n) + &
                   r4_h(i,j) * h_pha_tpxo(j1_h(i,j) + 1, i1_h(i,j)    , n)
           else
              h_pha_grid(i,j,n) = r1_h(i,j) * ( h_pha_tpxo(j1_h(i,j)    , i1_h(i,j)    , n) + tmp(1) ) + &
                   r2_h(i,j) * ( h_pha_tpxo(j1_h(i,j)    , i1_h(i,j) + 1, n) + tmp(2) ) + &
                   r3_h(i,j) * ( h_pha_tpxo(j1_h(i,j) + 1, i1_h(i,j) + 1, n) + tmp(3) ) + &
                   r4_h(i,j) * ( h_pha_tpxo(j1_h(i,j) + 1, i1_h(i,j)    , n) + tmp(4) )

           end if

           if (r1_h(i,j).eq.0.0 .and. r2_h(i,j).eq.0.0 .and. r3_h(i,j).eq.0.0 .and. r4_h(i,j).eq.0.0 ) then
              h_pha_grid(i,j,n)=1.0e35
              h_amp_grid(i,j,n)=1.0e35
           end if

           u_amp_grid(i,j,n) = r1_u(i,j) * u_amp_tpxo(j1_u(i,j)    , i1_u(i,j)    , n) + &
                r2_u(i,j) * u_amp_tpxo(j1_u(i,j)    , i1_u(i,j) + 1, n) + &
                r3_u(i,j) * u_amp_tpxo(j1_u(i,j) + 1, i1_u(i,j) + 1, n) + &
                r4_u(i,j) * u_amp_tpxo(j1_u(i,j) + 1, i1_u(i,j)    , n)

           npos=0; nneg=0; abs_sum_pha=0.0
           iindex(1)=i1_u(i,j); iindex(2)=i1_u(i,j)+1; iindex(3)=iindex(2); iindex(4)=iindex(1)
           jindex(1)=j1_u(i,j); jindex(2)=j1_u(i,j); jindex(3)=j1_u(i,j)+1; jindex(4)=jindex(3)
           tmp=0.0
           do count = 1,4
              if (u_pha_tpxo(jindex(count),iindex(count),n).ge.0.0)then
                 npos=npos+1;
                 tmp(count)=0.0
              else
                 nneg=nneg+1
                 tmp(count)=360.0
              end if
              abs_sum_pha=abs_sum_pha+abs(u_pha_tpxo(jindex(count),iindex(count),n))
           end do
           if (npos.eq.0 .or. nneg.eq.0 .or. (abs_sum_pha/count).lt.90.0) then
              u_pha_grid(i,j,n) = r1_u(i,j) * u_pha_tpxo(j1_u(i,j)    , i1_u(i,j)    , n) + &
                   r2_u(i,j) * u_pha_tpxo(j1_u(i,j)    , i1_u(i,j) + 1, n) + &
                   r3_u(i,j) * u_pha_tpxo(j1_u(i,j) + 1, i1_u(i,j) + 1, n) + &
                   r4_u(i,j) * u_pha_tpxo(j1_u(i,j) + 1, i1_u(i,j)    , n)
           else
              u_pha_grid(i,j,n) = r1_u(i,j) * ( u_pha_tpxo(j1_u(i,j)    , i1_u(i,j)    , n) + tmp(1) ) + &
                   r2_u(i,j) * ( u_pha_tpxo(j1_u(i,j)    , i1_u(i,j) + 1, n) + tmp(2) ) + &
                   r3_u(i,j) * ( u_pha_tpxo(j1_u(i,j) + 1, i1_u(i,j) + 1, n) + tmp(3) ) + &
                   r4_u(i,j) * ( u_pha_tpxo(j1_u(i,j) + 1, i1_u(i,j)    , n) + tmp(4) )

           end if

           v_amp_grid(i,j,n) = r1_v(i,j) * v_amp_tpxo(j1_v(i,j)    , i1_v(i,j)    , n) + &
                r2_v(i,j) * v_amp_tpxo(j1_v(i,j)    , i1_v(i,j) + 1, n) + &
                r3_v(i,j) * v_amp_tpxo(j1_v(i,j) + 1, i1_v(i,j) + 1, n) + &
                r4_v(i,j) * v_amp_tpxo(j1_v(i,j) + 1, i1_v(i,j)    , n)

           npos=0; nneg=0; abs_sum_pha=0.0
           iindex(1)=i1_v(i,j); iindex(2)=i1_v(i,j)+1; iindex(3)=iindex(2); iindex(4)=iindex(1)
           jindex(1)=j1_v(i,j); jindex(2)=j1_v(i,j); jindex(3)=j1_v(i,j)+1; jindex(4)=jindex(3)
           tmp=0.0
           do count = 1,4
              if (v_pha_tpxo(jindex(count),iindex(count),n).ge.0.0)then
                 npos=npos+1;
                 tmp(count)=0.0
              else
                 nneg=nneg+1
                 tmp(count)=360.0
              end if
              abs_sum_pha=abs_sum_pha+abs(v_pha_tpxo(jindex(count),iindex(count),n))
           end do
           if (npos.eq.0 .or. nneg.eq.0 .or. (abs_sum_pha/count).lt.90.0) then

              v_pha_grid(i,j,n) = r1_v(i,j) * v_pha_tpxo(j1_v(i,j)    , i1_v(i,j)    , n) + &
                   r2_v(i,j) * v_pha_tpxo(j1_v(i,j)    , i1_v(i,j) + 1, n) + &
                   r3_v(i,j) * v_pha_tpxo(j1_v(i,j) + 1, i1_v(i,j) + 1, n) + &
                   r4_v(i,j) * v_pha_tpxo(j1_v(i,j) + 1, i1_v(i,j)    , n)
           else
              v_pha_grid(i,j,n) = r1_v(i,j) * ( v_pha_tpxo(j1_v(i,j)    , i1_v(i,j)    , n) + tmp(1) ) + &
                   r2_v(i,j) * ( v_pha_tpxo(j1_v(i,j)    , i1_v(i,j) + 1, n) + tmp(2) ) + &
                   r3_v(i,j) * ( v_pha_tpxo(j1_v(i,j) + 1, i1_v(i,j) + 1, n) + tmp(3) ) + &
                   r4_v(i,j) * ( v_pha_tpxo(j1_v(i,j) + 1, i1_v(i,j)    , n) + tmp(4) )
           end if
        end do
     end do

  end do


  deallocate(work,rmask)
  where(h_pha_grid.lt.0.)
     h_pha_grid=h_pha_grid+360.0
  end where
  where(u_pha_grid.ge.360.)
     u_pha_grid=u_pha_grid-360.0
  end where
  where(u_pha_grid.lt.0.)
     u_pha_grid=u_pha_grid+360.0
  end where
  where(u_pha_grid.ge.360.)
     u_pha_grid=u_pha_grid-360.0
  end where
  where(v_pha_grid.lt.0.)
     v_pha_grid=v_pha_grid+360.0
  end where
  where(v_pha_grid.ge.360.)
     v_pha_grid=v_pha_grid-360.0
  end where
  write (*, '(/,a,/)') "Define output file"

  status = nf_create(trim(outfile), NF_CLOBBER, ncid);  call check_err(status)
  ! --- Define global attributes ---

  status = nf_put_att_text(ncid, NF_GLOBAL, "type", 18, "tidal forcing file")
 call check_err(status)
  status = nf_put_att_text(ncid, NF_GLOBAL, "history",len_trim("Created for "//trim(grid)) , "Created for "//trim(grid))
  call check_err(status)

  status = nf_def_dim(ncid, "imax", imax, imax_id); call check_err(status)
  status = nf_def_dim(ncid, "jmax", jmax, jmax_id); call check_err(status)
  status = nf_def_dim(ncid, "nc", nc, nc_id); call check_err(status)
  status = nf_def_dim(ncid, "nct", nct, nct_id); call check_err(status)
  status = nf_def_dim(ncid, "xi_rho", imax+2, xi_id); call check_err(status)
  status = nf_def_dim(ncid, "eta_rho", jmax+2, eta_id); call check_err(status)
  write (*, '(a)') "Dimensions successfully created"

  allocate(dimensions(2))
  dimensions(1)=nct_id; dimensions(2)=nc_id
  status=nf_def_var(ncid,"con", NF_CHAR,2,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",9,"tidal constituents"); call check_err(status)

  dimensions(1)=imax_id; dimensions(2)=jmax_id
  status = nf_def_var(ncid,"longitude", NF_REAL,2,dimensions,varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",9,"longitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"standard_name",9,"longitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",12,"degrees_east"); call check_err(status)
  status = nf_put_att_real(ncid,varid,"valid_range",NF_REAL,2,real((/ -180, 180 /)));call check_err(status)

  status = nf_def_var(ncid,'latitude',NF_REAL,2,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",8,"latitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"standard_name",8,"latitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",13,"degrees_north"); call check_err(status)
  status = nf_put_att_real(ncid,varid,"valid_range",NF_REAL,2,real((/ -90, 90 /))); call check_err(status)

  status = nf_def_var(ncid,'mask',NF_REAL,2,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",4,"mask"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"option_0",4,"land"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"option_1",5,"water"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",1,"1"); call check_err(status)
  dimensions(1)=xi_id;dimensions(2)=eta_id
  status = nf_def_var(ncid,"lon_rho", NF_REAL,2,dimensions,varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",9,"longitude of RHO-points"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",12,"degree_east"); call check_err(status)
  status = nf_def_var(ncid,'lat_rho',NF_REAL,2,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",8,"latitude of RHO-points"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",13,"degree_north"); call check_err(status)
  status = nf_def_var(ncid,'angle',NF_REAL,2,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",30,"angle between XI-axis and EAST"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",7,"radians"); call check_err(status)
  write (*, '(a)') "Grid variables successfully created"
  deallocate(dimensions)
  allocate(dimensions(3))
  dimensions(1)=imax_id; dimensions(2)=jmax_id; dimensions(3)=nc_id
  status = nf_def_var(ncid,'h_amp',NF_REAL,3,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",25 ,"Tidal elevation amplitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",5, "meter"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"field",17, "amplitude, scalar"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"coordinates",18, "longitude latitude"); call check_err(status)

  status = nf_def_var(ncid,'h_pha',NF_REAL,3,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",21 ,"Tidal elevation phase"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",11, "degree, GMT"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"field",13, "phase, scalar"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"coordinates",18, "longitude latitude"); call check_err(status)

  status = nf_def_var(ncid,'u_amp',NF_REAL,3,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",28 ,"Tidal WE velocity amplitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",14, "meter second-1"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"field",12, "vector, W->E"); call check_err(status)
   status = nf_put_att_text(ncid,varid,"coordinates",18, "longitude latitude"); call check_err(status)

  status = nf_def_var(ncid,'u_pha',NF_REAL,3,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",24 ,"Tidal WE velocity phase"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",11, "degree, GMT"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"field",13, "phase, scalar"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"coordinates",18, "longitude latitude"); call check_err(status)

  status = nf_def_var(ncid,'v_amp',NF_REAL,3,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",28 ,"Tidal SN velocity amplitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",14, "meter second-1"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"field",12, "vector, S->N"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"coordinates",18, "longitude latitude"); call check_err(status)

  status = nf_def_var(ncid,'v_pha',NF_REAL,3,dimensions,varid);  call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",24 ,"Tidal SN velocity phase"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",11, "degree, GMT"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"field",13, "phase, scalar"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"coordinates",18, "longitude latitude"); call check_err(status)
  write (*, '(a)') "Tidal variables successfully created"
  status = nf_enddef(ncid);  call check_err(status)

  do n= 1, nc
     where (topo_grid .eq. 0)
        h_amp_grid(:,:,n)=0
        h_pha_grid(:,:,n)=0
        u_amp_grid(:,:,n)=0
        u_pha_grid(:,:,n)=0
        v_amp_grid(:,:,n)=0
        v_pha_grid(:,:,n)=0
     end where
  end do

  status = nf_inq_varid(ncid,'con',varid); call check_err(status)
  status = nf_put_var_text(ncid,varid,con);  call check_err(status)


  status = nf_inq_varid(ncid,'h_amp',varid); call check_err(status)
  status = nf_put_vara_real(ncid,varid,(/ 1, 1, 1 /), (/ imax,jmax,nc /),h_amp_grid);  call check_err(status)

  status = nf_inq_varid(ncid,'h_pha',varid); call check_err(status)
  status = nf_put_vara_real(ncid,varid,(/ 1, 1, 1 /), (/ imax,jmax,nc /),h_pha_grid);  call check_err(status)

  status = nf_inq_varid(ncid,'u_amp',varid); call check_err(status)
  status = nf_put_vara_real(ncid,varid,(/ 1, 1, 1 /), (/ imax,jmax,nc /),u_amp_grid);  call check_err(status)

  status = nf_inq_varid(ncid,'u_pha',varid); call check_err(status)
  status = nf_put_vara_real(ncid,varid,(/ 1, 1, 1 /), (/ imax,jmax,nc /),u_pha_grid);  call check_err(status)

  status = nf_inq_varid(ncid,'v_amp',varid); call check_err(status)
  status = nf_put_vara_real(ncid,varid,(/ 1, 1, 1 /), (/ imax,jmax,nc /),v_amp_grid);  call check_err(status)

  status = nf_inq_varid(ncid,'v_pha',varid); call check_err(status)
  status = nf_put_vara_real(ncid,varid,(/ 1, 1, 1 /), (/ imax,jmax,nc /),v_pha_grid);  call check_err(status)
  write (*, '(a)') "Tidal variables written to file"

  status = nf_inq_varid(ncid,'mask',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,real(topo_grid));  call check_err(status)

  status = nf_inq_varid(ncid,'angle',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,angle); call check_err(status)

  status = nf_inq_varid(ncid,'longitude',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,lon_grid); call check_err(status)
  status = nf_inq_varid(ncid,'latitude',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,lat_grid); call check_err(status)
  status = nf_inq_varid(ncid,'lon_rho',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,lon_rho); call check_err(status)
  status = nf_inq_varid(ncid,'lat_rho',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,lat_rho); call check_err(status)
  write (*, '(a)') "Grid variables written to file"
  status = nf_close(ncid); call check_err(status)

  write (*, '(/,a,/)') "tpxo2grid: exiting normally"
300 continue
end program tpxo2grid
