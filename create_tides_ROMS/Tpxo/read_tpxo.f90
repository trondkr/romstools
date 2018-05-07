subroutine read_tpxo(F,FF,name,nx,ny,nc,ndim,infile,gridfile,type)

  ! Subroutine to read tpxo data, and make sure that grid i centered on greenwich
  implicit none
  include "netcdf.inc"

  real, dimension(ny,nx)    :: F
  real, dimension(ny,nx,nc) :: FF
  integer                   :: nx, ny, nc,ndim
  character (len=5)         :: name
  character (len=100)        :: gridfile
  character (len=100)        :: infile
  character (len=1)         :: type   ! Must be either z,u,v

  real, dimension(ny,nx)    :: G, lat, lon
  real, dimension(ny,nx,nc) :: GG
  integer                   :: ncid=0,varid=0, status=0
  character (len=5)         :: clon, clat
  integer                   :: i=0,j=0,n=0

  integer                   :: iszero=0, zeroind=0, zeroind2=0, i2=0
  integer                   :: oldind=0, newind=0
  real                      :: sum=0, min=0
  clon="lon_"//type
  clat="lat_"//type

  write(*, '(a8,a5,a6,a20)')  "Extract ", trim(name), " from ", trim(infile)

  status=nf_open(gridfile, NF_NOWRITE, ncid); call check_err(status)
  status=nf_inq_varid(ncid,trim(clat),varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,lat); call check_err(status)

  status=nf_inq_varid(ncid,trim(clon),varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,lon); call check_err(status)

  status=nf_close(ncid); call check_err(status)
  where (lon.gt.180.) 
     lon = lon-360
  end where

  ! Check for longitude 0
  iszero=0
  do i=1,nx
     if (lon(1,i).eq.0.0) then
        iszero=1
        zeroind=i
        exit
     end if
  end do
  min=1.0e37
  if (iszero.eq.0) then

     do i=1,nx
        i2 = i+1
        if (i2.gt.nx) i2=1

        if (abs(lon(1,i)-lon(1,i2)).gt.180) then
           sum=(lon(1,i)+lon(1,i2)+360.)/2.
        else
           sum=(lon(1,i)+lon(1,i2))/2.
        end if
        if (sum.lt.min) then
           min=sum
           iszero=1
           zeroind=i; zeroind2=i2
        end if

        exit
     end do
  end if

  oldind=zeroind
  newind=nint(real(nx)/2.)+1
  if  (ndim.eq.2) then
     status=nf_open(infile, NF_NOWRITE, ncid); call check_err(status)
     status=nf_inq_varid(ncid,trim(name),varid); call check_err(status)
     status=nf_get_var_real(ncid,varid,G); call check_err(status)
     status=nf_close(ncid); call check_err(status)
     do j=1,ny
        do i=1,nx
           if ( newind .gt. nx ) newind=1
           if ( oldind .gt. nx ) oldind=1

           F(j,newind)=G(j,oldind)
           newind=newind+1
           oldind=oldind+1
        end do
     end do
  elseif ( ndim.eq.3) then
     status=nf_open(infile, NF_NOWRITE, ncid); call check_err(status)
     status=nf_inq_varid(ncid,trim(name),varid); call check_err(status)
     status=nf_get_var_real(ncid,varid,GG); call check_err(status)
     status=nf_close(ncid); call check_err(status)
     do  n= 1,nc
        do j=1,ny
           do i=1,nx
              if ( newind .gt. nx ) newind=1
              if ( oldind .gt. nx ) oldind=1

              FF(j,newind,n)=GG(j,oldind,n)
              newind=newind+1
              oldind=oldind+1
           end do
        end do
     end do
  end if






end subroutine read_tpxo
