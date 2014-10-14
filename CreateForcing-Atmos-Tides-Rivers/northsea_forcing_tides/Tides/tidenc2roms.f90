! -----------------------------------------------------------------------
! tide2roms
! 
! Read netcdf files with tidal information
! for a polar stereographic domain and creates a NetCDF tide-file for ROMS
! ----------------------------------------------------------------------

! The program is written in Fortran 90
! It uses standard Fortran 77 NetCDF calls for NetCDF 
! partly directly and partly through the module ncout_util
!
! Written by Ann Kristin Sperrevik, met.no, 2010
!
! ----------------------------------------------------------------------



program tidenc2roms

  implicit none
  include "netcdf.inc"

  ! ------------------------------------------
  ! --- Constants 
  ! ------------------------------------------


  character(len=*), parameter :: outfile = "tide.nc"

  ! --- Constituents
  integer, parameter :: Nconstit = 8
  character(len=5), dimension(Nconstit), parameter :: constit =   &
       (/ 'K2   ', 'S2   ', 'M2   ', 'N2   ', 'K1   ', 'P1   ',   &
       'O1   ', 'Q1   ' /)!,  'MM   ', 'MF   ', 'M4   ' /)

  character(len=*), parameter :: infile =  "tpxo.nc"  

  real, dimension(Nconstit), parameter :: freq =      &
       (/ 0.0835614924, 0.0833333333,                 &
       0.0805114007, 0.0789992488,                 &
       0.0417807462, 0.0415525871,                 &
       0.0387306544, 0.0372185026 /) ! ,             &
  !          0.0015121518, 0.0030500918,                 &
  !          0.1610228013 /)
  real, dimension(Nconstit), parameter :: period =      &
       (/ 1/freq(1), 1/freq(2), 1/freq(3),              &
       1/freq(4), 1/freq(5), 1/freq(6),              &
       1/freq(7), 1/freq(8) /)!, 1/freq(9),           &
  !        1/freq(10), 1/freq(11) /)  

  ! --- NetCDF global title attribute
  character(len=*), parameter :: title = "Tidal forcing from met.no"

  ! --- Debug switch
  logical, parameter :: debug = .TRUE.


  ! --- Mathematical values
  real, parameter :: pi = 3.1415926
  real, parameter :: rad = pi/180.0
  real, parameter :: deg = 180.0/pi
  complex, parameter :: i = (0.0, 1.0)

  ! ---------------------------------------------
  ! --- Variables 
  ! ---------------------------------------------

  ! -- Variables for feltfile reading

  integer :: imax, jmax, nc, nct

  real, dimension(:,:), allocatable :: uamp, uphase
  real, dimension(:,:), allocatable :: vamp, vphase
  complex, dimension(:,:), allocatable :: A, B
  real, dimension(:,:), allocatable :: rA, rB
  real, dimension(:,:), allocatable :: fA, fB

  real, dimension(:,:), allocatable :: rmask, error
  real, dimension(:,:), allocatable :: angle      ! Local rotation angle  

  ! --- Output fields
  real, dimension(:,:,:), allocatable :: h_amp
  real, dimension(:,:,:), allocatable :: h_pha
  real, dimension(:,:,:), allocatable :: u_amp
  real, dimension(:,:,:), allocatable :: u_pha
  real, dimension(:,:,:), allocatable :: v_amp
  real, dimension(:,:,:), allocatable :: v_pha
  real, dimension(:,:), allocatable :: tide_Eamp
  real, dimension(:,:), allocatable :: tide_Ephase
  real, dimension(:,:), allocatable :: tide_Cmax
  real, dimension(:,:), allocatable :: tide_Cmin
  real, dimension(:,:), allocatable :: tide_Cangle
  real, dimension(:,:), allocatable :: tide_Cphase
  real, dimension(:,:), allocatable :: Areal, Aimag
  real, dimension(:,:), allocatable :: Breal, Bimag
  real, dimension(:,:), allocatable :: lat,lon, field, lon_rho, lat_rho
  character (len=4), dimension(:), allocatable :: con 
  ! --- NetCDF dimensions
  integer ::  xi_rho, eta_rho, tide_period

  ! --- misc NetCDF
  integer :: status    ! Error return code
  integer :: ncid      ! NetCDF file identifier
  integer :: varid     ! NetCDF variable identifier
  integer :: dimid
  integer :: imax_id, jmax_id

  ! --- Indices
  integer :: ii, jj
  integer :: k

  ! --- misc variables
  !real :: lon, lat
  !real :: lon0, lat0
  !real :: lon1, lat1
  character(len=200) :: string, constring
  integer, dimension(8) :: date
  integer, dimension(:), allocatable :: indx
  integer :: nvalue
  integer :: testdum
  real :: reflat
  character (len=20) :: str, str2
  ! ==================================================================

  write (*,'(/,a12,a)') "Input file: ", trim(infile)
  ! Get dimensions from input file
  status=nf_open(infile, NF_NOWRITE, ncid); call check_err(status)
  status=nf_inq_dimid(ncid,'imax', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid, dimid, imax); call check_err(status)
  status=nf_inq_dimid(ncid,'jmax', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid, dimid, jmax); call check_err(status)

  status=nf_inq_dimid(ncid,'nc', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid, dimid, nc); call check_err(status)
  status=nf_inq_dimid(ncid,'nct', dimid); call check_err(status)
  status = nf_inq_dimlen(ncid, dimid, nct); call check_err(status)
  !
  ! --- Allocate arrays ---
  !
  allocate(field(imax, jmax),lon(imax, jmax),lat(imax, jmax))
  allocate(h_amp(imax, jmax, nc))
  allocate(h_pha(imax, jmax, nc))
  allocate(u_amp(imax, jmax, nc))
  allocate(u_pha(imax, jmax, nc))
  allocate(v_amp(imax, jmax, nc))
  allocate(v_pha(imax, jmax, nc))
  allocate(tide_Eamp(  0:imax+1, 0:jmax+1))
  allocate(tide_Ephase(0:imax+1, 0:jmax+1))
  allocate(tide_Cmax(  0:imax+1, 0:jmax+1))
  allocate(tide_Cmin(  0:imax+1, 0:jmax+1))
  allocate(tide_Cangle(0:imax+1, 0:jmax+1))
  allocate(tide_Cphase(0:imax+1, 0:jmax+1))
  allocate(uamp(0:imax+1, 0:jmax+1), uphase(0:imax+1, 0:jmax+1))
  allocate(vamp(0:imax+1, 0:jmax+1), vphase(0:imax+1, 0:jmax+1))
  allocate(A(0:imax+1, 0:jmax+1), B(0:imax+1, 0:jmax+1))
  allocate(rA(0:imax+1, 0:jmax+1), rB(0:imax+1, 0:jmax+1))
  allocate(fA(0:imax+1, 0:jmax+1), fB(0:imax+1, 0:jmax+1))
  allocate(angle(0:imax+1,0:jmax+1),lon_rho(0:imax+1,0:jmax+1),lat_rho(0:imax+1,0:jmax+1))
  allocate(rmask(0:imax+1,0:jmax+1), error(0:imax+1,0:jmax+1))


  allocate(con(nc))
  status=nf_inq_varid(ncid,"con",varid); call check_err(status)
  status=nf_get_var_text(ncid,varid,con); call check_err(status)

  status=nf_inq_varid(ncid,"angle",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,field); call check_err(status)

  status=nf_inq_varid(ncid,"longitude",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,lon); call check_err(status)

  status=nf_inq_varid(ncid,"latitude",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,lat); call check_err(status) 

  status=nf_inq_varid(ncid,"lon_rho",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,lon_rho); call check_err(status)

  status=nf_inq_varid(ncid,"lat_rho",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,lat_rho); call check_err(status) 

  status=nf_inq_varid(ncid,"h_amp",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,h_amp); call check_err(status)

  status=nf_inq_varid(ncid,"h_pha",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,h_pha); call check_err(status)

  status=nf_inq_varid(ncid,"u_amp",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,u_amp); call check_err(status)

  status=nf_inq_varid(ncid,"u_pha",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,u_pha); call check_err(status)

  status=nf_inq_varid(ncid,"v_amp",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,v_amp); call check_err(status)

  status=nf_inq_varid(ncid,"v_pha",varid); call check_err(status)
  status=nf_get_var_real(ncid,varid,v_pha); call check_err(status)

  status=nf_close(ncid); call check_err(status)
  write (*, '(a26)' ) "Inputfile grid dimensions:"
  write (*, '(a7,i5, x, a11, i5,/)' ) "imax = ", imax, "and jmax = ", jmax
  write (*, '(a)') "Tidal constituents on inputfile:"
  print *, con

  ! Include only constituents given in constit 
  allocate(indx(nc))
  do k = 1, nc 
     str=con(k)
     call To_Upper(str)
     if (verify(trim(constit(1)),trim(str)).eq.0) indx(1)=k
     if (verify(trim(constit(2)),trim(str)).eq.0) indx(2)=k
     if (verify(trim(constit(3)),trim(str)).eq.0) indx(3)=k
     if (verify(trim(constit(4)),trim(str)).eq.0) indx(4)=k 
     if (verify(trim(constit(5)),trim(str)).eq.0) indx(5)=k
     if (verify(trim(constit(6)),trim(str)).eq.0) indx(6)=k
     if (verify(trim(constit(7)),trim(str)).eq.0) indx(7)=k
     if (verify(trim(constit(8)),trim(str)).eq.0) indx(8)=k
 
!!$     if (verify(trim(constit(9)),trim(str)).eq.0) indx(9)=k 
!!$     if (verify(trim(constit(10)),trim(str)).eq.0) indx(10)=k
!!$     if (verify(trim(constit(11)),trim(str)).eq.0) indx(11)=k
!!$     if (verify(trim(constit(12)),trim(str)).eq.0) indx(12)=k
!!$     if (verify(trim(constit(13)),trim(str)).eq.0) indx(13)=k
  end do
 
  do k=1, nconstit
     string=constit(k)

    read(string,'(a)') str2
     if (len_trim(constring).eq.len(constring)) then
        constring=trim(str2)
     else
        constring = trim(constring)//" "//trim(str2)
     end if
  end do
  write (*, '(/,a)') "Tidal constituents on outputfile:"
  write (*, '(a)') trim(constring)

  ! --------------------------------
  ! ---- Define NetCDF grid file ---
  ! --------------------------------
  write (*, '(/,a,a)') "Create and define output file: ", trim(outfile)
  status = nf_create(trim(outfile), NF_CLOBBER, ncid);  call check_err(status)
  ! --- Define global attributes ---

  status = nf_put_att_text(ncid, NF_GLOBAL, "title", len_trim(title), title)
  call check_err(status)

  status = nf_put_att_text(ncid, NF_GLOBAL, "type", 9, "Tidal forcing file")
  call check_err(status)

  write(string,'(A)') "Created by tide2roms from sequential files from met.no"
  status = nf_put_att_text(ncid, NF_GLOBAL, "history", len_trim(string), string)
  call check_err(status)
 
  status = nf_put_att_text(ncid, NF_GLOBAL, 'constituents',     &
       len_trim(constring), constring)
  call check_err(status)

  call date_and_time(values=date)
  status = nf_put_att_int(ncid, NF_GLOBAL, 'Creation_date', NF_INT, 3, date(1:3))

  status = nf_def_dim(ncid, "xi_rho", imax+2, imax_id); call check_err(status)
  status = nf_def_dim(ncid, "eta_rho", jmax+2, jmax_id); call check_err(status)
  status = nf_def_dim(ncid, "tide_period", Nconstit, tide_period); call check_err(status)

  status = nf_def_var(ncid,"tide_period", NF_REAL,1,(/ tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",12,"Tidal period"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",4,"hour"); call check_err(status)

  status = nf_def_var(ncid,"tide_Eamp", NF_REAL,3,(/ imax_id, jmax_id, tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",25,"Tidal elevation amplitude"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",5,"meter"); call check_err(status)

  status = nf_def_var(ncid,"tide_Ephase", NF_REAL,3,(/ imax_id, jmax_id, tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",22,"Tidal elevation phase"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",7,"degrees"); call check_err(status)

  status = nf_def_var(ncid,"tide_Cmax", NF_REAL,3,(/ imax_id, jmax_id, tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",21,"Maximum tidal current"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",14,"meter second-1"); call check_err(status)

  status = nf_def_var(ncid,"tide_Cmin", NF_REAL,3,(/ imax_id,jmax_id, tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",21,"Minimum tidal current"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",14,"meter second-1"); call check_err(status)

  status = nf_def_var(ncid,"tide_Cangle", NF_REAL,3,(/ imax_id,jmax_id, tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",20,"Tidal ellipsis angle"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",7,"degrees"); call check_err(status)

  status = nf_def_var(ncid,"tide_Cphase", NF_REAL,3,(/ imax_id,jmax_id, tide_period /),varid); call check_err(status)
  status = nf_put_att_text(ncid,varid,"long_name",21,"Tidal current phase"); call check_err(status)
  status = nf_put_att_text(ncid,varid,"units",7,"degrees"); call check_err(status)

  status = nf_enddef(ncid);  call check_err(status)

  status = nf_inq_varid(ncid,'tide_period',varid); call check_err(status)
  status = nf_put_var_real(ncid,varid,period);  call check_err(status)

 
  ! -------------------------------------
  ! --- Loop through the constituents ---
  ! -------------------------------------
  write (*, '(/,a)') "Loop through the constituents: "
  reflat=55;
  do k = 1, Nconstit
   
     write (*,'(/,a4,a13,f7.4,1x,a5)') trim(constit(k)), " with period ", period(k), "hours"
     call extend(h_amp(:,:,indx(k)), tide_Eamp,imax,jmax)
     call extend(h_pha(:,:,indx(k)), tide_Ephase,imax,jmax)
     call extend(u_amp(:,:,indx(k)), uamp,imax,jmax)
     call extend(u_pha(:,:,indx(k)), uphase,imax,jmax)
     call extend(v_amp(:,:,indx(k)), vamp,imax,jmax)
     call extend(v_pha(:,:,indx(k)), vphase,imax,jmax)

     write(*,'(a)') "Perform phase shift and nodal corrections"

     call shiftphase(tide_Eamp,tide_Ephase,uamp, uphase,vamp,vphase,imax,   &
          jmax,constit(k),reflat)

     allocate(Areal(0:imax+1,0:jmax+1), Aimag(0:imax+1,0:jmax+1))
     allocate(Breal(0:imax+1,0:jmax+1), Bimag(0:imax+1,0:jmax+1))

     ! --- Elevation, fill
     ! Eamp*cos(wt-Ephase) = 
     !    Eamp*cos(Ephase)*cos(wt) + Eamp*sin(Ephase)*sin(wt)
  

     Areal = tide_Eamp*cos(tide_Ephase*rad)
     Breal = tide_Eamp*sin(tide_Ephase*rad)

     where (tide_Eamp < 0.001)
        Areal = 1.0E36
        Breal = 1.0E36
     end where
     call fill(imax+2, jmax+2, 1, imax+2, 1, jmax+2, Areal,              &
          1.0E35, 0.01, 1.6, 100, rmask, error, nvalue) 
     call fill(imax+2, jmax+2, 1, imax+2, 1, jmax+2, Breal,              &
          1.0E35, 0.01, 1.6, 100, rmask, error, nvalue) 

     tide_Eamp   = sqrt(Areal*Areal + Breal*Breal)
     tide_Ephase = atan2(Breal, Areal)*deg
     where(tide_Ephase < 0)
        tide_Ephase = tide_Ephase + 360.0
     end where
     write(*,'(a)') "Convert to ellipse parameters"
     ! --- Convert to ellipse parameters
     A = 0.5*(uamp*exp(-i*uphase*rad) + i*vamp*exp(-i*vphase*rad))
     B = 0.5*(uamp*exp( i*uphase*rad) + i*vamp*exp( i*vphase*rad))

     ! --- Fill
     Areal = real(A)
     Aimag = imag(A)
     where (abs(A) < 0.001)
        Areal = 1.0E36
        Aimag = 1.0E36
     end where

     call fill(imax+2, jmax+2, 1, imax+2, 1, jmax+2, Areal,              &
          1.0E35, 0.01, 1.6, 100, rmask, error, nvalue) 
     call fill(imax+2, jmax+2, 1, imax+2, 1, jmax+2, Aimag,              &
          1.0E35, 0.01, 1.6, 100, rmask, error, nvalue) 

     A = Areal + i*Aimag

     Breal = real(B)
     Bimag = imag(B)
     where (abs(B) < 0.001)
        Breal = 1.0E36
        Bimag = 1.0E36
     end where

     call fill(imax+2, jmax+2, 1, imax+2, 1, jmax+2, Breal,              &
          1.0E35, 0.01, 1.6, 100, rmask, error, nvalue) 
     call fill(imax+2, jmax+2, 1, imax+2, 1, jmax+2, Bimag,              &
          1.0E35, 0.01, 1.6, 100, rmask, error, nvalue) 

     B = Breal + i*Bimag
     deallocate(Areal,Breal, Aimag, Bimag)
     rA = abs(A)
     fA = imag(log(A/rA))
     where(rA < 0.00001)
        fA = 0.0
     end where
     rB = abs(B) 
     fB = imag(log(B/rB))
     where(rB < 0.00001)
        fB = 0.0
     end where
     tide_Cmax = rA + rB
     tide_Cmin = rA - rB
     tide_Cangle = 0.5*(fA + fB)*deg
     write(*,'(a)') "Rotate so that Cangle is measured from the grid's east"
     ! rotate so that Cangle is measured from east
     tide_Cangle = tide_Cangle + angle

     where(tide_Cangle < 0)
        tide_Cangle = tide_Cangle + 360.0
     end where
     tide_Cphase = 0.5*(fB - fA)*deg
     where(tide_Cphase < 0)
        tide_Cphase = tide_Cphase + 360.0
     end where
     write(*,'(a)') "Write tidal forcing to file"
     status = nf_inq_varid(ncid, "tide_Eamp", varid);  call check_err(status)
 
     status = nf_put_vara_real(ncid, varid, (/ 1, 1, k /), (/ imax+2, jmax+2, 1 /), tide_Eamp(0:imax+1,0:jmax+1));  call check_err(status)
     status = nf_inq_varid(ncid, "tide_Ephase", varid);  call check_err(status)
     status = nf_put_vara_real(ncid, varid, (/ 1, 1, k /), (/ imax+2, jmax+2, 1 /), tide_Ephase(0:imax+1,0:jmax+1));  call check_err(status)
     status = nf_inq_varid(ncid, "tide_Cmax", varid);  call check_err(status)
     status = nf_put_vara_real(ncid, varid, (/ 1, 1, k /), (/ imax+2, jmax+2, 1 /), tide_Cmax(0:imax+1,0:jmax+1));  call check_err(status)
     status = nf_inq_varid(ncid, "tide_Cmin", varid);  call check_err(status)
     status = nf_put_vara_real(ncid, varid, (/ 1, 1, k /), (/ imax+2, jmax+2, 1 /), tide_Cmin(0:imax+1,0:jmax+1));  call check_err(status)

     status = nf_inq_varid(ncid, "tide_Cangle", varid);  call check_err(status)
     status = nf_put_vara_real(ncid, varid, (/ 1, 1, k /), (/ imax+2, jmax+2, 1 /), tide_Cangle(0:imax+1,0:jmax+1));  call check_err(status)
     status = nf_inq_varid(ncid, "tide_Cphase", varid);  call check_err(status)
     status = nf_put_vara_real(ncid, varid, (/ 1, 1, k /), (/ imax+2, jmax+2, 1 /), tide_Cphase(0:imax+1,0:jmax+1));  call check_err(status)

  end do ! k

  status = nf_close(ncid);  call check_err(status)
  write (*, '(/,a,/)') "tidenc2roms: exiting normally"
  ! -------------------------
  ! --- Internal routines ---
  ! -------------------------

end program tidenc2roms

  subroutine extend(A,B,imax,jmax)
    real, dimension(imax,jmax), intent(in) :: A
    real, dimension(0:size(A,1)+1, 0:size(A,2)+1), intent(out) :: B
    integer :: imax, jmax

    B(1:imax, 1:jmax) = A
    B(1:imax, 0)      = A(:,1)
    B(1:imax, jmax+1) = A(:,jmax)
    B(0, 1:jmax)      = A(1,:)
    B(imax+1, 1:jmax) = A(imax,:)
    B(0,0)            = A(1,1)
    B(0,jmax+1)       = A(1,jmax)
    B(imax+1,0)       = A(imax,1)
    B(imax+1,jmax+1)  = A(imax,jmax)
  end subroutine extend


  subroutine To_upper(str)
    character(*), intent(in out) :: str
    integer :: i

    do i = 1, len(str)
       select case(str(i:i))
       case("a":"z")
          str(i:i) = achar(iachar(str(i:i))-32)
       end select
    end do
  end subroutine To_upper

  subroutine To_lower(str)
    character(*), intent(in out) :: str
    integer :: i

    do i = 1, len(str)
       select case(str(i:i))
       case("A":"Z")
          str(i:i) = achar(iachar(str(i:i))+32)
       end select
    end do
  end subroutine To_Lower

