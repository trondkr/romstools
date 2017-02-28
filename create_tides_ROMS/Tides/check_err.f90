  subroutine check_err(status)
    include "netcdf.inc"
  ! NetCDF file operation error handling
    integer, intent(in) :: status
    if (status .ne. NF_NOERR) then
      print *, "***NetCDF error, program terminating"
      print *, "NetCDF error message:"
      print *, "  ", nf_strerror(status)
      stop
    endif
  end subroutine check_err
