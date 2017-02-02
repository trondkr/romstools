subroutine handle_err(status)
use netcdf
integer status
write(*,*) 'Error on netCDF, errno = ',status
print *, trim(nf90_strerror(status))
stop "stopped"
end
