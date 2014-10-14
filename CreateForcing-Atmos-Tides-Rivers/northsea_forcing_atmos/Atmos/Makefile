#                                 Makefile for producing ROMS forcing
FC = ftn
FFLAGS = -O
LIBPATH = $(NETCDF_DIR)/lib
#
INCLUDE = $(NETCDF_DIR)/include
#
LIBRARY = netcdf
#
TARGETS = all clean
EXECS   = get_forcing 
#
#
all: $(EXECS)

get_forcing: get_forcing.o
	$(FC) $(FFLAGS) -o $@ get_forcing.f90 handle_err.f90 fill.f90 spec_hum.f90 jd.f smooth_polar.f90 -I$(INCLUDE) -L$(LIBPATH) -l$(LIBRARY)

get_forcing.o: get_forcing.f90
handle_err.o: handle_err.f90
fill.o: fill.f90
spec_hum.o: spec_hum.f90
jd.o: jd.f
smooth_polar.o: smooth_polar.f90

clean:
	@for name in $(EXECS) ; do\
	(rm -f $$name $$name.o); \
	done

