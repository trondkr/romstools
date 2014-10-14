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
EXECS   = MakeRivers_Nordsjoen
#
#
all: $(EXECS)

MakeRivers_Nordsjoen: MakeRivers_Nordsjoen.o
	$(FC) $(FFLAGS) -o $@ MakeRivers_Nordsjoen.f90 -I$(INCLUDE) -L$(LIBPATH) -l$(LIBRARY)

MakeRivers_Nordsjoen.o: MakeRivers_Nordsjoen.f90

clean:
	@for name in $(EXECS) ; do\
	(rm -f $$name $$name.o); \
	done

