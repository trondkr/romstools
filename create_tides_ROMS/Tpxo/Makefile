#
#                                 Makefile for producing ROMS forcing
FC = gfortran
FFLAGS = -O3 
#
NC_CONFIG ?= nc-config
LIBPATH =  /usr/local/lib
#
INCLUDE =  /usr/local/include
#
LIBRARY = netcdff
#
TARGETS = all clean
EXECS   = tpxo2grid
#
#
all: $(EXECS)

tpxo2grid: tpxo2grid.o
	$(FC) $(FFLAGS) -o $@ tpxo2grid.f90 -O3 inout.f90 check_err.f90 fill.f90 fill_mod.f90 creep_fill.f90 interpweights.f90 seq_grid.f90 read_tpxo.f90 seq_grid2roms.f90 -I$(INCLUDE) -L$(LIBPATH) -l$(LIBRARY)

tpxo2grid.o: tpxo2grid.f90
inout.o: inout.f90
check_err.o: check_err.f90
fill.o: fill.f90
fill_mod.o: fill_mod.f90
creep_fill.o: creep_fill.f90
interpweights.o: interpweights.f90
seq_grid.o: seq_grid.f90
read_tpxo.o: read_tpxo.f90
seq_grid2roms.o: seq_grid2roms.f90

clean:
	@for name in $(EXECS) ; do\
	(rm -f $$name $$name.o); \

	done

