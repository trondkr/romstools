# Compiler

FC = ftn
FC77 = ftn
FFLAGS = -O

# ---------------

# Source for utilities
USRCF90 = readsup.f90 readln.f90 shiftphase.f90 check_err.f90 fill.f90 new_astro.f90 new_vuf.f90
USRCF = gday.f

# Program sources
PSRC = tidenc2roms.f90

PROGS = tidenc2roms

SRC = $(USRCF90) $(PSRC)

UOBJF90 = $(USRCF90:.f90=.o)
UOBJF = $(USRCF:.f=.o)
UOBJ = $(UOBJF90) $(UOBJF)

OBJF90 = $(SRC:.f90=.o)
OBJF = $(USRCF:.f=.o)
OBJ = $(OBJF90) $(OBJF)

# NetCDF
LIBPATH = $(NETCDF_DIR)/lib
INCLUDE = $(NETCDF_DIR)/include

all: $(PROGS)

.SUFFIXES:;
.SUFFIXES: .f90 .f .o

.f90.o:
	$(FC) -c $(FFLAGS) -I$(INCLUDE) $<

.f.o:
	$(FC77) -c $(FFLAGS) -I$(INCLUDE) $<

$(PROGS): $(OBJ)
	$(FC) $(FFLAGS) -o $(PROGS) $(OBJ) -L$(LIBPATH) -lnetcdf

clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *~
	-rm -f $(PROGS)

love:
	@echo "not war"
