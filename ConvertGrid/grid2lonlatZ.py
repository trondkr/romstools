# -*- coding: utf-8 -*-

# Prototype program for interpolasjon til lengde/bredde grid
# Kan leke med område med lon0, lon1 og oppløsning med lonstep
#
# Interpolerer seg litt inn over land.
# Dette kan forbedres ved å velge et output lengde/bredde gitter
# og lage manuelt en skikkelig landmaske.
#
# Tar bare overflatetemperatur.
# For produksjon bør også konvertere fra sigma-koordinater
# til faste dyp (hvilke?)
#
# Bjørn Å., 2008-09-10
#
# Modifisert BÅ 2010-04-07
# Fikset CF-kompatibilitet


import numpy as np
from mathgrid import *
import os
from os import listdir
from os.path import isfile, join
from netCDF4 import Dataset, MFDataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *


# -----------
# Functions
# -----------

def get_ij(longitude, latitude, lons, lats):
    """
    i, j = get_ij(longitude, latitude, grd)

    return the index of the closest point on the grid from the
    point (longitude,latitude) in degree
    """

    lons = lons[:,:] - longitude
    lats = lats[:,:] - latitude

    diff = (lons * lons) + (lats * lats)

    jindex, iindex = np.where(diff==diff.min())

    return iindex[0], jindex[0]


def sample2D(F, X, Y):
    """Bilinar sample of a 2D field

F = 2D array
X, Y = position in grid coordinates, scalars or compatible arrays

Note reversed axes, for integers i and j we have
    sample2D(F, i, j) = F[j,i]

Using linear interpolation
"""

    Z = np.add(X,Y)     # Test for compatibility
    if np.isscalar(Z):  # Both X and Y are scalars
        I = int(X)
        J = int(Y)
        P = X - I
        Q = Y - J
    else:
        # Make arrays of common shape
        X0 = X + np.zeros_like(Z)
        Y0 = Y + np.zeros_like(Z)
        I = X0.astype('int')
        J = Y0.astype('int')
        P = X0 - I
        Q = Y0 - J

    W00 = (1-P)*(1-Q)
    W01 = (1-P)*Q
    W10 = P*(1-Q)
    W11 = P*Q

    jmax, imax = F.shape

    # Condition for being outside the domain
    outside = np.logical_or(np.logical_or(I < 0, I > imax-2),
                            np.logical_or(J < 0, J > jmax-2))


    # Set all invalid indices to zero
    I = np.where(outside, 0, I)
    J = np.where(outside, 0, J)

    # grid cells on land give zero weigth
    W00 = sea_mask[J,I]     * W00
    W01 = sea_mask[J+1,I]   * W10
    W10 = sea_mask[J,I+1]   * W01
    W11 = sea_mask[J+1,I+1] * W11
    total = W00 + W01 + W10 + W11  # total weigth

    undef = np.logical_or(outside, total <= 0)
    # avoid division by zero
    total = np.where(undef, 0.001, total)
    #print J.shape
    #print I.shape
    #I = np.where(undef, 0, I)
    #J = np.where(undef, 0, J)
    A = np.where(undef, UNDEF,
        (W00*F[J,I] + W01*F[J+1,I] + W10*F[J,I+1] + W11*F[J+1,I+1]) / total)
    A = np.where(A < -1.0E10, UNDEF, A)
    A = np.where(A > 1.0E10, UNDEF, A)
    return A

# ---------------------------------------------
#                   MAIN
# ---------------------------------------------

mypath = "/work/users/trondk/NS8km/DELIVERY2/"
mypath="/Users/trondkr/Projects/is4dvar/grid2lonlat/"
infiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]


infiles = ['19940315_mm-IMR-MODEL-ROMS-NWS-20140430-fv02.1.nc']

# Det rektangulære griddet som data skal interpoleres til
lon0, lon1, lonstep = -12, 13.0, 0.1
lat0, lat1, latstep =  49.0, 63.0, 0.05

UNDEF = -32768
FIRST = True

# Dette kunne tas automatisk fra input fil
# På den annen side kan nå variable lett kuttes ut


# 2D variable, YX
varY = []  # konverterer som test på interpolasjonen

# 2D variable, YX
#varXY = ['lon', 'lat']  # konverterer som test på interpolasjonen
varXY=[]
# 3D variable, TYX
#varTYX =['Elev']
varTYX=['ssflux','shflux','zeta']

# 4D variable
varTZYX = ['vozocrtx','vomecrty','votemper','vosaline']

# Levels
levels = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500, 600,
          700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]

print "\nStarted converting netcdf files in polar stereographic projection to"
print "rectangular grid. Input files need to be on Z-level. Scrip must be run in"
print "directory where inputfiles are stored. Results will be stored in sub-diretory: RESULTS\n"

if not os.path.exists('RESULTS'):
        os.makedirs('RESULTS')

for infile in infiles:
    # Read data
    outfile = "RESULTS/rectangular_"+infile[:-3]+".nc"
    f0 = Dataset(infile)
    print "Opened input file: %s"%(infile)
    print "Results will be stored to: %s"%(outfile)

    # Create lookup-table on first loop of infiles
    if (FIRST):
        temp0 = np.array(f0.variables['votemper'][0,0,:,:])
        lonin=f0.variables["lon"][:]
        latin=f0.variables["lat"][:]

        sea_mask = where(temp0 == UNDEF, 0, 1)

        # Lager her en look-up tabell fra lengde, bredde til grid-koord. Dette oppsettet
        # forutsetter at du kjenner grid koordinatene i longe og breddegrad
        # 1. Definer output grid rektangulært
        print "\nCreating lookup-table for converting geographical coordinates to grid coordinates. Takes some time..."
        lon = np.arange(lon0, lon1+lonstep, lonstep)
        lat = np.arange(lat0, lat1+lonstep, latstep)
        lon2=np.zeros(np.shape(lon))
        lat2=np.zeros(np.shape(lat))
        # Rydd opp slik at matrisene ser bra ut uten for mangoe desimaler
        for l in xrange(len(lon)):
            lon2[l]="%3.3f"%(lon[l])
        for l in xrange(len(lat)):
            lat2[l]="%3.3f"%(lat[l])

        lon=lon2; lat=lat2
        llon, llat = meshgrid(lon, lat)

        # Grid koordinatene til det rejtangulære griddet vil lagres i X, Y
        X=np.zeros((np.shape(llon)))
        Y=np.zeros((np.shape(llon)))
        lx=0; ly=0

        # Loop over alle grid punkter (lengde,breddegrad) og finn tilsvarende
        # grid X,Y koordinater. Disse punktene benyttes i sample2D funskjonen for
        # bilineær interpolasjon
        for lx in xrange(len(llon[:,0])):
            for ly in xrange(len(llon[0,:])):
                X[lx,ly], Y[lx,ly] = get_ij(llon[lx,ly],llat[lx,ly],lonin,latin)
            ly=0

        FIRST = False

    # --------------------------
    # Create output NetCDF file
    # --------------------------
    if os.path.exists(outfile): os.remove(outfile)
    print "Interpolating variable:"
    f1 = Dataset(outfile, mode='w', format='NETCDF3_CLASSIC')

    # Global attributes
    f1.Conventions = "CF-1.6"
    f1.title = "NWS REA North Sea assimilation hindcast archive version 1.0. Monthly-mean (full water column) fields"
    f1.history = "Converted by grid2lonlatZ.py"
    f1.source = "IMR, ROMSv3.5, IS4DVAR, NorthSea-8km reanalysis"
    f1.institution = "Institute of Marine Research, Norway"
    f1.references = "http://www.imr.no"
    f1.product_version = "1.0"
    f1.contact = "servicedesk@myocean.eu.org"
    f1.netcdf_version_id = "netCDF-4"
    # Define dimensions
    f1.createDimension('time', None)
    f1.createDimension('lev', len(levels))
    f1.createDimension('lon', len(lon))
    f1.createDimension('lat', len(lat))

    v = f1.createVariable('time', 'd', ('time',))
    v0 = f0.variables['time']
    v.long_name = 'time'
    v.units = "Days since 1948-01-01 00:00:00"
    v.calendar = "Gregorian"

    ntimes = len(f0.dimensions['time'])
    v[:ntimes] = v0[:]

    v = f1.createVariable('lev', 'd', ('lev',))
    v.long_name = "depth"
    v.units = 'm'
    v.positive = 'down'
    v[:] = levels

    # Define coordinate variables
    v = f1.createVariable('lon', 'd', ('lon',))
    v.long_name = 'longitude'
    v.units = 'degrees_east'
    v[:] = lon

    v = f1.createVariable('lat', 'd', ('lat',))
    v.long_name = 'latitude'
    v.units = 'degrees_north'
    v[:] = lat

    # Data variable
    for var in varXY:
        print "Creating (x,y) variable", var
        v0 = f0.variables[var]
        v1 = f1.createVariable(var, 'f', ('lat', 'lon'), fill_value=UNDEF)
        try:
            v1.long_name = v0.long_name
            v1.units = v0.units
        except:
            pass

        v1._CoordinateAxes = "lat lon "


        # Modification
        if var == 'Topo':
            v1.standard_name = 'sea_floor_depth'


        F0 = np.array(v0[:,:])
        #undef = (F0==-32767)
        F0 = F0 * v0.scale_factor
        F1 = sample2D(F0, X, Y)
        # hadde trengt et valid-range atributt, U og V er ofte < null
        F1[F1<=-10.0] = UNDEF              # ikke helt bra, temperature kan være < 0
        v1[:,:] = F1.astype('float64')


    for var in varY:
        print var
        v0 = f0.variables[var]
        v1 = f1.createVariable(var, 'f', ('lev'),
                               fill_value=UNDEF)
        try:
            v1.long_name = v0.long_name
            v1.units = v0.units
        except:
            pass


    for var in varTYX:
        print var
        v0 = f0.variables[var]
        v1 = f1.createVariable(var, 'f', ('time', 'lat', 'lon'),
                               fill_value=UNDEF)

        if var == 'ssflux':
            v1.standard_name = "surface_salinity_flux"
            v1.field="ssflux, scalar, series"
        elif var == 'shflux':
            v1.standard_name = "surface_heat_flux"
            v1.field="shflux, scalar, series"
        elif var=="zeta":
             v1.standard_name = "sea_surface_height_above_geoid"
             v1.field="zeta, scalar, series"
        else:
            v1.standard_name=v0.standard_name

        v1.long_name = v0.long_name
        v1.units = v0.units
        v1.scale_factor=v0.scale_factor
        v1.add_offset=v0.add_offset
        v1.valid_min=v0.valid_min
        v1.valid_max=v0.valid_max
        v1.missing_value=v0.missing_value
        v1._CoordinateAxes = "time lat lon "

        ntimes = len(f0.variables['time'][:])
        for i in xrange(ntimes):
            F0 = np.array(v0[i,:,:])
            #undef = (F0==-32767)
           # F0 = F0 * v0.scale_factor
            F1 = sample2D(F0, X, Y)
            # hadde trengt et valid-range atributt, U og V er ofte < null
        #    F1[F1<=-10.0] = UNDEF
            v1[i,:,:] = F1.astype('float64')

    for var in varTZYX:
        print var
        v0 = f0.variables[var]
        v1 = f1.createVariable(var, 'f', ('time', 'lev', 'lat', 'lon'),
                               fill_value=UNDEF)

        v1.standard_name=v0.standard_name
        v1.long_name = v0.long_name
        if (var=="vosaline"):
            v1.units = "psu"
        else:
            v1.units = v0.units
        v1.scale_factor=v0.scale_factor
        v1.add_offset=v0.add_offset
        v1.valid_min=v0.valid_min
        v1.valid_max=v0.valid_max
        v1.missing_value=v0.missing_value
        v1._CoordinateAxes = "time lev lat lon "


        # Modifications
        if var == 'vozocrtx':
             v1.standard_name = "eastward_sea_water_velocity"
             v1.field="eastward velocity, scalar, series"
        elif var == 'vomecrty':
             v1.standard_name = "northward_sea_water_velocity"
             v1.field="northward velocity, scalar, series"
        elif var == 'vosaline':
             v1.units = "  "
             v1.standard_name = "sea_water_salinity"
             v1.long_name="salinity"
             v1.field="salinity, scalar, series"
        elif var == 'votemper':
             v1.units = "degree_Celsius"
             v1.standard_name = "sea_water_temperature"
             v1.long_name="temperature"
             v1.field="temperature, scalar, series"

        elif var == 'zeta':
            # v1.units = "degree_Celsius"
             v1.standard_name = "sea_surface_height"
            # v1.long_name="temperature"
             v1.field="ssflux, scalar, series"


        ntimes = len(f0.variables['time'][:])
        for l in xrange(ntimes):
            Fz = np.array(v0[l,:,:,:])
           # F0 = F0 * v0.scale_factor

            # Convert to Z levels
          #  Fz = multi_zslice(F0[::-1,:,:], z_rho, -np.asarray(levels))

            for k in xrange(len(levels)):
                A = Fz[k,:,:]
              #  A[A<=-30.0] = UNDEF
                F1 = sample2D(A, X, Y)
              #  F1[F1<=-30.0] = UNDEF
                v1[l,k,:,:] = F1.astype('float64')


    print "Finished\n"
    f0.close()
    f1.close()


