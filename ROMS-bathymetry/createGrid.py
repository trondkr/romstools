import os
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4

import pyroms
import pyroms_toolbox
from ROMS_bathy_smoother import *

# Grid dimension
Lm = 60
Mm = 80

lon0=98. ; lat0 = 25.
lon1=98. ; lat1 = -10.
lon2 = 125. ; lat2 = -10.
lon3 = 125. ; lat3 = 25.

lon_min = min(lon0, lon1, lon2, lon3)
lon_max = max(lon0, lon1, lon2, lon3)
lon_0 = 0.5*(lon_min + lon_max)

lat_min = min(lat0, lat1, lat2, lat3)
lat_max = max(lat0, lat1, lat2, lat3)
lat_0 = 0.5*(lat_min + lat_max)

map = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,\
  urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0, \
  resolution='i')

lonp = array([lon0, lon1, lon2, lon3])
latp = array([lat0, lat1, lat2, lat3])
beta = array([1, 1, 1, 1])

# Do this if you aren't going to move the grid corners interactively.
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)

# Do this if you are going to use the Boundary Interactor
#map.drawcoastlines()
#xp, yp = map(lonp, latp)
#bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mm+3,Lm+3), proj=map)
#hgrd=bry.grd

lonv, latv = map(hgrd.x_vert, hgrd.y_vert, inverse=True)
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)
for verts in map.coastsegs:
    hgrd.mask_polygon(verts)

# Edit the land mask interactively.
pyroms.grid.edit_mask_mesh(hgrd, proj=map)

#### Use the following to interpolate from etopo2 bathymetry.
# generate the bathy
# read in topo data (on a regular lat/lon grid)
# this topo come with basemap so you should have it on your laptop.
# just update datadir with the appropriate path
# you can get this data from matplolib svn with
# svn co https://matplotlib.svn.sourceforge.net/svnroot/matplotlib/trunk/htdocs/screenshots/data/"
datadir = '/home/frederic/python/basemap-0.99.4/examples/'
topo = np.loadtxt(os.path.join(datadir, 'etopo20data.gz'))
lons = np.loadtxt(os.path.join(datadir, 'etopo20lons.gz'))
lats = np.loadtxt(os.path.join(datadir, 'etopo20lats.gz'))

# shift data so lons go from -180 to 180 instead of 20 to 380.
topo,lons = shiftgrid(180.,topo,lons,start=False)

# keep only the US
topo = topo[270:525, 30:400]
lons = lons[30:400]
lats = lats[270:525]

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

# interpolate new bathymetry
lon, lat = meshgrid(lons, lats)
h = griddata(lon.flat,lat.flat,topo.flat,hgrd.lon_rho,hgrd.lat_rho)

# insure that depth is always deeper than hmin
h = pyroms_toolbox.change(h, '<', hmin, hmin)

# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

#### Use this for a flat bottom at 50 meters.
h = ones((Mm+2, Lm+2)) * 50.

# Either case.
hgrd.h = h

hc = 5.0
theta_b = 0.4
theta_s = 5.0
Tcline = 5
N = 42

vgrd = pyroms.vgrid.s_coordinate(h, theta_b, theta_s, Tcline, N)
grd_name = 'test'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename='grd.nc')

