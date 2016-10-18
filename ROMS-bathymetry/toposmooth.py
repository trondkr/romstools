import os
from numpy import *
from matplotlib.pyplot import *
from netCDF4 import Dataset
from pylab import *
from mpl_util import LevelColormap
import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import mpl_toolkits.basemap as mp
from bathy_smoother import *
import mpl_util
import laplace_filter
from bathy_smoother import *
import writeGrid

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2013, 7, 10)
__modified__ = datetime.datetime(2015, 3, 12)
__version__  = "1.0"
__status__   = "Development, 10.7.2013, 2.2.2015, 12.3.2015"

doc=""" TOPOSMOOTH.PY
    This script takes a grid created and edited using editNS8KMmask.py or createGrid.py
    and interpolates the Etopo1 bathymetry to the grid and saves as grd.vgrid.h. The data
    is slightly smoothed using Laplacian filter.

"""
def findSubsetIndices(min_lat,max_lat,min_lon,max_lon,lats,lons):

    """Array to store the results returned from the function"""
    res=np.zeros((4),dtype=np.float64)
    minLon=min_lon; maxLon=max_lon

    distances1 = []; distances2 = []
    indices=[]; index=1

    for point in lats:
        s1 = max_lat-point # (vector subtract)
        s2 = min_lat-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    distances1 = []; distances2 = []; index=1

    for point in lons:
        s1 = maxLon-point # (vector subtract)
        s2 = minLon-point # (vector subtract)
        distances1.append((np.dot(s1, s1), point, index))
        distances2.append((np.dot(s2, s2), point, index-1))
        index=index+1

    distances1.sort()
    distances2.sort()
    indices.append(distances1[0])
    indices.append(distances2[0])

    """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
    minJ=indices[1][2]
    maxJ=indices[0][2]
    minI=indices[3][2]
    maxI=indices[2][2]

    res[0]=minI; res[1]=maxI; res[2]=minJ; res[3]=maxJ;
    return res

def createNiceMap(grd,h):
    fig=plt.figure()
    levels = [10,25,50,100,250,500,1000,2500,5000]

    mymap = Basemap(llcrnrlon=-18.0,
                      llcrnrlat=46.0,
                      urcrnrlon=25.5,
                      urcrnrlat=67.5,
                      resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)

   # mymap.drawcoastlines()
   # mymap.drawcountries()
   # mymap.fillcontinents(color='grey')
    mymap.drawmeridians(np.arange(grd.hgrid.lon_rho.min(),grd.hgrid.lon_rho.max(),10),labels=[0,0,0,1])
    mymap.drawparallels(np.arange(grd.hgrid.lat_rho.min(),grd.hgrid.lat_rho.max(),4),labels=[1,0,0,0])

    x,y = mymap(grd.hgrid.lon_rho,grd.hgrid.lat_rho)

    CS1 = mymap.contourf(x,y,h,levels,
                           cmap=mpl_util.LevelColormap(levels,cmap=cm.Blues),
                           extend='upper',
                           alpha=1.0,
                           origin='lower',
                           rasterized=True)
    CS1.axis='tight'
    plotfile='figures/map_KINO.pdf'
    plt.savefig(plotfile, dpi='200')

    #plt.show()

def getEtopo1(lonStart,lonEnd,latStart,latEnd):
    """Get the etopo2 data"""
    etopo1name='/Users/trond/GMT/Etopo/ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lon = etopo1.variables["lon"][:]
    lat = etopo1.variables["lat"][:]
    print lonStart,lonEnd,latStart,latEnd
    res = findSubsetIndices(latStart-5,latEnd+5,lonStart-40,lonEnd+10,lat,lon)

    lonsEtopo1,latsEtopo1=np.meshgrid(lon[res[0]:res[1]],lat[res[2]:res[3]])
    bathyEtopo1 = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    print "Extracted data for area: (%s,%s) to (%s,%s)"%(lonsEtopo1.min(),latsEtopo1.min(),lonsEtopo1.max(),latsEtopo1.max())

    bathySmoothed = laplace_filter.laplace_filter(bathyEtopo1,M=None)
    etopo1.close()

    return bathySmoothed, lon[res[0]:res[1]], lat[res[2]:res[3]]

""" MAIN """

"""Get the grid file defined in /Users/trondkr/Projects/KINO/map/gridid.txt"""
grd = pyroms.grid.get_ROMS_grid('KINO1600M')

"""Get the Etopo1 raw data"""
print "Extracting Etopo1 data for region"
bathyEtopo1, lonsEtopo1, latsEtopo1 = getEtopo1(grd.hgrid.lon_rho.min(),
                                                grd.hgrid.lon_rho.max(),
                                                grd.hgrid.lat_rho.min(),
                                                grd.hgrid.lat_rho.max())

""" Fix minimum depth in Etopo1 matrix"""
print "Changing the minimum depth in Etopo1 file (remove land values)"
hmin = -10
topo = pyroms_toolbox.change(bathyEtopo1, '>', hmin, hmin)
print "Max and min of topo",np.max(topo),np.min(topo)
"""Ocean values are originally negative in etopo1 file but we switch to positive after having removed the
land values"""
topo=-topo


""" interpolate new bathymetry """
print "Interpolating new bathymetry"

#print topo.shape, len(lonsEtopo1), len(latsEtopo1)
h = mp.interp(topo,
            lonsEtopo1,
            latsEtopo1,
            grd.hgrid.lon_rho,
            grd.hgrid.lat_rho,
            checkbounds=False,
            masked=False, order=1)


""" Insure that depth is always deeper than hmin"""
print "Started changing minimum depth"
h = pyroms_toolbox.change(h, '<', abs(hmin), abs(hmin))
print "Depth range in bathymetry: %4.1f - %4.1f"%(topo.min(), topo.max())

""" Check bathymetry roughness """
print "Checking rougness"
RoughMat = bathy_tools.RoughnessMatrix(grd.vgrid.h, grd.hgrid.mask_rho)
print '1a: Max Roughness value in file is: ', RoughMat.max()

RoughMat = bathy_tools.RoughnessMatrix(h, grd.hgrid.mask_rho)
print '1b: Max Roughness value in new bathymetry is: ', RoughMat.max()

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.12
h = bathy_smoothing.smoothing_Positive_rx0(grd.hgrid.mask_rho, h, rx0_max)
RoughMat = bathy_tools.RoughnessMatrix(h, grd.hgrid.mask_rho)
print '2: Max Roughness value in new bathymetry after smoothing is: ', RoughMat.max()
print "Depth range in bathymetry after smoothing: %4.1f - %4.1f"%(h.min(), h.max())

replace = np.where(grd.hgrid.mask_rho == 0)

h[replace] = 0
hgrd = grd.hgrid
# Make sure we remove fjord that have been removed in landmask
grd.vgrid.h = h #*grd.hgrid.mask_rho

theta_b = 0.1
theta_s = 7.0
Tcline = 250
N = 40
grd_name = 'KINO1600M'

vgrd = pyroms.vgrid.s_coordinate(h, theta_b, theta_s, Tcline, N, hraw=grd.vgrid.h)
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

f=Dataset("kino_1600m_30072015_vf20.nc")
visc_factor = f.variables["visc_factor"][:]
diff_factor = f.variables["diff_factor"][:]

print "Range visc_factor %s - %s"%(np.min(visc_factor),np.max(visc_factor))
print "Range diff_factor %s - %s"%(np.min(diff_factor),np.max(diff_factor))

writeGrid.write_ROMS_grid(grd, visc_factor, diff_factor, filename='kino_1600m_02122015_vf20.nc')
#pyroms.grid.write_ROMS_grid(grd, filename='kino_1600m_30072015_vf20_v2.nc')


""" Plot the interpolated bathymetry and the land mask"""
print "Plotting results"
#fig=plt.figure()
#ax = fig.add_subplot(2,1,1)
#pcolor(visc_factor)
#axis('tight')
#ax = fig.add_subplot(2,1,2)
#pcolor(diff_factor)
#axis('tight')
#show()

createNiceMap(grd,grd.vgrid.h)
