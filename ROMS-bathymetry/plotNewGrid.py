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

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2015, 3, 12)
__modified__ = datetime.datetime(2015, 3, 12)
__version__  = "1.0"
__status__   = "Development, 10.7.2013"

def getSelectedStations():

  # Preset stations in MyOcean grid - these are grid indices:
  xpos=[655,504,339,158]
  ypos=[124,444,112,209]
  return xpos, ypos
  
def ingrid(lon, lat, lon_bd,lat_bd):
    return mpl.mlab.inside_poly(zip(lon, lat),  zip(lon_bd, lat_bd))

def getPolygon(grid_lon,grid_lat):

    lon_bd = np.concatenate((grid_lon[:,0],grid_lon[-1,:],grid_lon[::-1,-1], grid_lon[0,::-1] ))
    lat_bd = np.concatenate((grid_lat[:,0],grid_lat[-1,:],grid_lat[::-1,-1], grid_lat[0,::-1] ))

    """ Save the polygon as array to plot later"""
    polygon_data=np.empty((2,len(lon_bd)))
    for k in xrange(len(lon_bd)):
        polygon_data[0,k]=lon_bd[k]
        polygon_data[1,k]=lat_bd[k]

    return polygon_data

def drawGridBoundaries(ax,map,polygon_data,mycolor):
    print "Plotting grid boundaries"
    polygon_data_xy=[]
    for i in xrange(len(polygon_data[0,:])):
        myx,myy=map(polygon_data[0,i],polygon_data[1,i])
        vertices=[myx,myy]
        polygon_data_xy.append(vertices)

    polygon_data_xy=np.asarray(polygon_data_xy)
    patch = Polygon(array(polygon_data_xy), facecolor='none',
        edgecolor=mycolor, linewidth=4)
    ax.add_patch(patch)

def getGrid(filename):
    cdf=Dataset(filename)
    grid_lon=cdf.variables["lon_rho"][:]
    grid_lat=cdf.variables["lat_rho"][:]
    grid_h=cdf.variables["h"]
    print "Grid dimmensions: %s and %s"%(grid_lon.shape,grid_lat.shape)
    cdf.close()

    grid_lon_min=grid_lon.min()
    grid_lon_max=grid_lon.max()

    grid_lat_min=grid_lat.min()
    grid_lat_max=grid_lat.max()
    print "Grid domain longitude from %s to %s"%(grid_lon_min,grid_lon_max)
    print "Grid domain latitude from %s to %s"%(grid_lat_min,grid_lat_max)
    print "----\n"
    return grid_lon, grid_lat,grid_h


def createNiceMap(grd,h,polygon_data_KINO,polygon_data_NS8KM,plotSelectedStations):
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_subplot(111)
    levels = [10,25,50,100,250,500,1000,2500,5000]

    mymap = Basemap(llcrnrlon=-18.0,
                      llcrnrlat=46.0,
                      urcrnrlon=25.5,
                      urcrnrlat=67.5,
                      resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)

    mymap.drawcoastlines()
    mymap.drawcountries()
    mymap.fillcontinents(color='grey')
    mymap.drawmeridians(np.arange(grd.hgrid.lon_rho.min(),grd.hgrid.lon_rho.max(),10),labels=[0,0,0,1])
    mymap.drawparallels(np.arange(grd.hgrid.lat_rho.min(),grd.hgrid.lat_rho.max(),4),labels=[1,0,0,0])

    x,y = mymap(grd.hgrid.lon_rho,grd.hgrid.lat_rho)
    print np.min(grd.hgrid.lon_rho),np.max(grd.hgrid.lon_rho)
    print np.min(grd.hgrid.lat_rho),np.max(grd.hgrid.lat_rho)
    
    CS1 = mymap.contourf(x,y,h,levels,
                           cmap=mpl_util.LevelColormap(levels,cmap=cm.Blues),
                           extend='upper',
                           alpha=1.0,
                           origin='lower',
                           rasterized=True)
    CS1.axis='tight'
    
    """Draw grid boundaries"""
    drawGridBoundaries(ax,mymap,polygon_data_NS8KM,mycolor="red")
    drawGridBoundaries(ax,mymap,polygon_data_KINO,mycolor="magenta")

    if (plotSelectedStations):

        xpos, ypos = getSelectedStations()
        xloc,yloc = mymap(grd.hgrid.lon_rho[ypos,xpos],grd.hgrid.lat_rho[ypos,xpos])
    
        mymap.plot(xloc,yloc,marker="o", color="red", markersize=10, linewidth=0)

    plotfile='figures/map_NS8KM_and_KINO1600M.pdf'
    plt.savefig(plotfile, dpi='200')

    plt.show()

""" MAIN """

"""Get the grid file defined in /Users/trond/GMT/pyroms-master/pyroms/pyroms/gridid.txt"""
grd = pyroms.grid.get_ROMS_grid('KINO1600M')
polygon_data_NS8KM = getPolygon(grd.hgrid.lon_rho,grd.hgrid.lat_rho)

"""Read the grid info from the grid file"""
filename="kino_1600m_18062015.nc"
lon_rho,lat_rho,grid_h = getGrid(filename)

plotSelectedStations=True

"""Calculate the x,y grid coordinates"""
(Mp,Lp)=lon_rho.shape
X=np.arange(0,Mp,1)
Y=np.arange(0,Lp,1)

roms_Xgrid,roms_Ygrid=np.meshgrid(Y,X)

"""Loop over all times and store to file or make map"""
polygon_data_KINO = getPolygon(lon_rho,lat_rho)


""" Plot the interpolated bathymetry and the land mask"""

#show()
createNiceMap(grd,grd.vgrid.h,polygon_data_KINO,polygon_data_NS8KM,plotSelectedStations)