
"""
This script is used to debug crash/restart files from
the new NS8km model.

Create the resulting images into an animation using:
ffmpeg -r 3 -sameq -i %03d.jpeg NS8KM_1989to1993.mp4

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from netCDF4 import Dataset
import datetime

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 8, 20)
__modified__ = datetime.datetime(2014, 2, 11)
__version__  = "1.0"
__status__   = "Development, 20.8.2012, 9.7.2013, 16.7.2013, 11.2.2014"

doc ="""This script reads the output of running ROMS and plots defined variables
at z-levels using pyroms.
"""


def contourMap(myCDFU,myCDFV):
   
    """Get the apporporiate grid"""
    tlat      = myCDFU.variables["lat"][:]
    tlon      = myCDFU.variables["lon"][:]
   
    mydataU = np.squeeze(myCDFU.variables['u_eastward'][:]) * 100.
    mydataV = np.squeeze(myCDFV.variables['v_northward'][:]) * 100.
    print "U:",mydataU.min(), mydataU.max()
    print "V:",mydataV.min(), mydataV.max()
    
    
    mydata = np.sqrt(mydataU**2 + mydataV**2)
    
    plt.figure(figsize=(10,10), frameon=False)

    map = Basemap(llcrnrlon=-18.0,
                  llcrnrlat=46.0,
                  urcrnrlon=25.5,
                  urcrnrlat=67.5,
                  resolution='l',projection='tmerc',lon_0=0,lat_0=50,area_thresh=200.)
    delta=20
    
    levels = np.arange(-60,60,10)
    x, y = map(tlon,tlat)

    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()
    map.drawmapboundary()

  #  map.drawmeridians(np.arange(0,360,1))
  #  map.drawparallels(np.arange(0,90,1))

   # mymaskeddata = np.ma.masked_values(mydata,mymask)
    
    CS1 = map.contourf(x,y,mydataV,levels,cmap=cm.get_cmap('RdBu_r',len(levels)-1) )#,alpha=0.5)
    delta=2
    
   # map.quiver(tlon[::delta,::delta],
   #            tlat[::delta,::delta],
   #            mydataU[::delta,::delta],
   #            mydataV[::delta,::delta],
   #            scale=2500,color='#000099',
   #            latlon=True)
    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)

    plotfile='figures/quiver_climatology.png'
    
    plt.savefig(plotfile,bbox_inches='tight')
    plt.show()
    plt.close()


"""" ------------------------------------------------------------------
     MAIN
     Trond Kristiansen, 09.07.2013, 11.02.2014
     Trond.Kristiansen@imr.no
     ------------------------------------------------------------------
"""

doc="""This script reads the output from running createAveragesOnHexagon_part2.py"""

depth='10'
season='all'
myfileU = 'northsea_8km_timmean_u_eastward_depth_'+str(depth)+'_season_'+str(season)+'.nc'
myfileV = 'northsea_8km_timmean_v_northward_depth_'+str(depth)+'_season_'+str(season)+'.nc'
              
myCDFU=Dataset(myfileU)
myCDFV=Dataset(myfileV)

contourMap(myCDFU, myCDFV)
                              
myCDFU.close()
myCDFV.close()           