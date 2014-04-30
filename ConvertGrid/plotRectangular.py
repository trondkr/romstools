
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from netCDF4 import Dataset
import datetime
import os

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 8, 20)
__modified__ = datetime.datetime(2014, 2, 11)
__version__  = "1.0"
__status__   = "Development, 20.8.2012, 9.7.2013, 16.7.2013, 11.2.2014"

doc = """This script reads the output of running grid2lonloatZ.py"""


def contourMap(myCDF,myvar,mylev):

    """Get the apporporiate grid"""
    tlat      = myCDF.variables["lat"][:]
    tlon      = myCDF.variables["lon"][:]
    levels    = myCDF.variables["lev"][:]
    times      = myCDF.variables["time"][:]

    mylevindex=0
    for lev in xrange(len(levels)):
        if (float(levels[lev]) == float(mylev)):
            mylevindex=lev
            print "Will plot depth: %s m (index=%s)"%(mylev,lev)

    mydata=np.squeeze(myCDF.variables[myvar][0,mylevindex,:,:])
    refdate=datetime.datetime(1948,1,1)
    current=refdate + datetime.timedelta(days=times[0])
    print "Date of file is: %s"%(current)
    plt.figure(figsize=(10,10), frameon=False)

    map = Basemap(llcrnrlon=-18.0,
                  llcrnrlat=46.0,
                  urcrnrlon=25.5,
                  urcrnrlat=67.5,
                  resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)
    delta=20.
    levels = np.arange(mydata.min(),mydata.max(),(mydata.max()-mydata.min())/delta)
    levels=np.arange(272,288,0.5)
    print "Range of input data: %s - %s"%(mydata.min(), mydata.max())
    print np.shape(tlon), np.shape(tlat), np.shape(mydata)

    tlons, tlats=np.meshgrid(tlon,tlat)
    x,y = map(tlons,tlats)

    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()
    map.drawmapboundary()

  #  map.drawmeridians(np.arange(0,360,1))
  #  map.drawparallels(np.arange(0,90,1))

    CS1 = map.contourf(x,y,mydata,levels,cmap=cm.get_cmap('RdBu_r',len(levels)-1) ,alpha=1.0)
    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)

    if current.month < 10:
        currentdate=str(current.year)+'0'+str(current.month)
    else:
        currentdate=str(current.year)+str(current.month)

    if not os.path.exists('figures'):
        os.makedirs('figures')
    plotfile='figures/rectangular_'+str(myvar)+'_yyyymm_'+str(currentdate)+'lev_'+str(mylev)+'.pdf'
    print plotfile
    plt.savefig(plotfile, bbox_inches='tight')
   # plt.show()
    plt.close()


#------------------------------------------------------------------
#     MAIN
#     Trond Kristiansen, 09.07.2013, 11.02.2014
#     Trond.Kristiansen@imr.no
#     ------------------------------------------------------------------


mylevels=['0','100','300','800','2000']

#mylevels=['0']
myvars=['votemper']


for myvar in myvars:
    for mylev in mylevels:
        myfile= 'test19940315_mm-IMR-MODEL-ROMS-NWS-20140403-fv02.1.nc'
        myCDF=Dataset(myfile)

        contourMap(myCDF,myvar,mylev)

        myCDF.close()

