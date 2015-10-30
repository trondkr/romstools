
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
import matplotlib as mpl
import mpl_util

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 8, 20)
__modified__ = datetime.datetime(2015, 10, 29)
__version__  = "1.0"
__status__   = "Development, 20.8.2012"



def contourMap(myCDF, bathymetry,lonrho,latrho,myvar,survey,Nobs,obs,first):

    print "Plotting values for index %s to %s"%(sum(Nobs[0:survey]),sum(Nobs[0:survey])+obs)

    """Get the apporporiate grid"""
    lat      = myCDF.variables["obs_lat"][sum(Nobs[0:survey]):sum(Nobs[0:survey])+obs]
    lon      = myCDF.variables["obs_lon"][sum(Nobs[0:survey]):sum(Nobs[0:survey])+obs]
    xgrid    = myCDF.variables["obs_Xgrid"][sum(Nobs[0:survey]):sum(Nobs[0:survey])+obs]
    ygrid    = myCDF.variables["obs_Ygrid"][sum(Nobs[0:survey]):sum(Nobs[0:survey])+obs]
    obs_time = int(np.squeeze(myCDF.variables["obs_time"][sum(Nobs[0:survey]):sum(Nobs[0:survey])+1]))
    mydata   = myCDF.variables[myvar][sum(Nobs[0:survey]):sum(Nobs[0:survey])+obs]

    if sum(np.isnan(mydata)) > 0:
        print "DATA contains NAN!!!!!"
    if first is True:
        print "Max and min latitude: %s - %s"%(lat.max(), lat.min())
        print "Max and min longitude: %s - %s\n"%(lon.max(), lon.min())
        first=False

    """ Current time:"""
    refDate=datetime.datetime(1948,1,1,0,0,0)
    currentDate=refDate + datetime.timedelta(days=obs_time)
    print "Plotting timestep: %s"%(currentDate)

    """Start creating figure"""
    plt.figure(figsize=(10,10), frameon=False)

    map = Basemap(llcrnrlon=-18.0,
                  llcrnrlat=46.0,
                  urcrnrlon=25.5,
                  urcrnrlat=67.5,
                  resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)

    x, y = map(lon,lat)
    x2, y2 = map(lonrho,latrho)

    ax1=plt.subplot(1,1,1)
    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()
    map.drawmapboundary()
    print "Min %s and max %s of data"%(np.ma.min(mydata),np.ma.max(mydata))
    levels=np.arange(-2,18,0.5)

    map.scatter(x,y,s=4,c=mydata,cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r), edgecolor='none',vmin=-2, vmax=18)
    c = plt.colorbar(orientation='horizontal')
    c.set_label("Obs. SST")

   # h=bathymetry.flatten()
   #
   # ax1=plt.subplot(2,2,2)
   #
   # map.scatter(x,y,s=4,c=bathymetry/np.ma.max(bathymetry),cmap=cm.get_cmap('jet'), edgecolor='none')
   ## c = plt.colorbar(orientation='horizontal')
   # c.set_label("Bathymetry")
   #
   # ax2=plt.subplot(2,2,3)
   # scatter(xgrid,ygrid,s=4,c=mydata/np.ma.max(mydata),cmap=cm.get_cmap('jet'), edgecolor='none')
   # ax2.axis('tight')
   #
   # ax2=plt.subplot(2,2,4)
   # scatter(xgrid,ygrid,s=4,c=bathymetry/np.ma.max(bathymetry),cmap=cm.get_cmap('jet'), edgecolor='none')
   #
   # ax2.axis('tight')
    plotfile='figures/obs_'+str(myvar)+'_'+str(currentDate)+'.jpeg'

    plt.savefig(plotfile)
   # plt.show()


infile="NS8KM_AVHRR_obsSST_2012_to_2013.nc"
gridfile="/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m.nc"

myvar="obs_value"
first=False

try:
    myCDF=Dataset(infile)
except IOError:
    print 'Could not open file %s'%(infile)

try:
    myGrid=Dataset(gridfile)
except IOError:
    print 'Could not open file %s'%(gridfile)

Nobs=myCDF.variables["Nobs"][:]
bathymetry=(np.squeeze(myGrid.variables["h"][:]))

lonrho=myGrid.variables["lon_rho"][:]
latrho=myGrid.variables["lat_rho"][:]
maskrho=myGrid.variables["mask_rho"][:]

xgrid2=np.arange(0,len(lonrho[:,0]),1)

ygrid2=np.arange(0,len(latrho[:,0]),1)

bathymetry=np.ma.masked_where(maskrho==0,bathymetry)

maxsurvey=Nobs.shape[0]

surveys=np.arange(0,maxsurvey,1)

increment=1; counter=0
for survey,obs in zip(surveys,Nobs):

    if (counter==increment):
        print "Plotting : %s for Nobs: %s and survey %s"%(myvar,obs,survey)

        contourMap(myCDF,bathymetry,lonrho,latrho,myvar,survey,Nobs,obs,first)
        counter=0
    counter+=1