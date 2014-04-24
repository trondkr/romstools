
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


def contourMap(myWOA,myCDF,mylevel,myvar,month,currentIndex,mytype):
   
    """Get the apporporiate grid"""
    tlat      = myCDF.variables["lat"][:]
    tlon      = myCDF.variables["lon"][:]
    
    if mytype!="surfacelongtermmean" and myvar not in ['zeta','ssflux','shflux']:
        mydata=np.squeeze(myCDF.variables[myvar][currentIndex,0,:,:])
    else:
        mydata=np.squeeze(myCDF.variables[myvar][currentIndex,:,:])
        
    if mytype!="trend" or mytype!="surfacelongtermmean":
        if myvar in ["temp","salt"]:
            depthWOA=np.squeeze(myWOA.variables["depth"][:])
            depthindex=0
        
            for c in xrange(len(depthWOA)):
                if (int(depthWOA[c]) == int(mylevel)):
                    depthindex=c
                    print "Found depth level in WOA file: %s %s"%(depthWOA[depthindex],mylevel)
            
            mydataWOA=np.squeeze(myWOA.variables[myvar][month-1,depthindex,:,:])
        else:
            mydataWOA=mydata*0.0
    else:
        mydataWOA=mydata*0.0
    
    """Calculate the difference between REA and WOA"""    
    mydata = mydata - mydataWOA
    
    plt.figure(figsize=(10,10), frameon=False)

    map = Basemap(llcrnrlon=-18.0,
                  llcrnrlat=46.0,
                  urcrnrlon=25.5,
                  urcrnrlat=67.5,
                  resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)
    delta=20
    levels = np.arange(-2,21,0.5)

    if myvar=="temp":
        levels = np.arange(-3,3,0.2)

    if myvar=="salt":
        levels =  np.arange(-2,2,0.1)
        
    if mytype=="trend" or mytype=="surfacelongtermmean":
        print mydata.max(), mydata.min()
        levels = np.arange(mydata.min(),mydata.max(),(mydata.max()-mydata.min())/20.)
    x, y = map(tlon,tlat)

    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()
    map.drawmapboundary()

  #  map.drawmeridians(np.arange(0,360,1))
  #  map.drawparallels(np.arange(0,90,1))

   # mymaskeddata = np.ma.masked_values(mydata,mymask)
    if mytype=="surfacelongtermmean":
        CS1 = map.contourf(x,y,mydata,levels,cmap=cm.get_cmap('RdBu',len(levels)-1) )#,alpha=0.5)
    else:
        CS1 = map.contourf(x,y,mydata,levels,cmap=cm.get_cmap('RdBu_r',len(levels)-1) )#,alpha=0.5)
    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)

    plotfile='figures/diff_rea_woa_'+str(myvar)+'_time_'+str(month)+'_depth_'+str(mylevel)+'.pdf'
    if mytype=="trend":
        plotfile='figures/trend_'+str(myvar)+'_time_depth_'+str(mylevel)+'.pdf'
    if mytype=="surfacelongtermmean":
        plotfile='figures/longtermmean_'+str(myvar)+'_time_depth_'+str(mylevel)+'.pdf'
    
    
    plt.savefig(plotfile, bbox_inches='tight')
   # plt.show()
    plt.close()


"""" ------------------------------------------------------------------
     MAIN
     Trond Kristiansen, 09.07.2013, 11.02.2014
     Trond.Kristiansen@imr.no
     ------------------------------------------------------------------
"""

doc="""This script reads the output from running createAveragesOnHexagon.py"""

woa="nordsjoen_8km_clim_WOAMONTHLY_climatology_zlevels.nc"

myvars=['salt','temp']
mylevels=['0','100','300','800','2000']
#mylevels=['0']
myprefix='northsea_8km'
mystats=['ymonavg']
myvars=['zeta']
mytype='trend'
#mytype='normal'
#mytype='surfacelongtermmean'

if mytype=="normal":
    for mystat in mystats:
        for myvar in myvars:
            for mylevel in mylevels:
                myfile= str(myprefix)+'_'+str(mystat)+'_'+str(myvar)+'_depth_'+str(mylevel)+'.nc'
                print myfile
                myCDF=Dataset(myfile)
                myWOA=Dataset(woa)
                times=(myCDF.variables["time"][:])
                refDate=datetime.datetime(1948, 1, 1, 0, 0, 0)
                currentIndex=0
                for time in times:
                    currentDate=refDate + datetime.timedelta(days=time)
                    print "Plotting time-step: %s - %s" % (time, currentDate)
    
                    contourMap(myWOA,myCDF,mylevel,myvar,currentDate.month,currentIndex,mytype)
                    currentIndex+=1
           
                myCDF.close()

"""Plot the trends"""
if mytype=="surfacelongtermmean":
    
    for myvar in ['ssflux','shflux']:
        myfile= str(myprefix)+'_timmean_'+str(myvar)+'.nc'
        myCDF=Dataset(myfile)
          
        contourMap(myCDF,myCDF,'surface',myvar,0,0,mytype)
               
        myCDF.close()
            
"""Plot the trends"""
if mytype=="trend":
    for myvar in myvars:
        for mylevel in mylevels:
            if myvar in ['ssflux','shflux','zeta']:
                 myfile= str(myprefix)+'_trend_'+str(myvar)+'_trend.nc'
            else:
                myfile= str(myprefix)+'_trend_'+str(myvar)+'_depth_'+str(mylevel)+'_trend.nc'
            print myfile
            myCDF=Dataset(myfile)
          
            contourMap(myCDF,myCDF,mylevel,myvar,0,0,mytype)
               
            myCDF.close()