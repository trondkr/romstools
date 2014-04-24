

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from netCDF4 import Dataset
import datetime
import pandas as pd

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 8, 20)
__modified__ = datetime.datetime(2014, 2, 11)
__version__  = "1.0"
__status__   = "Development, 20.8.2012, 9.7.2013, 16.7.2013, 11.2.2014"


# This script reads the output from running createAveragesOnHexagon_part2.py and plots the
# result as depth integreated time-series. More specifilly, it plots all the integrated timeseries of
# domain average of temperature monthly fields as a function of time in
# different layers: (L1) 0-100m (L2) 100-300m (L3) 300-1000m (L4) 1000-3000m
# (TOT) all water column


def calculateWeights(depths,level1,level2):
    print "Finiding weights for depth levels: %s to %s"%(level1,level2)
    
    allwgts=[]
    alldepths=[]
    totalwgt=abs(int(level2)-int(level1))
    
    for depth in depths:
        if (int(depth)==int(level1)):
            startdepth=depth
            alldepths.append(startdepth)
        elif (int(level1) < int(depth) <= int(level2)):
            wgt = float(int(depth) - int(startdepth))/(totalwgt*1.0)
            allwgts.append(wgt)
            alldepths.append(depth)
            startdepth=depth
            
    print "Weights sum up to: %f"%(sum(np.asarray(allwgts)))
    return np.asarray(allwgts),np.asarray(alldepths)
  

   
def createTimeseries(alldata,times,refDate):
  
    myvalues=[]
    mydates=[]
    for time,myvalue in zip(times,alldata):
       
        currentDate=refDate + datetime.timedelta(days=time)
        
       # if currentDate.year==1996 or (currentDate.year==1995 and currentDate.month==12):
       #     if myvar=="temp":
       #         myvalue+=3
                
        mydates.append(currentDate)
        myvalues.append(myvalue)
        
    ts=pd.Series(myvalues,mydates)

    return ts     
            
         
def plotTimeseries(myaxes,ts,myvar,depthlevels,final,counter):
    
    ts_annual = ts.resample("A")
    ts_quarterly = ts.resample("Q")
    ts_monthly = ts.resample("M")
      
  #  ts_quarterly.plot(style="orange",linewidth=5, label="Seasonal")
  #  ts_annual.plot(style="r",linewidth=6, label="Annual")
  #  ts.plot(style="b", marker='o', linewidth=1,label="Monthly raw",alpha=0.5)
    
    ts_monthly.plot(ax=myaxes,style="b", marker='o', linewidth=4,label="Monthly")
    if myvar=="salt":
      #  if counter==2:
       #      ylabel('Salinity')
        myaxes.set_ylim(34.8,35.2)
  #  legend(loc='best')
    if myvar=="temp":
        
      #  if counter==2:
      #      ylabel('Temperature')
        myaxes.set_ylim(-1.0,14)
        plt.yticks(np.arange(-1,14+1, 4))
    if myvar=="salt":
        plt.yticks(np.arange(34.8,35.2+0.1, 0.2))
    #    ylabel('Salinity')
   
    if final is True:
        plotfile='figures/timeseries_'+str(myvar)+'_alldepths.pdf'
        plt.savefig(plotfile,dpi=300)
        print 'Saved figure file %s\n'%(plotfile)
      #  plt.show()
    
# ------------------------------------------------------------------
#      MAIN
#     Trond Kristiansen, 09.07.2013, 11.02.2014
#     Trond.Kristiansen@imr.no
#     ------------------------------------------------------------------
#
# This script reads the output from running createAveragesOnHexagon_part2.py and plots the
# result as depth integreated time-series. More specifilly, it plots all the integrated timeseries of
# domain average of temperature monthly fields as a function of time in
# different layers: (L1) 0-100m (L2) 100-300m (L3) 300-1000m (L4) 1000-3000m
# (TOT) all water column

myvars=['temp','salt']

mylevels=[0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500,
           600, 700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]

myprefix='northsea_8km'


# Define different layers: (L1) 0-100m (L2) 100-300m (L3) 300-1000m (L4) 1000-3000m
# (TOT) all water column

level1s=[0,  100,  300, 1000,    0]
level2s=[100,300, 1000, 3000, 3000]   
   
for myvar in myvars:
    f, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(len(level1s), sharex=True, sharey=True)
    counter=0
    final=False
    first=True     
    # Calculate depth weighted and area integrated values as a function of time
    for level1,level2,myaxes in zip(level1s,level2s,[ax1,ax2,ax3,ax4,ax5]):
        
        depthlevels="%s - %s"%(level1,level2)
        print "Running: %s"%(depthlevels)
        allwgts,alldepths = calculateWeights(mylevels,level1,level2)
     
        for depth,wgt in zip(alldepths,allwgts):
            myfile= str(myprefix)+'_fldmean_'+str(myvar)+'_depth_'+str(depth)+'.nc'
           
            myCDF=Dataset(myfile)
            times=(myCDF.variables["time"][:])
            refDate=datetime.datetime(1948, 1, 1, 0, 0, 0)
            currentIndex=0
           
            mydata=np.squeeze(myCDF.variables[myvar][:])
           
            mydata=mydata*wgt
            if first is True:
                alldata=np.zeros((np.shape(mydata)))
                first=False
            alldata=alldata+mydata
                
            myCDF.close()
        counter+=1
        
        if counter==len(level1s): final=True
        ts = createTimeseries(alldata,times,refDate)
        plotTimeseries(myaxes,ts,myvar,depthlevels,final,counter)
        first=True