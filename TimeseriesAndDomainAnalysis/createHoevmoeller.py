

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
# results as timeseries (Hoevmoeller diagram).


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
        
        mydates.append(currentDate)
        myvalues.append(myvalue)
        
    ts=pd.Series(myvalues,mydates)

    return ts     
            
         
def plotTimeseries(mytimes,myyears,times,mydata,myvar,depthlevels,mytype):
    
    depthlevels=-np.asarray(depthlevels)
    ax = figure().add_subplot(111)
    
    y,x = np.meshgrid(depthlevels,times)
    
    mydata=np.rot90(mydata)
    if myvar=='temp':
        levels = np.arange(-2,16,1)
    if myvar=='salt':
        levels = np.arange(mydata.min(),mydata.max()+0.1,0.05)
    if mytype=="T-CLASS3-IC_CHANGE":
             levels = np.arange(mydata.min(),mydata.max()+0.1,0.1)
    if mytype=="S-CLASS3-IC_CHANGE":
             levels = np.arange(mydata.min(),mydata.max()+0.1,0.05)
             
    print mydata.min(), mydata.max()
    cs=contourf(x,y,mydata,levels,cmap=cm.get_cmap('RdBu_r',len(levels)-1))
    plt.colorbar(cs)
    xticks(mytimes,myyears,rotation=-90)

    plotfile='figures/'+str(mytype)+'_'+str(myvar)+'_alldepths.pdf'
    plt.savefig(plotfile,dpi=300)
    print 'Saved figure file %s\n'%(plotfile)
    #plt.show()
    
# ------------------------------------------------------------------
#     MAIN
#     Trond Kristiansen, 09.07.2013, 11.02.2014
#     Trond.Kristiansen@imr.no
#     ------------------------------------------------------------------

# This script reads the output from running createAveragesOnHexagon_part2.py and plots the
# results as timeseries (Hoevmoeller diagram).

myvars=['temp','salt']

mydepths=[0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500,
           600, 700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]

myprefix='northsea_8km'

for myvar in myvars:
   
    first=True
    counter=0
    myyears=[]
    mytimes=[]
    """Calculate depth weighted and area integrated values as a function of time"""
    for depth in mydepths:
        
        print "Running: %s"%(depth)
       
        myfile= str(myprefix)+'_fldmean_'+str(myvar)+'_depth_'+str(depth)+'.nc'
           
        myCDF=Dataset(myfile)
        times=(myCDF.variables["time"][:])
        refDate=datetime.datetime(1948, 1, 1, 0, 0, 0)
        oldyear=-9
        if first is True:
            for t in times:
                current=refDate + datetime.timedelta(days=t)
                if current.month==6 and oldyear!=current.year:
                    myyears.append(current.year)
                    mytimes.append(t)
                    oldyear=current.year
                
                
        mydata=np.squeeze(myCDF.variables[myvar][:])
        ts = createTimeseries(mydata,times,refDate)
            
        if first is True:
            allData=np.zeros((len(mydepths),np.shape(mydata)[0]))
         
            first=False
        
        allData[counter,:]   = mydata     
        myCDF.close()
        counter+=1
    plotTimeseries(mytimes,myyears,times,allData,myvar,mydepths,"hoevmoller")
    
    # Calculate value relative to starting values
    allDataChange=np.zeros((np.shape(allData)))
    
    for i in xrange(len(mydepths)):
        y=0
        for t in times:
            current=refDate + datetime.timedelta(days=t)
            allDataChange[i,y]=allData[i,y]-allData[i,current.month-1]
            y+=1
          
    if myvar=="temp":
        plotTimeseries(mytimes,myyears,times,allDataChange,myvar,mydepths,"T-CLASS3-IC_CHANGE")
    if myvar=="salt":
        plotTimeseries(mytimes,myyears,times,allDataChange,myvar,mydepths,"S-CLASS3-IC_CHANGE")
    
    first=True