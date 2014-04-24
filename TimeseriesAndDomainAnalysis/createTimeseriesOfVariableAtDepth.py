from pylab import *
import os, sys, string
import paramiko
import os
import pandas as pd
import datetime as datetime
import pandas as pd
from netCDF4 import Dataset

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2014, 1, 22)
__modified__ = datetime.datetime(2014, 1, 22)
__version__  = "1.0"
__status__   = "Production"


# This script reads the output from running createAveragesOnHexagon_part2.py and plots the
# result as time-series. Specify what variable you want and what depth level below.

   
def createTimeseries(myCDF,myvar,times,refDate):
  
    if myvar in ["temp","salt"]:
        mydata=np.squeeze(myCDF.variables[myvar][:,0,:,:])
    else:
        mydata=np.squeeze(myCDF.variables[myvar][:,:,:])
    myvalues=[]
    mydates=[]
    for x in xrange(len(times[:])):
        time=times[x]
        myvalue=mydata[x]
        if myvar=="zeta":
            """meter to centimeter"""
            myvalue=myvalue*100.
        if myvar=="ssflux":
            """m s-1 to mm year-1"""
            myvalue=myvalue #*1000*86400*365.
            
        currentDate=refDate + datetime.timedelta(days=time)
      #  if currentDate.year==1996 or (currentDate.year==1995 and currentDate.month==12):
      #      if myvar=="temp":
      #          myvalue+=3
          
       # if currentDate.year==1994:  
        mydates.append(currentDate)
        myvalues.append(myvalue)
        
    ts=pd.Series(myvalues,mydates)
    
    print ts
    return ts
         
def plotTimeseries(ts,myvar):
    
    ts_annual = ts.resample("A")
    ts_quarterly = ts.resample("Q")
    ts_monthly = ts.resample("M")
      
  #  ts_quarterly.plot(style="orange",linewidth=5, label="Seasonal")
  #  ts_annual.plot(style="r",linewidth=6, label="Annual")
  #  ts.plot(style="r", marker='o', linewidth=1,label="Monthly raw",alpha=0.5)
    
    ts_monthly.plot(style="b", marker='o', linewidth=2,label="Monthly")
    
    #legend(loc='best')
    if myvar=="temp":
        ylabel('Temperature ($^\circ$C)')
    if myvar=="salt":
        ylabel('Salinity (psu)')
    if myvar=="shflux":
        ylabel('Watt m-2')
    if myvar=="zeta":
        ylabel('SLA (cm)')
      
    plotfile='figures/timeseries_'+str(myvar)+'.pdf'
    plt.savefig(plotfile,dpi=300)
    print 'Saved figure file %s\n'%(plotfile)
    plt.show()
    
def main():
    
    
    # ------------------------------------------------------------------
    #     MAIN
    #     Trond Kristiansen, 09.07.2013, 11.02.2014
    #     Trond.Kristiansen@imr.no
    #     ------------------------------------------------------------------
    #
    #
    # This script reads the output from running createAveragesOnHexagon_part2.py and plots the
    # result as time-series. Specify what variable you want and what depth level below.

    myvars=['ssflux','zeta']
    myvars=['u_eastward','v_northward']
   # myvars=['ssflux','shflux']
    myvars=['temp']
    myprefix='northsea_8km'
    mylevel='0'
    
   
    for myvar in myvars:
          
        if myvar in ["temp","salt",'u_eastward','v_northward']:
            myfile= str(myprefix)+'_fldmean_'+str(myvar)+'_depth_'+str(mylevel)+'.nc'
        elif myvar in ["shflux","ssflux","zeta"]:
            myfile= str(myprefix)+'_fldmean_'+str(myvar)+'.nc'
        print myfile
                
        myCDF=Dataset(myfile)
        times=(myCDF.variables["time"][:])
        refDate=datetime.datetime(1948, 1, 1, 0, 0, 0)
       
        ts = createTimeseries(myCDF,myvar,times,refDate)
        
        plotTimeseries(ts,myvar)   
        myCDF.close()
      
if __name__=="__main__":
    main()
    
    
    