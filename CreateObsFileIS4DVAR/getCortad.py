import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
import matplotlib.pyplot as plt
from pylab import *
import mpl_util
import time

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 12, 3)
__modified__ = datetime.datetime(2012, 12, 3)
__version__  = "1.0"
__status__   = "Development, 03.12.2012"

def getCORTADtime():

    base="http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version4/"
    file1="cortadv4_row01_col07.nc"
    filename1=base+file1
    cdf1=Dataset(filename1)

    print "Time: Extrating timedata from openDAP: %s"%(filename1)

    time=np.squeeze(cdf1.variables["time"][:])
    cdf1.close()
    return time

def openCoRTAD(maxTries,delay):
    """ Info on the different tiles used to identoify a region is found here:
    http://www.nodc.noaa.gov/SatelliteData/Cortad/TileMap.jpg"""

    base="http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version4/"
    hostname="data.nodc.noaa.gov"
    
   # response = os.system("ping -c 1 " + hostname)
    #while maxTries >0:
    #    maxTries -= 1
       # if maxTries<0:
       #     print "Check for network connection with openDAP server. Tries %s times without successful connection!"
       #     sys.exit()
                
       # if response == 0:
       #     print "NOAA is up: %s"%(response)
    
    file1="cortadv4_row01_col07.nc"
    file2="cortadv4_row01_col08.nc"
    file3="cortadv4_row00_col07.nc"
    file4="cortadv4_row00_col08.nc"

    filename1=base+file1
    filename2=base+file2
    filename3=base+file3
    filename4=base+file4

    cdf1=Dataset(filename1)
    cdf2=Dataset(filename2)
    cdf3=Dataset(filename3)
    cdf4=Dataset(filename4)

    return cdf1,cdf2,cdf3,cdf4
       # else:
       #     time.sleep (delay)
       #     continue
                

def extractCoRTADLongLat(maxTries,delay):
    """Routine that extracts the longitude and latitudes for the
    combination of tiles. This is only necessary to do once so it is separated
    from the extraction of SST."""
    cdf1,cdf2,cdf3,cdf4=openCoRTAD(maxTries,delay)

    longitude1=np.squeeze(cdf1.variables["lon"][:])
    latitude1=np.squeeze(cdf1.variables["lat"][:])
    
    longitude2=np.squeeze(cdf2.variables["lon"][:])
    latitude2=np.squeeze(cdf2.variables["lat"][:])

    longitude3=np.squeeze(cdf3.variables["lon"][:])
    latitude3=np.squeeze(cdf3.variables["lat"][:])

    longitude4=np.squeeze(cdf4.variables["lon"][:])
    latitude4=np.squeeze(cdf4.variables["lat"][:])

    cdf1.close();cdf2.close();cdf3.close();cdf4.close();

    longitude=concatenate((longitude1,longitude2),axis=1)
    latitude=concatenate((latitude3,latitude1),axis=0)

    """ We have to flip this array so that we have increasing latitude
     values required by np.interp function. This means we also have to
     flip the input SST array"""
    latitude=np.flipud(latitude)
    lons,lats=np.meshgrid(longitude,latitude)

    print "Extracted longitude-latitude for CoRTAD region"
    print "Long min: %s Long max: %s"%(longitude.min(),longitude.max())
    print "Lat min: %s Lat max: %s"%(latitude.min(),latitude.max())
    print "------------------------------\n"
    return lons, lats, longitude, latitude

def extractCORTADSST(name,t,cdf1,cdf2,cdf3,cdf4):
    """Routine that extracts the SST values for
    the specific tiles and time-period (t)"""
    
    filledSST1=cdf1.variables["FilledSST"][t,:,:]
    filledSST2=cdf2.variables["FilledSST"][t,:,:]
    filledSST3=cdf3.variables["FilledSST"][t,:,:]
    filledSST4=cdf4.variables["FilledSST"][t,:,:]
    offset=cdf1.variables["FilledSST"].__getattribute__('add_offset')
    
    filledMaskedSST1=filledSST1 - offset
    filledMaskedSST2=filledSST2 - offset
    filledMaskedSST3=filledSST3 - offset
    filledMaskedSST4=filledSST4 - offset

    """Now we have all the data in 4 different arrays that we need to concentate.
    First we add the horisontal tiles, and finally we stack the two horisontal ones on top
    of each other."""
    filledMaskedSST_lower=concatenate((filledMaskedSST1,filledMaskedSST2),axis=1)
    filledMaskedSST_upper=concatenate((filledMaskedSST3,filledMaskedSST4),axis=1)
    filledMaskedSST_all=concatenate((filledMaskedSST_upper,filledMaskedSST_lower),axis=0)

    """Flip the SST array to be consistent with order of latitude array"""
    filledMaskedSST_all=np.flipud(filledMaskedSST_all)

    """ Scale and offset is autmoatically detected and edited by netcdf, but
    we need to mask the values that are not filled."""
    filledMaskedSST_final=ma.masked_less(filledMaskedSST_all,-2.)

    print "Min and max of SST: %s - %s"%(filledMaskedSST_final.min(),filledMaskedSST_final.max())
    print "------------------------------\n"

    return filledMaskedSST_final


        


if __name__ == "__main__":

    main()
