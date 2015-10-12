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
__modified__ = datetime.datetime(2014, 11, 10)
__version__  = "1.0"
__status__   = "Development, 03.12.2012, 10.11.2014"

def getCORTADtime():

    base="http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version5/"

    file1="cortadv5_FilledSST.nc"


    filename1=base+file1
    cdf1=Dataset(filename1)

    print "Time: Extrating timedata from openDAP: %s"%(filename1)
    
    time=np.squeeze(cdf1.variables["time"][:])
    
    cdf1.close()
    return time

def openCoRTAD(maxTries,delay):
    """ Info on the different tiles used to identify a region is found here:
    http://www.nodc.noaa.gov/SatelliteData/Cortad/TileMap.jpg"""

    base="http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version5/"
    #base="http://data.nodc.noaa.gov/cortad/Version5/"

    file1="cortadv5_FilledSST.nc"

    hostname="data.nodc.noaa.gov"
    
   # response = os.system("ping -c 1 " + hostname)
    #while maxTries >0:
    #    maxTries -= 1
       # if maxTries<0:
       #     print "Check for network connection with openDAP server. Tries %s times without successful connection!"
       #     sys.exit()
                
       # if response == 0:
       #     print "NOAA is up: %s"%(response)

    filename1=base+file1
    cdf1=Dataset(filename1)

    return cdf1
       # else:
       #     time.sleep (delay)
       #     continue
                

def extractCoRTADLongLat(maxTries,delay,minLon,maxLon,minLat,maxLat):
    """Routine that extracts the longitude and latitudes for the
    combination of tiles. This is only necessary to do once so it is separated
    from the extraction of SST."""
    cdf1 = openCoRTAD(maxTries,delay)

    longitude=np.squeeze(cdf1.variables["lon"][:])
    latitude=np.squeeze(cdf1.variables["lat"][:])

    cdf1.close();

    """ We have to flip this array so that we have increasing latitude
     values required by np.interp function. This means we also have to
     flip the input SST array"""
   # latitude=np.flipud(latitude)
    lons,lats=np.meshgrid(longitude,latitude)

    print "Full data range:"
    print "Long min: %s Long max: %s"%(longitude.min(),longitude.max())
    print "Lat min: %s Lat max: %s"%(latitude.min(),latitude.max())

    print "Find indexes for SST that define the grid domain"

    res = findSubsetIndices(minLat,maxLat,minLon,maxLon,latitude,longitude)

    # Note that at this point the 2D arrays are flipped in order so [lat,lon]
    print "Wanted: %3.3f %3.3f"%(minLon,minLat)
    print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[res[0]],latitude[res[2]])
    print "2D test: %3.3f %3.3f\n"%(lons[res[2],res[0]],lats[res[2],res[0]])

    print "Wanted: %3.3f %3.3f"%(minLon,maxLat)
    print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[res[0]],latitude[res[3]])
    print "2D test: %3.3f %3.3f\n"%(lons[res[3],res[0]],lats[res[3],res[0]])

    print "Wanted: %3.3f %3.3f"%(maxLon,minLat)
    print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[res[1]],latitude[res[2]])
    print "2D test: %3.3f %3.3f\n"%(lons[res[2],res[1]],lats[res[2],res[1]])

    print "Wanted: %3.3f %3.3f"%(maxLon,maxLat)
    print 'Index corresponds to longitude [%3.3f] and latitude [%3.3f]'%(longitude[res[1]],latitude[res[3]])
    print "2D test: %3.3f %3.3f\n"%(lons[res[3],res[1]],lats[res[3],res[1]])


    print "Extracted longitude-latitude for CoRTAD region"
    print "------------------------------\n"

    return longitude, latitude, lons, lats, res


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
    minJ=int(indices[1][2])
    maxJ=int(indices[0][2])
    minI=int(indices[3][2])
    maxI=int(indices[2][2])
    res[0]=minI; res[1]=maxI; res[2]=minJ; res[3]=maxJ;
    return res



def extractCORTADSST(name,t,cdf,indexes):
    # Routine that extracts the SST values for the specific tiles and time-period (t)

    filledSST=cdf.variables["FilledSST"][t,indexes[3]:indexes[2],indexes[0]:indexes[1]]
    offset=cdf.variables["FilledSST"].__getattribute__('add_offset')
    # Kelvin to Celsius
    filledSST = filledSST-offset
    #filledMaskedSST=np.where(filledSST > 270, filledSST  - offset, filledSST)

    """ Scale and offset is autmoatically detected and edited by netcdf, but
    we need to mask the values that are not filled."""
    #filledMaskedSST_final=ma.masked_less(filledMaskedSST,-2.)
    #filledMaskedSST_final=ma.masked_greater(filledMaskedSST,25.)

    print "Min and max of SST: %s - %s"%(filledSST.min(),filledSST.max())
    print "------------------------------\n"

    return filledSST


        


if __name__ == "__main__":

    main()
