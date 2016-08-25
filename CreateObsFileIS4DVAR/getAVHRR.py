import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
import matplotlib.pyplot as plt
from pylab import *
import mpl_util


__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2015, 10, 19)
__modified__ = datetime.datetime(2015, 10, 19)
__version__  = "1.0"
__status__   = "Development, 19.10.2015"

def getLonLat(currentDate):

    if (currentDate.day >= 10):
        d="%s"%(currentDate.day)
    else:
        d="0%s"%(currentDate.day)
    if (currentDate.month >= 10):
        m="%s"%(currentDate.month)
    else:
        m="0%s"%(currentDate.month)

    base="/Users/trondkr/Projects/NOWMAPS/AVHRR/%s/"%(currentDate.year)
   # base="ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/NetCDF/%s/AVHRR/"%(currentDate.year)
    
    file="avhrr-only-v2.%s%s%s.nc"%(currentDate.year,m,d)
    
    filename=base+file
    print "test", filename

    cdf=Dataset(filename)
    longitude=np.squeeze(cdf.variables["lon"][:])
    latitude=np.squeeze(cdf.variables["lat"][:])
    sst=np.squeeze(cdf.variables["sst"][:])

    cdf.close();

    return longitude,latitude,sst

def openAVHRR(currentDate,indexes):
    if (currentDate.day >= 10):
        d="%s"%(currentDate.day)
    else:
        d="0%s"%(currentDate.day)
    if (currentDate.month >= 10):
        m="%s"%(currentDate.month)
    else:
        m="0%s"%(currentDate.month)

    base="/Users/trondkr/Projects/NOWMAPS/AVHRR/%s/"%(currentDate.year)
   # base="ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/NetCDF/%s/AVHRR/"%(currentDate.year)
    file="avhrr-only-v2.%s%s%s.nc"%(currentDate.year,m,d)
    print file
    filename=base+file
    cdf=Dataset(filename)

    print "Time: Extrating timedata from netCDF4 file: %s date: %s"%(filename,currentDate)
    
    time=np.squeeze(cdf.variables["time"][:])
    longitude=np.squeeze(cdf.variables["lon"][:])
 
    sst=cdf.variables["sst"][0,0,:,:]
    
    mydata,longitude = convertDataTo180(sst,longitude)
    longitude=longitude[int(indexes[0]):int(indexes[1])]
    sst=mydata[int(indexes[2]):int(indexes[3]),int(indexes[0]):int(indexes[1])]
    
    cdf.close()

    return time, sst, longitude


def convertDataTo180(mydata,longitude):
    i0 = np.argmin(np.fabs(longitude[:]-180))
    mydataout = np.ma.zeros(np.squeeze(mydata.shape,mydata.dtype))
    lonsout = np.ma.zeros((len(longitude[:])),longitude.dtype)
    lonsout[0:len(longitude[:])-i0] = longitude[i0:]-360
    lonsout[len(longitude[:])-i0:] = longitude[1:i0+1]
    mydataout[:,0:len(longitude[:])-i0]  = mydata[:,i0:]
    mydataout[:,len(longitude[:])-i0:] = mydata[:,1:i0+1]

    return mydataout, lonsout


def extractAVHRRLongLat(minLon,maxLon,minLat,maxLat,currentDate):
    print "extract avhhr long lat"
    """Routine that extracts the longitude and latitudes for the
    combination of tiles. This is only necessary to do once so it is separated
    from the extraction of SST."""
    longitude,latitude,sst = getLonLat(currentDate)

    mydata,longitude = convertDataTo180(sst,longitude)
    
    """ We have to flip this array so that we have increasing latitude
     values required by np.interp function. This means we also have to
     flip the input SST array"""
    #latitude=np.flipud(latitude)
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


    print "Extracted longitude-latitude for AVHRR region"
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




if __name__ == "__main__":

    main()
