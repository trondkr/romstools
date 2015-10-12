
import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from pylab import *
import mpl_util
import matplotlib as mpl
import mpl_toolkits.basemap as mp
import getCortad
import writeObsfile
import time

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 4, 8)
__modified__ = datetime.datetime(2013, 11, 17)
__version__  = "1.0"
__status__   = "Development, 08.04.2012, 30.11.2012, 03.12.2012, 07.06.2013, 17.11.2013"

def help():

    """
    This is a program for creating an observation file used by IS4DVAR
    for assimilation simulations of ROMS. The program extracts SST data
    from the CoRTAD database using OpenDAP. The North Sea domain consists
    of 4 different CoRTAD tiles that need to be extracted separately, then assembled into
    one big array which is further modified by the boundaries of the ROMS
    grid file. Finally, the values within the grid domain are written to
    netCDF4 file that keeps track of all of the observations. Each observation
    has a time-stamp, a longitude-latitude value and an observation value (SST value).

    This program can also plot the grid domain boundaries using a polygon
    created by reading the grid file, plot the bathymtetry from Etopo1, and plot
    the extracted SST values for the domain. In addition, the grid points in the
    grid file can be plotted at strides to visualize the layout of grid cells
    relative to the SST values.

    Trond Kristiansen, 08.04.2012, 30.11.2012, 03.12.2012
    """



def addEtopo1Bathymetry(ax,map):
    """Get the etopo1 data"""
    etopo1name='/Users/trond/Projects/arcwarm/maps/ETOPO1_Ice_g_gmt4.grd'
    etopo1 = Dataset(etopo1name,'r')

    lons = etopo1.variables["lon"][:]
    lats = etopo1.variables["lat"][:]

    res = findSubsetIndices(latStart-5,latEnd+5,lonStart-40,lonEnd+10,lats,lons)

    lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]])
    print "Extracted data for area %s : (%s,%s) to (%s,%s)"%(name,lon.min(),lat.min(),lon.max(),lat.max())
    bathy = etopo1.variables["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])]
    bathySmoothed = laplaceFilter.laplace_filter(bathy,M=None)

    levels=[-6000,-5000,-3000, -2000, -1500, -1000,-500, -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, -25, -15, -10, -5, 0]

    x, y = map(lon,lat)
    CS1 = map.contourf(x,y,bathySmoothed,levels,
                           cmap=mpl_util.LevelColormap(levels,cmap=cm.Blues_r),
                           #cmap=cm.get_cmap('Greys_r',len(levels)),
                           extend='upper',
                           alpha=1.0,
                           origin='lower',
                           rasterized=True)

    CS1.axis='tight'

def addGridToMap(map,grid_lon,grid_lat,doEveryTenthX):
    doEveryTenthY=doEveryTenthX
    doEveryTenthCounterX=0; doEveryTenthCounterY=0
    print "Starting plotting station grid on map"
    for j in range(len(grid_lon[0,:])):
        if doEveryTenthCounterY==doEveryTenthY:
            for jj in range(len(grid_lon[:,0])):
                if doEveryTenthCounterX==doEveryTenthX:
                    xpt,ypt = map(grid_lon[jj,j],grid_lat[jj,j])
                #    map.plot([xpt],[ypt],marker="o",color="white", markersize=10)
                    map.plot([xpt],[ypt],marker="o",color="red", markersize=3.8)
                    doEveryTenthCounterX=0
                    doEveryTenthCounterY=0
                doEveryTenthCounterX+=1

        doEveryTenthCounterY+=1


def drawSST(ax,map,longSST,latSST,filledSST,myalpha):
    # Input arrays has to be 2D
    print "Drawing SST: min %s and max %s"%(filledSST.min(), filledSST.max())
    x2, y2 = map(longSST,latSST)
    levels=np.arange(-2,18,0.5)
    #levels=np.arange(filledSST.min(), filledSST.max(),0.2)

    if myalpha > 0.99:
        #CS2 = map.contourf(x2,y2,filledSST,levels,cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r),extend='both',alpha=myalpha)
        CS2 = map.pcolor(x2,y2,filledSST)
     
        plt.colorbar(CS2,orientation='vertical',extend='both', shrink=0.5)
    else:
        CS2 = map.contourf(x2,y2,filledSST,levels,cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r),extend='both',alpha=myalpha)
        #CS2 = map.contourf(x2,y2,filledSST,levels,cmap=mpl_util.LevelColormap(levels,cmap=cm.Greys),extend='upper',alpha=myalpha)



def makeMap(figureNumber,lonStart,lonEnd,latStart,latEnd,name,SSTi,lon_rho,lat_rho,polygon_data,currentDate,origSST,lonSST2D,latSST2D):
    """ Standard map settings-----------------------------------"""
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_subplot(111)
    doBathymetry=False
    if lonStart< 0 and lonEnd < 0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    elif lonStart > 0 and lonEnd > 0:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=((lonEnd)+(lonStart))/2.0

    print 'Creating map with center longitude ',lon_0

    map = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,\
            llcrnrlon=lonStart,urcrnrlon=lonEnd,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=200.,projection='lcc',\
            lat_1=latStart,lon_0=lon_0)

    map.drawcoastlines()
    map.drawcountries()
 #   map.fillcontinents(color='grey')
    map.drawmeridians(np.arange(lonStart,lonEnd,10),labels=[0,0,0,1])
    map.drawparallels(np.arange(latStart,latEnd,4),labels=[1,0,0,0])

    """ Add specifics to map--------------------------------------"""
    if doBathymetry is True:
        addEtopo1Bathymetry

    """Plot grid positions as dots"""
    #addGridToMap(map,grid_lon,grid_lat,doEveryTenthX=10)

    """Draw original SST"""
    drawSST(ax,map,lonSST2D,latSST2D,origSST,0.5)

    """Draw interpolated  SST"""
    drawSST(ax,map,lon_rho,lat_rho,SSTi,1.0)

    """Draw grid boundaries"""
    drawGridBoundaries(ax,map,polygon_data)

    """Save the map to file and finish ---------------------------"""
    plt.title('Grid: %s SST time:%s'%(name,currentDate))
    if not os.path.exists("figures"): os.mkdir("figures")
    plotfile='figures/SST_northsea_'+str(figureNumber)+'.png'
    print "Saving map to file %s"%(plotfile)
    plt.savefig(plotfile)
    #plt.show()


def drawGridBoundaries(ax,map,polygon_data):
    print "Plotting grid boundaries"
    polygon_data_xy=[]
    for i in xrange(len(polygon_data[0,:])):
        myx,myy=map(polygon_data[0,i],polygon_data[1,i])
        vertices=[myx,myy]
        polygon_data_xy.append(vertices)

    polygon_data_xy=np.asarray(polygon_data_xy)
    patch = Polygon(array(polygon_data_xy), facecolor='none',
        edgecolor=(.9, .3, .3, 1), linewidth=4)
    ax.add_patch(patch)

def getGrid(filename):
    cdf=Dataset(filename)
    grid_lon=cdf.variables["lon_rho"][:]
    grid_lat=cdf.variables["lat_rho"][:]
    mask_rho=cdf.variables["mask_rho"][:]
    grid_h=cdf.variables["h"]
    print "Grid dimmensions: %s and %s"%(grid_lon.shape,grid_lat.shape)
    cdf.close()

    grid_lon_min=grid_lon.min()
    grid_lon_max=grid_lon.max()

    grid_lat_min=grid_lat.min()
    grid_lat_max=grid_lat.max()
    print "Grid domain longitude from %s to %s"%(grid_lon_min,grid_lon_max)
    print "Grid domain latitude from %s to %s"%(grid_lat_min,grid_lat_max)
    print "----\n"
    return mask_rho, grid_lon, grid_lat,grid_h


def ingrid(lon, lat, lon_bd,lat_bd):
    return mpl.mlab.inside_poly(zip(lon, lat),  zip(lon_bd, lat_bd))


def getPolygon(lonSST,latSST,grid_lon,grid_lat):

    lon_bd = np.concatenate((grid_lon[:,0],grid_lon[-1,:],grid_lon[::-1,-1], grid_lon[0,::-1] ))
    lat_bd = np.concatenate((grid_lat[:,0],grid_lat[-1,:],grid_lat[::-1,-1], grid_lat[0,::-1] ))

    """ Save the polygon as array to plot later"""
    polygon_data=np.empty((2,len(lon_bd)))
    for k in xrange(len(lon_bd)):
        polygon_data[0,k]=lon_bd[k]
        polygon_data[1,k]=lat_bd[k]

    return polygon_data


def main():
    startTime = time.time()
    """Define the start and end date you want data extracted for:"""
    startYear=2009
    startMonth=10
    endYear=2012
    endMonth=12
    maxTries=3
    delay=10
    firstIteration=True
    lastIteration=False
    createFigure=False
    figureNumber=0
    USENETCDF4=True    # if false then use NETCDF3_CLASSIC


    """Name of output file to be created"""
    outputFile="NS8KM_obsSST_%s_to_%s.nc"%(startYear,endYear)
    if os.path.exists(outputFile): os.remove(outputFile)

    """Read the grid info from the grid file"""
    filename="/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m_v3.nc"
    mask_rho, lon_rho,lat_rho,grid_h = getGrid(filename)

    """Calculate the x,y grid coordinates"""
    (Mp,Lp)=lon_rho.shape
    X=np.arange(0,Mp,1)
    Y=np.arange(0,Lp,1)

    roms_Xgrid,roms_Ygrid=np.meshgrid(Y,X)

    """CoRTAD time is days since 1980/12/31 12:00:00"""
    mytime=getCortad.getCORTADtime()
    refDate=datetime.datetime(1981,12,31,12,0,0)

    """Have to convert the day of observation to the relative time used by ROMS
    which is 1948/1/1:00:00:00"""
    refDateROMS=datetime.datetime(1948,1,1,0,0,0)
    delta=refDate-refDateROMS
    daysSince1948to1980=delta.days

    """Find the start and end indexes to extract"""
    foundStart=False; foundEnd=False; startIndex=-9; endIndex=-9
    for index in xrange(len(mytime)):
        currentDate = refDateROMS + datetime.timedelta(days=float(mytime[index])+daysSince1948to1980)
        if foundStart is False:
            if currentDate.year==startYear:
                if currentDate.month==startMonth:
                    foundStart=True
                    startIndex=index
                    print "\n-----------------------------------------------"
                    print "Start date %s at index %s"%(currentDate,startIndex)

        if foundEnd is False:
            if currentDate.year==endYear:
                if currentDate.month==endMonth:
                    foundEnd=True
                    endIndex=index
                    print "FIXME : HARDCODING LAST INDEX !!!!!!!!!!!!!!!!!\n\n\n"
                    endIndex=1616
                    currentDate = refDateROMS + datetime.timedelta(days=float(mytime[endIndex])+daysSince1948to1980)
                    print "FIXME : HARDCODING LAST INDEX !!!!!!!!!!!!!!!!!\n\n\n"

                    print "End date %s at index %s"%(currentDate,endIndex)
    times=[i for i in range(startIndex,endIndex,1)]
    print "Created array of %s time-steps to iterate and extract data from"%(len(times))
    print "-----------------------------------------------\n"

    """Get the lomgitude-latitudes of the combination of tiles"""
    longitude, latitude, lonSST, latSST, indexes = getCortad.extractCoRTADLongLat(maxTries,
                                                                        delay,
                                                                        lon_rho.min(),
                                                                        lon_rho.max(),
                                                                        lat_rho.min(),
                                                                        lat_rho.max())
    indexes=np.asarray(indexes,dtype=np.int32)
    latitude  = np.flipud(latitude[indexes[3]:indexes[2]])
    longitude = longitude[indexes[0]:indexes[1]]

    """Loop over all times and store to file or make map"""
    polygon_data = getPolygon(lonSST[indexes[3]:indexes[2],indexes[0]:indexes[1]],
                              latSST[indexes[3]:indexes[2],indexes[0]:indexes[1]],
                              lon_rho,lat_rho)
    survey_time=[]

    for t in xrange(len(times)):
        
        """Open the files and check that NOAA is online"""
        cdf = getCortad.openCoRTAD(maxTries,delay)
        currentDate=refDateROMS + datetime.timedelta(days=int(mytime[times[t]])+daysSince1948to1980)

        
        """ Get the data for the current time"""
        filledSST = getCortad.extractCORTADSST("North Sea",times[t],cdf,indexes)

        """Interpolate the original values to the grid. This is the data that will be saved to file"""

        SSTi = mp.interp(np.flipud(filledSST),longitude,latitude,
                             lon_rho,lat_rho,checkbounds=False,masked=True,order=1)

        SSTi = np.where(SSTi < -0.5, -0.5, SSTi)

        SSTi = SSTi*mask_rho

        igood=np.nonzero(SSTi)
        numberOfobs=len(SSTi[igood])

        obs_lon=lon_rho[igood]
        obs_lat=lat_rho[igood]
        obs_value=SSTi[igood]
        obs_Xgrid=roms_Xgrid[igood]
        obs_Ygrid=roms_Ygrid[igood]
        Nobs=numberOfobs
        survey_time.append(int(mytime[times[t]])+daysSince1948to1980)
       
        obs_time=[]
        for ot in xrange(numberOfobs):
            obs_time.append(int(mytime[times[t]])+daysSince1948to1980)
            if ot==0:
                print refDateROMS + datetime.timedelta(days=int(mytime[times[t]])+daysSince1948to1980), int(mytime[times[t]])+daysSince1948to1980

        print "Found %s observations for %s"%(numberOfobs, currentDate)

        """Create map where the colored data shows the interpolated values and the
            grey colored data are the original data"""
        """Define the max and minimim area to crate map for (not used to create obs file)"""
        lat_start=43; lat_end=71.5; lon_start=-20; lon_end=35

        if createFigure is True:
            makeMap(figureNumber,lon_start,lon_end,lat_start,lat_end,filename,SSTi,lon_rho,lat_rho,polygon_data,currentDate,
                        filledSST,lonSST[indexes[3]:indexes[2],indexes[0]:indexes[1]],
                                  latSST[indexes[3]:indexes[2],indexes[0]:indexes[1]])
            figureNumber+=1

        """ Finished, now cleanup and make sure everything are arrays"""
        obs_time=np.asarray(obs_time)


        """Finally write the results to file"""

        """Temp variables not used until lastIteration is set to True, but required for function call"""
        obs_flag = 6; is3d = 1; survey =0; Nstate = 7
        if firstIteration is True:
            print "Writing data of TYPE: %s to file (6=Temperature)"%(obs_flag)

        unos = np.ones(len(obs_value))
        obs_type = obs_flag*unos
        obs_error = unos   # error eqaul one scale later
        obs_Zgrid = 0*unos
        obs_depth = 35*unos #If positive has to be the sigma level, if negative depth in meters
        obs_variance=np.asarray(np.ones(Nstate))


        print "Min and max of SST to file: %s - %s"%(obs_value.min(),obs_value.max())


        writeObsfile.writeData(outputFile,obs_lat,obs_lon,obs_value,Nobs,survey_time,obs_time,obs_Xgrid,obs_Ygrid,
                                   firstIteration,lastIteration,
                                   obs_flag,obs_type,obs_error,obs_Zgrid,obs_depth,obs_variance,
                                   survey,is3d,Nstate,USENETCDF4)
        firstIteration=False
        """Close the opendap files"""
        cdf.close();
    

    """Cleanup and write final dimensions and variables"""
    lastIteration=True
    """ some extra variables """

    obs_flag = 6       # for temperature data
    is3d = 1
    survey=len(survey_time)
    survey_time=np.asarray(survey_time)
    survey_time=survey_time.flatten()
    Nstate = 7;

    writeObsfile.writeData(outputFile,obs_lat,obs_lon,obs_value,Nobs,survey_time,obs_time,obs_Xgrid,obs_Ygrid,
                               firstIteration,lastIteration,
                               obs_flag,obs_type,obs_error,obs_Zgrid,obs_depth,obs_variance,
                               survey,is3d,Nstate,USENETCDF4)

    endTime=time.time()
    print "\n--------------------------------------------------------------"
    print "Program ended successfully after %s seconds"%(endTime-startTime)
    print "\n--------------------------------------------------------------\n"
if __name__ == "__main__":

    main()
