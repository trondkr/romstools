import numpy as np
from netCDF4 import Dataset
import datetime
import os
import grid
import tools

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2014, 2, 20)
__modified__ = datetime.datetime(2014, 2, 20)
__version__  = "1.0"
__status__   = "Development, 20.2.2014"


doc ="""
     CALCULATEFLUXES.PY

     This program reads all of the ROMS output files (history or average files) located in a specified folder and
     calculates the volume fluxes through specified transects. The transects can be many and needs to be defined in the
     function defineTransects(). The start and end points are in grid coordinates (x, y). Some examples for usage are
     given in the function for typical transects in teh North Sea.

     REQUIREMENTS TO RUN:
     python modules: numpy, basemap, netCDF4, f2py (included in numpy)
     Fortran compiler

     Hexagon specific: Everything is installed just add
     "module load python" and "module load PrgEnv-gnu" in your .bashrc file.
     These modules provide access to the python modules and the Fortran compiler

     module swap PrgEnv-pgi PrgEnv-gnu
     module unload notur
     f2py --verbose  -c -m iso iso.f90
     f2py --verbose  -c -m obs_interp obs_interp.f

"""

def contourTransectData(mydepths,mylongs,data,sectionName,currentDate,mytype):
    if mytype=="U" or mytype=="V":
        levels=np.arange(-0.6,0.6,0.02)
    elif mytype=="S":
        levels=np.arange(34.8,35.5,0.05)
    elif mytype=="T":
        levels=np.arange(-2,12,0.5)
    print "Data range:",(data.min(), data.max())
    plt.figure(figsize=(10,6))
    cs1=contourf(mylongs,mydepths,data,levels,cmap=cm.get_cmap('RdYlBu_r',len(levels)-1),alpha=1.0,extend='both')
    cs2=contour(mylongs,mydepths,data,cs1.levels[::1],colors = ('k',),hold='on',linewidths = (0.2,))
    colorbar(cs1,extend='both')

    plotfile='figures/transect_'+str(sectionName)+'_'+str(currentDate)+'_var_'+str(mytype)+'.png'
    plt.savefig(plotfile,dpi=300)

    #plt.show()

def defineTransects():

    # NAME OF TRANSECT     (startI,startJ) - (endI, endJ) where I,J are x,y in grid coordinates
    # 1. Faroe Islands - Shetland (162,160) - (162,119)
    # 2. Shetland - Norway (Shelf) (167,118) - (167,83)
    # 3. Shetland - Norway (Trench) (167,83) - (167,62)
    # 4. Orkney - Shetland (131,115) - (155,115)
    # 5. Aberdeen - German Bight (124,101) - (124,7)
    # 6. Dover - Calais (60,30) - (60,25)
    # 7. Torungen - Hirtshals (177,44) - (177,29)
    # 8. Skagen - Gothenburg (181,27) - (190,27)

    startI = [162,167,167,131,124,60,177,181]
    startJ = [160,118,83,115,101,30,44,27]
    endI   = [162,167,167,155,124,60,177,190]
    endJ   = [119,83,62,115,7,25,29,27]
    transectNames = ["Faroe Islands - Shetland","Shetland - Norway (Shelf)", "Shetland - Norway (Trench)",
                     "Orkney - Shetland","Aberdeen - German Bight","Dover - Calais",
                     "Torungen - Hirtshals","Skagen - Gothenburg"]
    
    return startI, startJ, endI, endJ, transectNames

if __name__ == "__main__":

    # EDIT ----------------------------------------------------------------------------------
    # Add the necessary information for the path to your grid file and the files you want to
    # calculate fluxes through. These can be history or average files output from ROMS.
    # This script reads the content (unix command: ls) of a directory and stores all
    # files in a list. Therefore, the files of the directory should only be ROMS output files.
    #
    # Find all files to be used to calculate the volume fluxes
    # and sort them according to name (which here also sorts according to time)

    mypath="/work/users/trondk/NS8km/FINAL/"
    firstfile=mypath+'ns8km_avg_16438.nc'
    ncfiles=os.listdir(mypath)
    ncfiles.sort()
    grdfile='/work/users/trondk/NS8km/FORCING/GRID/nordsjoen_8km_grid_hmax20m_v4.nc'
    doHeat=False # flag to calculate heat fluxes
    doSalt=False # flag to calculate salt fluxes
    plotU=False; plotT=False; plotS=False; plotV=False

    # STOP EDIT -------------------------------------------------------------------------------


    # Create the grid object:
    grd = grid.get_ROMS_grid('NS8KM',hist_file=firstfile,grid_file=grdfile)
       
    # Define lists containing grid positions where transects should
    # be calculated:
    startI, startJ, endI, endJ, transectNames = defineTransects()

    # Loop over all transects and perform the volume flux calculations and
    # write result to file:
    
    
    for istart,jstart,iend,jend,transectName in zip(startI,startJ,endI,endJ,transectNames):
        
        
        outputfile="/work/users/trondk/NS8km/FLUX/results/volumefluxes_"+str(transectName)+".csv"
        if os.path.exists(outputfile): os.remove(outputfile)
        print "Opened outputfile: %s"%(outputfile)
         
        for ncfile in ncfiles:
            
            if (int(ncfile[-8:-3]) >= 16438):
               
                ncfile=mypath+ncfile
                
                out=open(outputfile,'a')
                nc = Dataset(ncfile)
             
                T = nc.variables['temp'][:]
                S = nc.variables['salt'][:]
                U = nc.variables['u'][:]
                V = nc.variables['v'][:]
               
                times=(nc.variables["ocean_time"][:])
                refDate=datetime.datetime(1948,1,1,0,0,0)
            
                minDepth=0
                maxDepth=5000
                currentDate=refDate + datetime.timedelta(seconds=times[0])
              
              
                for tstep in xrange(1):
                    currentDate=refDate + datetime.timedelta(seconds=times[tstep])
                    print "Date: %s"%(currentDate)
                    if currentDate.year >= 1993:
                        # Volume flux Net
                        transpu, transpv = tools.section_transport_z(U[tstep,:,:,:], V[tstep,:,:,:], grd,istart,iend,jstart,jend,h1=minDepth,h2=maxDepth)
                        
                        # Volume flux IN
                        posU=np.copy(U[tstep,:,:,:])
                        posU[ posU < 0 ] = 0
                        
                        posV=np.copy(V[tstep,:,:,:])
                        posV[ posV < 0 ] = 0
                       
                        postranspu, postranspv = tools.section_transport_z(posU, posV, grd,istart,iend,jstart,jend,h1=minDepth,h2=maxDepth)
                        
                        # Volume flux OUT
                        negU=np.copy(U[tstep,:,:,:])
                        negU[ negU > 0 ] = 0
                        
                        negV=np.copy(V[tstep,:,:,:])
                        negV[ negV > 0 ] = 0
                       
                        negtranspu, negtranspv = tools.section_transport_z(negU, negV, grd,istart,iend,jstart,jend,h1=minDepth,h2=maxDepth)
                        
                        # --------------------------
                        
                        if (doHeat):
                            # Heat flux
                            transpuT, transpvT = tools.section_tracer_transport_z(U[tstep,:,:,:], V[tstep,:,:,:], T[tstep,:,:,:],grd, istart,iend,jstart, jend,h1=minDepth,h2=maxDepth)
                
                        
                        # Get average values of tracers T and S through section. This function returns a tuple of data, mask, position
                        transectT = tools.transect(T[tstep,:,:,:], istart, iend, jstart, jend, grd, Cpos='rho', vert=False, spval=1e37)
                        mytransectT=np.ma.masked_where(abs(transectT[0])>50,transectT[0])
                        transectS = tools.transect(S[tstep,:,:,:], istart, iend, jstart, jend, grd, Cpos='rho', vert=False, spval=1e37)
                        
                        mytransectS=np.ma.masked_where(abs(transectS[0])>50,transectS[0])
                        
                        # PLOTTING of TRANSECTS--------------------------------
                        if plotU:    
                            transectU = tools.transect(U[tstep,:,:,:], istart, iend, jstart, jend, grd, Cpos='u', vert=False, spval=1e37)
                            mytransectU=np.ma.masked_where(abs(transectU[0])>50,transectU[0])
                        if plotV:    
                            transectV = tools.transect(V[tstep,:,:,:], istart, iend, jstart, jend, grd, Cpos='v', vert=False, spval=1e37)
                            mytransectV=np.ma.masked_where(abs(transectV[0])>50,transectV[0])
                
                
                        if plotT is True:
                            mydepths=transectT[1]
                            mylongs=transectT[2]
                            mylats=transectT[3]
                            contourTransectData(mydepths,mylongs,mytransectT,sectionName,currentDate,"T")
                        if plotU is True:
                            mydepths=transectU[1]
                            mylongs=transectU[2]
                            mylats=transectU[3]
                            contourTransectData(mydepths,mylongs,-mytransectU,sectionName,currentDate,"U")
                        if plotS is True:
                            mydepths=transectS[1]
                            mylongs=transectS[2]
                            mylats=transectS[3]
                            contourTransectData(mydepths,mylongs,mytransectS,sectionName,currentDate,"S")
                
                        # ------------------------------------------------------- #
                        
                        print "Mean temperature      : %5.2f degC" % (np.mean(mytransectT))
                        print "Mean salinity         : %5.2f psu"  % (np.mean(mytransectS))
                
                        flux=np.sqrt(transpu*transpu + transpv*transpv)
                        print "Net volume flux           : %5.2f Sv"   %(float(flux)/1.e6)
                        print "Net northwards flux       : %5.2f Sv"   %(float(transpv)/1.e6)
                        print "Net southwards flux       : %5.2f Sv"   %(float(transpu)/1.e6) # Have to flip the sign!!!!!
                        
                        print "Positive northwards flux       : %5.2f Sv"   %(float(postranspv)/1.e6)
                        print "Positive southwards flux       : %5.2f Sv"   %(float(postranspu)/1.e6) # Have to flip the sign!!!!!
                
                        print "Negative northwards flux       : %5.2f Sv"   %(float(negtranspv)/1.e6)
                        print "Negative southwards flux       : %5.2f Sv"   %(float(negtranspu)/1.e6) # Have to flip the sign!!!!!
                
                
                        rho = 1025.0
                        #print "Salt flux             : %5.2f Sv" % ((rho * flux * 0.001*np.mean(mytransectS))/1.e6)
            
                        heat_cap = 4000.0 # J/(kg K)
                        
                            # units
                            # fluxT = C m3/s
                            # rho  = kg/m3
                            # heat_cap = J /Kg C
                            # Heat flux = (kg/m3)*(J/KgK)*(Km3/s) = J/s = Watt = TW /1.e12
                        if (doHeat):
                            print "Heat flux             : %5.2f TW"   % ((rho * heat_cap * fluxT)/1.0E12)
                            
                        outdata="%s, %s, %10.2f, %10.2f, %10.2f, %10.2f, %10.2f, %10.2f, %10.2f, %10.2f\n"%(times[tstep],
                                                            currentDate,
                                                            np.mean(mytransectT),
                                                            np.mean(mytransectS),
                                                            float(postranspu)/1.e6,
                                                            float(postranspv)/1.e6,
                                                            float(negtranspu)/1.e6,
                                                            float(negtranspv)/1.e6,
                                                            float(transpu)/1.e6,
                                                            float(transpv)/1.e6)
                        
                        out.writelines(outdata)
                    out.close()
                    nc.close()
                
main()
