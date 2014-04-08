import os
import datetime as datetime
from netCDF4 import Dataset
from subprocess import call
from itertools import chain

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2014, 4, 2)
__modified__ = datetime.datetime(2014, 4, 2)
__version__  = "1.0.0"
__status__   = "Production"


# DELETETIMESTEPSINNETCDFILES.PY
#
# This script reads netcdf files and investigates whether some timesteps occurr more than one time. This
# can happen when you concatenate many netcdf files using ncrcat (NCO - http://nco.sourceforge.net/nco.html)
# This script is not fast but it works.
#
# The example input file has the following format: nordsjoen_8km_Pair.nc
#
# Required: nco (specifically the ncrcat routine)
# Hexagon: use "load module nco"
#
# ------------------------------------

myvars=['Pair','Uwind','Vwind','cloud','Tair','Qair','rain']

myprefix='nordsjoen_8km_'

for myvar in myvars:
    filename="%s%s_2000_2009.nc"%(myprefix,myvar)
    newfilename="%s%s_2000_2009_edited.nc"%(myprefix,myvar)

    print "Working on file: %s"%(filename)

    cdf=Dataset(filename)
    timevar="%s_time"%(myvar)
    if myvar=="Pair":  timevar="pair_time"
    if myvar=="cloud":  timevar="cloud_time"
    if myvar in ["Uwind","Vwind"]:  timevar="wind_time"

    mytimes=cdf.variables[timevar][:]

    times={}
    mykeep=[]
    mykeepdays=[]
    mydelete=[]
    counter=0
    mylistofTMPfiles=[]

    for mytime in mytimes:

        if mytime not in times:
            times[mytime]=counter
            mykeep.append(counter)
            mykeepdays.append(mytime)
        else:
            mydelete.append(mytime)
        counter+=1
    refdate=datetime.datetime(1948,1,1,0,0)

    for deletetime in mydelete:

        print "Delete %s -> %s"%(deletetime, refdate + datetime.timedelta(days=float(deletetime)))
    print "Keeping %s timesteps in new netcdf files"%(len(times.keys()))

    mylist=[int(i) for i in mykeep]
    it=chain(mylist)
    opt=""

    while True:
         try:
             elem = it.next()
             opt=opt+str(elem)
         except StopIteration:
             print "Last element was:", elem, ". Finished creating list"
             break
         opt=opt+","
    opt=opt[0:-1]

    if len(mylist)>0:
        options="-d %s,%s %s %s"%(timevar,opt,filename,newfilename)
        command = "ncks %s"%(options)

        for timestep in mykeep:
            outfile="subset/tmp_%s"%(timestep)
            mylistofTMPfiles.append(outfile)

            command="ncks -d %s,%s %s %s"%(timevar,timestep,filename,outfile)
            print "running command: %s"%(command)
            call(command,shell=True)

        k=0
        for outfile in mylistofTMPfiles:
            if k==0:
                myfiles="%s"%(outfile)
            else:
                myfiles=myfiles + " %s"%(outfile)
            k+=1

        command="ncrcat %s %s"%(myfiles,newfilename)
        print command
        call(command,shell=True)

        # Delete all the subset files after successful ncrcat
        command="rm subset/*"
        call(command,shell=True)

    else:
        print "No timesteps to remove in file: %s"%(filename)
    print "\n\n-----"