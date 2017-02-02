from netCDF4 import Dataset, datetime, date2num,num2date
import numpy as np
import os
from subprocess import call

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2017, 2, 1)
__modified__ = datetime(2017, 2, 1)
__version__ = "1.0"
__status__ = "Development, modified on 01.02.2015"

myprefix='ws4km'
myvars=['Pair','Uwind','Vwind','cloud','Tair','Qair','rain']
years=[year for year in xrange(1989,2015,1)]

command = "export SKIP_SAME_TIME=1"
call(command,shell=True)

for myvar in myvars:
	newfilename="%s_%s_%s_%s.nc"%(myprefix,myvar,years[0],years[-1])
	mylistoffiles=[]
	if os.path.exists(newfilename): os.remove(newfilename)
	for year in years:
		filename="%s_%s_%s.nc"%(myprefix,myvar,year)
		mylistoffiles.append(filename)
		

	if len(mylistoffiles)>0:
		options=""
		for f in mylistoffiles:
			options=options+"%s "%(f)
		
		options=options+newfilename

		command = "cdo mergetime %s"%(options)
#		print "%s"%(command)
		call(command,shell=True)
		print "Result written to file: %s"%(newfilename)
		print "\nFinished with %s -----"%(myvar)
