import paramiko
import os
import pandas as pd
import datetime as datetime

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2014, 2, 11)
__modified__ = datetime.datetime(2014, 2, 11)
__version__  = "1.0.0"
__status__   = "Production"


"""Start EDIT"""

myvars=['temp','salt']
mylevels=[0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500,
           600, 700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]

myprefix='regscen_hindcast'
apattern = 'AA_10km_z.nc_*'
outfile  = 'allfiles.nc'

myhostname="hexagon.bccs.uib.no"
myusername="trondk"
mypassword="passord"
remotedir="/work/users/trondk/REGSCEN/"
remoteSTORAGEdir="/work/jonal/AA_10km/monthly/"
localdir="/Users/trondkr/Projects/RegScen/Analysis/"
first=False

"""Stop EDIT"""

#
# Calculations done:
# - Monthly climatology of salinity and temperature at selected depths
# - timeseries of SST and SSS
# - trend of SST and SSS
#

ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(myhostname,username=myusername,password=mypassword)  

# Prepare the environment
command = 'cd %s'%(remotedir)
stdin, stdout, stderr = ssh.exec_command(command)

ftp = ssh.open_sftp()

for myvar in myvars:
    # Extract all data of variable 'myvar'
    command = 'cdo -select,name=%s %s%s %s%s'%(myvar,remoteSTORAGEdir,apattern,remotedir,outfile)
    stdin, stdout, stderr = ssh.exec_command(command)
    exit_status = stdout.channel.recv_exit_status()
  
    
    for mylevel in mylevels:
        print "Extracting data for depth level: %s"%(mylevel)
        
        # Extract data at seperate depth levels of variable 'myvar'
        outfileLevel  = str(myprefix)+'_'+str(myvar)+'_depth_'+str(mylevel)+'.nc'
        command = 'cdo sellevel,%s %s%s %s%s'%(mylevel,remotedir,outfile,remotedir,outfileLevel)
        stdin, stdout, stderr = ssh.exec_command(command)
        exit_status = stdout.channel.recv_exit_status()
        
        print "Starting statistical calculations"
            
        # Calculate monthly stat at seperate depth levels of variable 'myvar'
        outfileAverageLevel  = str(myprefix)+'_ymonavg_'+str(myvar)+'_depth_'+str(mylevel)+'.nc'
        command = 'cdo ymonavg %s%s %s%s'%(remotedir,outfileLevel,remotedir,outfileAverageLevel)
        stdin, stdout, stderr = ssh.exec_command(command)
        exit_status = stdout.channel.recv_exit_status()
        
        remotefile=remotedir+outfileAverageLevel
        print "Calculations done and downloading file: %s"%(remotefile)      
        localfile=localdir+outfileAverageLevel
        ftp.get(remotefile, localfile)
        print "\n"

        # Calculate trend of SST and SSS
        outfileAverageLevelB  = str(myprefix)+'_trend_'+str(myvar)+'_depth_'+str(mylevel)+'_trend.nc'
        outfileAverageLevelA  = str(myprefix)+'_trend_'+str(myvar)+'_depth_'+str(mylevel)+'_intercept.nc'
        command = 'cdo trend %s%s %s%s %s%s'%(remotedir,outfileLevel,remotedir,outfileAverageLevelA,remotedir,outfileAverageLevelB)
        stdin, stdout, stderr = ssh.exec_command(command)
        exit_status = stdout.channel.recv_exit_status()
        
        remotefile=remotedir+outfileAverageLevelB
        print "Calculations done and downloading file: %s"%(remotefile)      
        localfile=localdir+outfileAverageLevelB
        ftp.get(remotefile, localfile)
        print "\n"
             
        
ftp.close()
ssh.close()

print "Program finished"