import paramiko
import os
import pandas as pd
import datetime as datetime

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2014, 2, 11)
__modified__ = datetime.datetime(2014, 3, 17)
__version__  = "1.0.0"
__status__   = "Production"


"""Start EDIT"""

myvars=['temp','salt']

myprefix='regscen_hindcast'
apattern = 'AA_10km_z.nc_*'
outfile2  = 'allfiles.nc'
outfile  = 'allfiles_region.nc'

myhostname="hexagon.bccs.uib.no"
myusername="trondk"
mypassword="passord"
remotedir="/work/users/trondk/REGSCEN/"
remoteSTORAGEdir="/work/jonal/AA_10km/monthly/"
localdir="/Users/trondkr/Projects/RegScen/Analysis/"
first=False

"""

    This script calculates the domain averages of salinity and temperature at all depth
    levels available. These timeseries will be weighted later to calculate the domain averaged
    monthly fields within depth layers. In addition, timeseries and domain averages of SSH and SHFLUX
    is calculated togetehr with the trend of SSH (zeta).

    Calculations done:
    - Domain averages of salinity and temperature at all depth levels
    - timeseries and domain averages of zeta and shflux
    - trend of zeta (sea level height)

"""

alldepths=[0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500,
           600, 700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]

"""Stop EDIT"""
    
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(myhostname,username=myusername,password=mypassword)  

# Prepare the environment
command = 'cd %s'%(remotedir)
stdin, stdout, stderr = ssh.exec_command(command)

ftp = ssh.open_sftp()

command = 'find %s -name %s'%(remoteSTORAGEdir,apattern)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

filelist = stdout.read().splitlines()

for myvar in myvars:
    
    # Extract all data of variable 'myvar'
    command = 'cdo -select,name=%s %s%s %s%s '%(myvar,remoteSTORAGEdir,apattern,remotedir,outfile)
    stdin, stdout, stderr = ssh.exec_command(command)
    exit_status = stdout.channel.recv_exit_status()
    
    
    for mylevel in alldepths:
        mylevel=str(mylevel)
        print "Extracting data for depth level: %s"%(mylevel)
        
        # Extract data at seperate depth levels of variable 'myvar'
        outfileLevel  = str(myprefix)+'_'+str(myvar)+'_depth_'+str(mylevel)+'.nc'
        command = 'cdo sellevel,%s %s%s %s%s'%(mylevel,remotedir,outfile,remotedir,outfileLevel)
        stdin, stdout, stderr = ssh.exec_command(command)
        exit_status = stdout.channel.recv_exit_status()
             
        # Calculate domain average at seperate depth levels of variable 'myvar'
        outfileAverageLevel  = str(myprefix)+'_fldmean_'+str(myvar)+'_depth_'+str(mylevel)+'.nc'
        command = 'cdo fldmean %s%s %s%s'%(remotedir,outfileLevel,remotedir,outfileAverageLevel)
        stdin, stdout, stderr = ssh.exec_command(command)
        exit_status = stdout.channel.recv_exit_status()
        
        remotefile=remotedir+outfileAverageLevel
        print "Calculations done and downloading file: %s"%(remotefile)      
        localfile=localdir+outfileAverageLevel
        ftp.get(remotefile, localfile)
        print "\n"
        

# Calculate trend of SLA
# Extract all data of variable zeta
command = 'cdo -sellonlatbox,-30,15,60,90 -select,name=zeta %s%s %s%s'%(remoteSTORAGEdir,apattern,remotedir,outfile)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

outfileAverageLevelB  = str(myprefix)+'_trend_zeta_trend.nc'
outfileAverageLevelA  = str(myprefix)+'_trend_zeta_intercept.nc'
command = 'cdo trend %s%s %s%s %s%s'%(remotedir,outfile,remotedir,outfileAverageLevelA,remotedir,outfileAverageLevelB)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

remotefile=remotedir+outfileAverageLevelB
print "ZETA Calculations done and downloading file: %s"%(remotefile)      
localfile=localdir+outfileAverageLevelB
ftp.get(remotefile, localfile)
print "\n"
    
# Calculate domain average of ZETA
outfileAverage  = str(myprefix)+'_fldmean_zeta.nc'
command = 'cdo fldmean %s%s %s%s'%(remotedir,outfile,remotedir,outfileAverage)

stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

remotefile=remotedir+outfileAverage
print "Calculations done and downloading file: %s"%(remotefile)      
localfile=localdir+outfileAverage
ftp.get(remotefile, localfile)
print "\n"    
       

# SHFLUX
# Extract all data of variable ssflux
command = 'cdo -sellonlatbox,-30,15,60,90 -select,name=shflux %s%s %s%s'%(remoteSTORAGEdir,apattern,remotedir,outfile)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

# Calculate domain average of shflux - timeseries
outfileAverage  = str(myprefix)+'_fldmean_shflux.nc'
command = 'cdo fldmean %s%s %s%s'%(remotedir,outfile,remotedir,outfileAverage)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

remotefile=remotedir+outfileAverage
localfile=localdir+outfileAverage
ftp.get(remotefile, localfile)
print "\n"

# Calculate time average of surface fresh water flux areamap
outfileAverage  = str(myprefix)+'_timmean_shflux.nc'
command = 'cdo timmean %s%s %s%s'%(remotedir,outfile,remotedir,outfileAverage)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

remotefile=remotedir+outfileAverage
print "SHFLUX Calculations done and downloading file: %s"%(remotefile)      
localfile=localdir+outfileAverage
ftp.get(remotefile, localfile)
print "\n"              
                      
ftp.close()
ssh.close()

print "Program finished"