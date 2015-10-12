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

myprefix = 'ocean_avg_'
apattern = 'ocean_avg_00*'
outfile  = 'allfiles.nc'

myhostname="hexagon.bccs.uib.no"
myusername="trondk"
mypassword="ivNavac92171"
remotedir="/work/users/trondk/KINO/FORWARD/Run/"
remoteSTORAGEdir="/work/users/trondk/KINO/FORWARD/Run/"
localdir="/Users/trondkr/Projects/KINO/Analysis/"
first=False

"""This script calculates the domain averages of u_eastward an v_northward velocity at specified depth
level (depth="10").

Calculations done:
 - timeseries of u_eastward and v_northward velocity at defined depth level
 - climatology of u_eastward and v_northward currents at 20 m depth

NOTE:
If you don't have u_eastward and v_northward components you can use u and v velcoties and rotate the currents with the
angle variable afterwards. You have to change the variable names from u_eastward to u and v_northward to v below.

"""
os.popen("rm northsea_8km_timmean_u_eastward.nc")
os.popen("rm northsea_8km_timmean_v_northward.nc")
os.popen("ls -lrth")
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(myhostname,username=myusername,password=mypassword)  

"""Prepare the environment"""
command = 'cd %s'%(remotedir)
stdin, stdout, stderr = ssh.exec_command(command)

ftp = ssh.open_sftp()
        
season="all"      
depth="40"

outfile='ueastward_allfiles.nc'
"""Extract all data of variable u_eastward"""
command = 'cdo -select,name=u_eastward %s%s %s%s'%(remoteSTORAGEdir,apattern,remotedir,outfile)
print command
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

"""Extract u_eastward at depth"""
outfileLevel  = str(myprefix)+'_u_eastward_depth_'+str(depth)+'_season_'+str(season)+'.nc'
command='rm %s'%(outfileLevel)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()
       
command = 'cdo -sellevel,%s  %s%s %s%s'%(depth,remotedir,outfile,remotedir,outfileLevel)
print command
stdin, stdout_cdo, stderr = ssh.exec_command(command)
exit_status = stdout_cdo.channel.recv_exit_status()
        

"""Calculate domain average of u_eastward at depth - timeseries"""
outfileAverage  = str(myprefix)+'_timmean_u_eastward_depth_'+str(depth)+'_season_'+str(season)+'.nc'
command='rm %s'%(outfileAverage)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()
    
command = 'cdo timmean %s%s %s%s'%(remotedir,outfileLevel,remotedir,outfileAverage)
print command
stdin, stdout_cdo, stderr = ssh.exec_command(command)
exit_status = stdout_cdo.channel.recv_exit_status()

remotefile=remotedir+outfileAverage
print "U_EASTWARD Calculations done and downloading file: %s"%(remotefile)      
localfile=localdir+outfileAverage
print remotefile, localfile
if os.path.exists(localfile): os.remove(localfile)
ftp.get(remotefile, localfile)
print "\n"

""" V_NORTHWARD """                
outfile='v_northward_allfiles.nc'
"""Extract all data of variable u_eastward"""
command = 'cdo -select,name=v_northward %s%s %s%s'%(remoteSTORAGEdir,apattern,remotedir,outfile)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()


"""Extract v_northward at 20 m depth"""
outfileLevel  = str(myprefix)+'_v_northward_depth_'+str(depth)+'_season_'+str(season)+'.nc'
command='rm %s'%(outfileLevel)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()
  
command = 'cdo -sellevel,%s  %s%s %s%s'%(depth,remotedir,outfile,remotedir,outfileLevel)
stdin, stdout_cdo, stderr = ssh.exec_command(command)
exit_status = stdout_cdo.channel.recv_exit_status()
        
"""Calculate domain average of v_northwardat 20 m depth - timeseries"""
outfileAverage  = str(myprefix)+'_timmean_v_northward_depth_'+str(depth)+'_season_'+str(season)+'.nc'
command='rm %s'%(outfileAverage)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()
  
command = 'cdo timmean %s%s %s%s'%(remotedir,outfileLevel,remotedir,outfileAverage)
stdin, stdout_cdo, stderr = ssh.exec_command(command)
exit_status = stdout_cdo.channel.recv_exit_status()

remotefile=remotedir+outfileAverage
print "V_NORTHWARD Calculations done and downloading file: %s"%(remotefile)      
localfile=localdir+outfileAverage
if os.path.exists(localfile): os.remove(localfile)
ftp.get(remotefile, localfile)
print "\n"
                
        
ftp.close()
ssh.close()

print "Program finished"