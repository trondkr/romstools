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
mypassword="password"
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
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(myhostname,username=myusername,password=mypassword)  

"""Prepare the environment"""
command = 'cd %s'%(remotedir)

stdin, stdout, stderr = ssh.exec_command(command)


ftp = ssh.open_sftp()
        
season="all"      
depth="40"

command = 'find %s -name %s'%(remoteSTORAGEdir,apattern)
stdin, stdout, stderr = ssh.exec_command(command)
exit_status = stdout.channel.recv_exit_status()

filelist = stdout.read().splitlines()
filelist.sort()
print "\STARTING extracting data:\n"
for myfile in filelist:
	head,tail = os.path.split(myfile)
	
	for myvar in ['u_eastward','v_northward']:
	
		"""Extract all data of variable u_eastward"""
		print "Variable: %s"%(myvar)
		command = 'cdo -seltimestep,1/24 -select,name=%s %s %s%s'%(myvar,myfile,remotedir,"tmp")
		print command
		stdin, stdout, stderr = ssh.exec_command(command)
		exit_status = stdout.channel.recv_exit_status()

		outfileLevel  = str(tail[0:-3])+'_'+str(myvar)+'_depth_'+str(depth)+'.nc'
		command='rm %s'%(outfileLevel)
		stdin, stdout, stderr = ssh.exec_command(command)
		exit_status = stdout.channel.recv_exit_status()
	       
		command = 'cdo -sellevidx,%s  %s%s %s%s'%(depth,remotedir,"tmp",remotedir,outfileLevel)
		print command
		stdin, stdout_cdo, stderr = ssh.exec_command(command)
		exit_status = stdout_cdo.channel.recv_exit_status()
	     
		remotefile=remotedir+outfileLevel
		print "%s Calculations done and downloading file: %s"%(myvar,remotefile)      
		localfile=localdir+outfileLevel
		print remotefile, localfile
		if os.path.exists(localfile): os.remove(localfile)
		ftp.get(remotefile, localfile)
		print "\n"

		# Clean up
		command = 'rm %s'%(remotefile)
		print command
		stdin, stdout_cdo, stderr = ssh.exec_command(command)
		exit_status = stdout_cdo.channel.recv_exit_status()
	    
		outfileLevel=""
		remotefile=""
		localfile=""


ftp.close()
ssh.close()

print "Program finished"