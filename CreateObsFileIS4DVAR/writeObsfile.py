
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os
import datetime

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 12, 30)
__modified__ = datetime.datetime(2013, 2, 12)
__version__  = "1.0"
__status__   = "Development, 12.2.2013"

"""When calling this function some of the variables are optional and will not be written until lastIteration=True"""

def writeData(outputFile,obs_lat,obs_lon,obs_value,Nobs,survey_time,obs_time,obs_Xgrid,obs_Ygrid,
                               firstIteration,lastIteration,
                               obs_flag,obs_type,obs_error,obs_Zgrid,obs_depth,obs_variance,
                               survey,is3d,Nstate,USENETCDF4):


   if USENETCDF4 is True:
      myZLIB=True
      myformat='NETCDF4'
   else:
      myZLIB=False
      myformat='NETCDF3_CLASSIC'

   if firstIteration is True:

      f1 = Dataset(outputFile, mode='w', format=myformat)
      f1.description="This is a obs file for SST"
      f1.history = 'Created ' + time.ctime(time.time())
      f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
      f1.type='NetCDF4 using program createMapNS.py'

      f1.options='Program requires: getCortad.py and writeObsfile.py'

      """ Define dimensions """
      f1.createDimension('one', 1)
      f1.createDimension('state_variable', Nstate)
      f1.createDimension('datum', None)

      v_spherical = f1.createVariable('spherical', 'c', ('one',),zlib=myZLIB)
      v_spherical.long_name = 'grid type logical switch'
      v_spherical.option_T  = "spherical"
      v_spherical.option_F  = "Cartesian"
      v_spherical[:]        = "T"

      v_obs_type = f1.createVariable('obs_type', 'i', ('datum',),zlib=myZLIB)
      v_obs_type.long_name = 'model state variable associated with observation'
      v_obs_type.opt_1 ='free-surface'
      v_obs_type.opt_2 ='vertically integrated u-momentum component';
      v_obs_type.opt_3 ='vertically integrated v-momentum component';
      v_obs_type.opt_4 ='u-momentum component'
      v_obs_type.opt_5 ='v-momentum component'
      v_obs_type.opt_6 ='potential temperature'
      v_obs_type.opt_7 ='salinity'
      v_obs_type[:]    = obs_type

      v_time = f1.createVariable('obs_time', 'd', ('datum',),zlib=myZLIB)
      v_time.long_name = 'Time of observation'
      v_time.units     = 'days'
      v_time.field     = 'time, scalar, series'
      v_time.calendar  = 'standard'
      v_time[:]        = obs_time


      v_obs_lon = f1.createVariable('obs_lon', 'd', ('datum',),zlib=myZLIB)
      v_obs_lon.long_name = 'Longitude of observation'
      v_obs_lon.units     = 'degrees_east'
      v_obs_lon.min       = -180
      v_obs_lon.max       = 180
      v_obs_lon[:]        = obs_lon

      v_obs_lat = f1.createVariable('obs_lat', 'd', ('datum',),zlib=myZLIB)
      v_obs_lat.long_name = 'Latitude of observation'
      v_obs_lat.units     = 'degrees_north'
      v_obs_lat.min       = -90
      v_obs_lat.max       = 90
      v_obs_lat[:]        = obs_lat

      v_obs_depth = f1.createVariable('obs_depth', 'd', ('datum',),zlib=myZLIB)
      v_obs_depth.long_name = 'Depth of observation'
      v_obs_depth.units     = 'meter'
      v_obs_depth.minus     = 'downwards'
      v_obs_depth[:]        = obs_depth

      v_obs_error = f1.createVariable('obs_error', 'd', ('datum',),zlib=myZLIB)
      v_obs_error.long_name = 'Observation error covariance'
      v_obs_error.units     = 'squared state variable units'
      v_obs_error[:]        = obs_error

      v_obs_val = f1.createVariable('obs_value', 'd', ('datum',),zlib=myZLIB)
      v_obs_val.long_name = 'Observation value'
      v_obs_val.units     = 'state variable units'
      v_obs_val[:]        = obs_value

      v_obs_xgrid = f1.createVariable('obs_Xgrid', 'd', ('datum',),zlib=myZLIB)
      v_obs_xgrid.long_name = 'x-grid observation location'
      v_obs_xgrid.units     = 'nondimensional'
      v_obs_xgrid.left      = "INT(obs_Xgrid(datum))"
      v_obs_xgrid.right     = "INT(obs_Xgrid(datum))+1"
      v_obs_xgrid[:]        = obs_Xgrid

      v_obs_ygrid = f1.createVariable('obs_Ygrid', 'd', ('datum',),zlib=myZLIB)
      v_obs_ygrid.long_name = 'y-grid observation location'
      v_obs_ygrid.units     = 'nondimensional'
      v_obs_ygrid.top       = "INT(obs_Ygrid(datum))+1"
      v_obs_ygrid.bottom    = "INT(obs_Ygrid(datum))"
      v_obs_ygrid[:]        = obs_Ygrid

      v_obs_zgrid = f1.createVariable('obs_Zgrid', 'd', ('datum',),zlib=myZLIB)
      v_obs_zgrid.long_name = 'z-grid observation location'
      v_obs_zgrid.units     = 'nondimensional'
      v_obs_zgrid.up        = "INT(obs_Zgrid(datum))+1"
      v_obs_zgrid.down      = "INT(obs_Zgrid(datum))"
      v_obs_zgrid[:]        = obs_Zgrid

      f1.close()

   if firstIteration is False and lastIteration is False:

      f1 = Dataset(outputFile, mode='a', format=myformat)

      t0 = time.time()
      """Find index for ading new info to arrays (same for all variables)"""
      myshape=f1.variables["obs_Zgrid"][:].shape
      indexStart=myshape[0]
      indexEnd=obs_Zgrid.shape[0]+myshape[0]

      f1.variables["obs_type"][indexStart:indexEnd] = obs_type
      f1.variables["obs_time"][indexStart:indexEnd] = obs_time
      f1.variables["obs_lon"][indexStart:indexEnd] = obs_lon
      f1.variables["obs_lat"][indexStart:indexEnd] = obs_lat
      f1.variables["obs_depth"][indexStart:indexEnd] = obs_depth
      f1.variables["obs_error"][indexStart:indexEnd] = obs_error
      f1.variables["obs_value"][indexStart:indexEnd] = obs_value
      f1.variables["obs_Xgrid"][indexStart:indexEnd] = obs_Xgrid
      f1.variables["obs_Ygrid"][indexStart:indexEnd] = obs_Ygrid
      f1.variables["obs_Zgrid"][indexStart:indexEnd] = obs_Zgrid
    
      t1 = time.time()
      print "array append created in %s seconds"%(t1-t0)
      f1.close()

   if lastIteration is True:

      f1 = Dataset(outputFile, mode='a', format=myformat)

      f1.createDimension('survey', survey)

      v_obs = f1.createVariable('Nobs', 'i', ('survey',),zlib=myZLIB)
      v_obs.long_name = 'number of observations with the same survey time'
      v_obs.field     = 'scalar, series'
      v_obs[:]        = Nobs

      v_time = f1.createVariable('survey_time', 'i', ('survey',),zlib=myZLIB)
      v_time.long_name = 'Survey time'
      v_time.units     = 'day'
      v_time.field     = 'time, scalar, series'
      v_time.calendar  = 'standard'
      v_time[:]        = survey_time

      v_obs_var = f1.createVariable('obs_variance', 'd', ('state_variable',),zlib=myZLIB)
      v_obs_var.long_name = 'global time and space observation variance'
      v_obs_var.units     = 'squared state variable units'
      v_obs_var[:]        = obs_variance

      f1.close()
