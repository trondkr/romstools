# encoding: utf-8

import sys
import os
import numpy as np
from mpl_toolkits.basemap import Basemap
from datetime import datetime
import netCDF4 as netCDF

import pyroms
from pyroms.hgrid import *
from pyroms.vgrid import *
from pyroms.grid import *
from pyroms import io

def write_ROMS_grid(grd, visc_factor, diff_factor, filename='roms_grd.nc'):
    """
    write_ROMS_grid(grd, filename)

    Write ROMS_CGrid class on a NetCDF file.
    """

    Mm, Lm = grd.hgrid.x_rho.shape

    
    # Write ROMS grid to file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF4')
    nc.Description = 'ROMS grid'
    nc.Author = 'Trond Kristiansen'
    nc.Created = datetime.now().isoformat()
    nc.type = 'ROMS grid file'

    nc.createDimension('xi_rho', Lm)
    nc.createDimension('xi_u', Lm-1)
    nc.createDimension('xi_v', Lm)
    nc.createDimension('xi_psi', Lm-1)
    
    nc.createDimension('eta_rho', Mm)
    nc.createDimension('eta_u', Mm)
    nc.createDimension('eta_v', Mm-1)
    nc.createDimension('eta_psi', Mm-1)

    nc.createDimension('xi_vert', Lm+1)
    nc.createDimension('eta_vert', Mm+1)

    nc.createDimension('bath', None)

    if hasattr(grd.vgrid, 's_rho') is True and grd.vgrid.s_rho is not None:
        N, = grd.vgrid.s_rho.shape
        nc.createDimension('s_rho', N)
        nc.createDimension('s_w', N+1)

    def write_nc_var(var, name, dimensions, long_name=None, units=None):
        nc.createVariable(name, 'f8', dimensions)
        if long_name is not None:
            nc.variables[name].long_name = long_name
        if units is not None:
            nc.variables[name].units = units
        nc.variables[name][:] = var
        print ' ... wrote ', name

    if hasattr(grd.vgrid, 's_rho') is True and grd.vgrid.s_rho is not None:
        write_nc_var(grd.vgrid.theta_s, 'theta_s', (), 'S-coordinate surface control parameter')
        write_nc_var(grd.vgrid.theta_b, 'theta_b', (), 'S-coordinate bottom control parameter')
        write_nc_var(grd.vgrid.Tcline, 'Tcline', (), 'S-coordinate surface/bottom layer width', 'meter')
        write_nc_var(grd.vgrid.hc, 'hc', (), 'S-coordinate parameter, critical depth', 'meter')
        write_nc_var(grd.vgrid.s_rho, 's_rho', ('s_rho'), 'S-coordinate at RHO-points')
        write_nc_var(grd.vgrid.s_w, 's_w', ('s_w'), 'S-coordinate at W-points')
        write_nc_var(grd.vgrid.Cs_r, 'Cs_r', ('s_rho'), 'S-coordinate stretching curves at RHO-points')
        write_nc_var(grd.vgrid.Cs_w, 'Cs_w', ('s_w'), 'S-coordinate stretching curves at W-points')

    write_nc_var(grd.vgrid.h, 'h', ('eta_rho', 'xi_rho'), 'bathymetry at RHO-points', 'meter')
    #ensure that we have a bath dependancy for hraw
    if len(grd.vgrid.hraw.shape) == 2:
        hraw = np.zeros((1, grd.vgrid.hraw.shape[0], grd.vgrid.hraw.shape[1]))
        hraw[0,:] = grd.vgrid.hraw
    else:
        hraw = grd.vgrid.hraw
    write_nc_var(hraw, 'hraw', ('bath', 'eta_rho', 'xi_rho'), 'raw bathymetry at RHO-points', 'meter')
    write_nc_var(grd.hgrid.f, 'f', ('eta_rho', 'xi_rho'), 'Coriolis parameter at RHO-points', 'second-1')
    write_nc_var(1./grd.hgrid.dx, 'pm', ('eta_rho', 'xi_rho'), 'curvilinear coordinate metric in XI', 'meter-1')
    write_nc_var(1./grd.hgrid.dy, 'pn', ('eta_rho', 'xi_rho'), 'curvilinear coordinate metric in ETA', 'meter-1')
    write_nc_var(grd.hgrid.dmde, 'dmde', ('eta_rho', 'xi_rho'), 'XI derivative of inverse metric factor pn', 'meter')
    write_nc_var(grd.hgrid.dndx, 'dndx', ('eta_rho', 'xi_rho'), 'ETA derivative of inverse metric factor pm', 'meter')
    write_nc_var(grd.hgrid.xl, 'xl', (), 'domain length in the XI-direction', 'meter')
    write_nc_var(grd.hgrid.el, 'el', (), 'domain length in the ETA-direction', 'meter')

    write_nc_var(grd.hgrid.x_rho, 'x_rho', ('eta_rho', 'xi_rho'), 'x location of RHO-points', 'meter')
    write_nc_var(grd.hgrid.y_rho, 'y_rho', ('eta_rho', 'xi_rho'), 'y location of RHO-points', 'meter')
    write_nc_var(grd.hgrid.x_u, 'x_u', ('eta_u', 'xi_u'), 'x location of U-points', 'meter')
    write_nc_var(grd.hgrid.y_u, 'y_u', ('eta_u', 'xi_u'), 'y location of U-points', 'meter')
    write_nc_var(grd.hgrid.x_v, 'x_v', ('eta_v', 'xi_v'), 'x location of V-points', 'meter')
    write_nc_var(grd.hgrid.y_v, 'y_v', ('eta_v', 'xi_v'), 'y location of V-points', 'meter')
    write_nc_var(grd.hgrid.x_psi, 'x_psi', ('eta_psi', 'xi_psi'), 'x location of PSI-points', 'meter')
    write_nc_var(grd.hgrid.y_psi, 'y_psi', ('eta_psi', 'xi_psi'), 'y location of PSI-points', 'meter')
    write_nc_var(grd.hgrid.x_vert, 'x_vert', ('eta_vert', 'xi_vert'), 'x location of cell verticies', 'meter')
    write_nc_var(grd.hgrid.y_vert, 'y_vert', ('eta_vert', 'xi_vert'), 'y location of cell verticies', 'meter')

    if hasattr(grd.hgrid, 'lon_rho'):
        write_nc_var(grd.hgrid.lon_rho, 'lon_rho', ('eta_rho', 'xi_rho'), 'longitude of RHO-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_rho, 'lat_rho', ('eta_rho', 'xi_rho'), 'latitude of RHO-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_u, 'lon_u', ('eta_u', 'xi_u'), 'longitude of U-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_u, 'lat_u', ('eta_u', 'xi_u'), 'latitude of U-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_v, 'lon_v', ('eta_v', 'xi_v'), 'longitude of V-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_v, 'lat_v', ('eta_v', 'xi_v'), 'latitude of V-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_psi, 'lon_psi', ('eta_psi', 'xi_psi'), 'longitude of PSI-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_psi, 'lat_psi', ('eta_psi', 'xi_psi'), 'latitude of PSI-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_vert, 'lon_vert', ('eta_vert', 'xi_vert'), 'longitude of cell verticies', 'degree_east')
        write_nc_var(grd.hgrid.lat_vert, 'lat_vert', ('eta_vert', 'xi_vert'), 'latitude of cell verticies', 'degree_north')

    nc.createVariable('spherical', 'c')
    nc.variables['spherical'].long_name = 'Grid type logical switch'
    nc.variables['spherical'][:] = grd.hgrid.spherical
    print ' ... wrote ', 'spherical'

    write_nc_var(grd.hgrid.angle_rho, 'angle', ('eta_rho', 'xi_rho'), 'angle between XI-axis and EAST', 'radians')

    write_nc_var(grd.hgrid.mask_rho, 'mask_rho', ('eta_rho', 'xi_rho'), 'mask on RHO-points')
    write_nc_var(grd.hgrid.mask_u, 'mask_u', ('eta_u', 'xi_u'), 'mask on U-points')
    write_nc_var(grd.hgrid.mask_v, 'mask_v', ('eta_v', 'xi_v'), 'mask on V-points')
    write_nc_var(grd.hgrid.mask_psi, 'mask_psi', ('eta_psi', 'xi_psi'), 'mask on psi-points')

    if visc_factor != None:
        write_nc_var(visc_factor, 'visc_factor', ('eta_rho', 'xi_rho'), 'horizontal viscosity sponge factor')
    if diff_factor != None:
        write_nc_var(diff_factor, 'diff_factor', ('eta_rho', 'xi_rho'), 'horizontal diffusivity sponge factor')
    
    nc.close()
