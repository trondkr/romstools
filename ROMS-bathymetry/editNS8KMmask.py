import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

"""This script is to be used interactively in ipython"""
grd = pyroms.grid.get_ROMS_grid('NS8KM')

mymap = Basemap(llcrnrlon=-18.0,
                  llcrnrlat=46.0,
                  urcrnrlon=25.5,
                  urcrnrlat=67.5,
                  resolution='i',projection='tmerc',lon_0=0,lat_0=50,area_thresh=50.)

coast = pyroms.utility.get_coast_from_map(mymap)
pyroms.grid.edit_mask_mesh_ij(grd.hgrid, coast=coast)


"""After the edits have been done , run this command (interactively)"""
pyroms.grid.write_ROMS_grid(grd, filename='grid_py.nc')
