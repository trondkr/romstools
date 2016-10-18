import os
from numpy import *
from matplotlib.pyplot import *
from netCDF4 import Dataset
from pylab import *
import pyroms
import pyroms_toolbox
from mpl_toolkits.basemap import Basemap, shiftgrid
import mpl_toolkits.basemap as mp
import mpl_util
import string

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2015, 7, 19)
__modified__ = datetime.datetime(2015, 7, 19)
__version__  = "1.0"
__status__   = "Development, 19.7.2015"

doc=""" 

Adds changes created using editMask.py saved to mask_change.txt 

"""

f=Dataset("kino_1600m_20072015.nc", mode='a')
mask_rho=f.variables["mask_rho"][:,:]

infile="/Users/trondkr/Projects/KINO/GRID/editgrid/mask_change.txt"
fmask=open(infile,"r")

lines=fmask.readlines()
for line in lines:
    l=string.split(line," ")
  
    i=int(float(l[0].strip()))
    j=int(float(l[1].strip()))
    m=float(l[2].strip())
    print "Changing %s %s from %s to %s"%(j,i,mask_rho[j,i],m)
    mask_rho[j,i]=m

print np.shape(f.variables["mask_rho"][:,:]),np.shape(mask_rho)
f.variables["mask_rho"][:,:]=mask_rho

pcolor(mask_rho)
plt.show()
f.close()
