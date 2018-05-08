Romstools
=========
[![GitHub version](https://badge.fury.io/gh/trondkr%2Fromstools.svg)](http://badge.fury.io/gh/trondkr%2Fromstools)

Romstools is a toolbox that contains a variety of programs useful for Regional Ocean Model (ROMS) developers. Our goal with this and other toolboxes available on Github is to provide everything required to setup, run, and analyse the ROMS model (ROMS - http://myroms.org/). Currently, the following tools are included but split into separate packages/folders that can be used indepdentely.

<ul>
<li><a href="https://github.com/trondkr/romstools/tree/master/CreateObsFileIS4DVAR"><strong>CreateObsFileIS4DVAR</strong></a> - Generate observation file from SST required to run IS4DVAR assimilation with ROMS </li>
  
<li><a href="https://github.com/trondkr/romstools/tree/master/VolumeFlux"><strong>VolumeFlux</strong></a> - ROMS volume flux calculations for transects between (x,y) grid points </li>

<li><a href="https://github.com/trondkr/model2roms"><strong>model2roms</strong></a> - Easily create BRY, INIT ,and CLIM files required to run ROMS using the model2roms toolbox </li>

<li><a href="https://github.com/trondkr/romstools/tree/master/Tools"><strong>Tools</strong></a> - Selection of useful scripts that can come in handy when working with NetCDF files and ROMS </li>

<li><a href="https://github.com/trondkr/romstools/tree/master/create_atmos_ROMS"><strong>create_atmos_ROMS</strong></a> - Create ERA INTERIM forcing for ROMS </li>

<li><a href="https://github.com/trondkr/romstools/tree/master/create_tides_ROMS"><strong>create_tides_ROMS</strong></a>  - Create tidal forcing for ROMS using TPXO9.1 </li>

<li><a href="https://github.com/trondkr/romstools/tree/master/Postpro"><strong>Postpro</strong></a> - Program for converting all of your sigma level ROMS output files into Z-level files </li>
</ul>


<h3> General requirements</h3>
<ul>
<li>Python installation with numpy, netcdf4, basemap, and matplotlib</li>
<li>In some cases a Fortran compiler combined with F2PY (part of numpy) is required to create python modules (gfortran is great)</li>
<li>A full Python distribution package such as <a href="https://store.continuum.io/cshop/anaconda/">Anaconda</a> is reccomended</li>
<li><a href="https://www.unidata.ucar.edu/downloads/netcdf/index.jsp">Fortran</a> NetCDF interface (required to compile create_tides)</li>
<li> To run interaxctively on Hexagon (HPC) remember to: <b>module unload xtpe-interlagos</b> </li>
</ul>

<h3>Tools</h3>

The tools section contains various programs to calculate long term averages, and trends based on a number of input
files that typically cover a long time period. For example, if you want to analyse the modeled multi-year averages
of short wave fluxes for the North Sea you can do so by running the <em>createAveragesOnHexagon.py</em> script.
The script first calculates the averages on the supercomputer before downloading the resulting NetCDF4 files to
local computer. Finally, you can plot the results stored in the NetCDF4 file by running <em>createMaps.py</em>. The result
will look something like this:
![Shflux averaged 1993-2009](http://www.trondkristiansen.com/wp-content/gallery/romstools/longtermmean_shflux_time_depth_surface.jpg)

<h3> Generate observation file for IS4DVAR </h3>

This is a set of scripts that automatically downloads Sea Surface Temperature (SST) data for your grid for the time period you are interested in and stored the data as an observation file. This file requires a special format and is read by ROMS when assimilating SST data using the Incremental Strong constraint 4DVAR assimilation technique. The SST data are downloaded using openDAP from the CoRTAD SST archive (http://www.nodc.noaa.gov/sog/Cortad/) or the OI AVHRR SST (http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html).

<h3> Postpro </h3>

This is a very useful set of tools for interpolating all of your ROMS results files into Z-level. To use this package: 

<ul>
<li>Create symbolic links in the Files4Stat folder to your result files</li>
<li>Edit the execute_Postpro.sh and execute it to generate input file for roms2z ./execute_Postpro.sh</li>
<li>Run the ROMS_Postpro.job script (if not on super computer create a script file out of this job script)</li>
<li>Split the resulting large z-level-file into monthly average files by running the following: </li>
</ul>
```bash
cd /work/trondk
./ncfile_subsetting.sh 1948 01 01 northsea_8km_z.nc time d
for (( yy=1993; yy<=2013; yy++ )); do
for mm in 01 02 03 04 05 06 07 08 09 10 11 12; do
ncrcat -O northsea_8km_z.nc_${yy}${mm}* northsea_8km_z.nc_${yy}${mm}
cdo timmean northsea_8km_z.nc_${yy}${mm} northsea_8km_z_mean.nc_${yy}${mm}
rm northsea_8km_z.nc_${yy}${mm}
echo “Finished with year $yy and month $mm”
done
done
rm northsea_8km_z.nc_${yy}${mm}????
```

<h3> Volume flux calculations </h3>

To calculate the fluxes you run the script <em>calculateFluxes.py</em> either as a standalone python program or executed
as a job-script. If you run the program on a super-computer (e.g. Hexagon) you may want to run the script as a job
and use the script <em>runJob.sh</em> to queue your job (Hexagon: qsub runJob.sh).

This toolbox is a mixture of Fortrand and Python tools where the core programs are taken from the excellent pyroms
toolbox also available on github: https://github.com/kshedstrom/pyroms. However, to avoid having to install the
entire toolbox to calculate the fluxes, this smaller toolbox was created.

![Calculated volume fluxes for various sections in the North Sea between 1993-2009](http://www.trondkristiansen.com/wp-content/gallery/romstools/ns8km_vflux_volumeflux.png)

<h4> Define transects and depth ranges </h4>
You can define a list of transects you want the volume transport calculated for in the function
<em>defineTransects()</em>.  The output from running is a comma separated value file containing the positive,
negative, and net transport through the transect. You can also calculate the transport for the e.g. just the upper
500 meters of the water column by defining the minimum and maximum depths:

```Python
minDepth=0
maxDepth=500
```

<h4> Requirements </h4>
This toolbox still requires you to have a Fortran compiler to compile the Fortran programs and generate Python modules.
On Hexagon this is done by loading the gnu modules (for Fortran compiler compatible with Python and numpy). In the
terminal window type (note: you do not issue the two first commands 'module swap' and 'module unload' unless you work on super computers):

```bash
module swap PrgEnv-pgi PrgEnv-gnu
module unload notur
f2py --verbose  -c -m iso iso.f90
f2py --verbose  -c -m obs_interp obs_interp.f
```
This should provide you with two python modules (obs_interp.so and iso.so) which you can try to import to python with:

```bash
Python 2.7.2 (default, Mar 22 2012, 12:32:11)
[GCC 4.6.1 20110627 (Cray Inc.)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import iso
>>> print iso.__doc__
```

<h3> Create Atmospheric forcing - create_atmos_ROMS</h3>
This toolset contains necessary files to create forcing files for the atmosphere based on ERA INTERIM data. 

Get the toolbox and move into folder ```create_atmos_ROMS```, edit the file ```make_atmos.sh``` so that the following variables are correct according to your setup:
``` bash
gridfile=/work/shared/nn9297k/FAABolous/WS4KM_grd.nc
atmeraintdir=/work/shared/norkyst/NorKyst-800m_Forcing/Atmos/ERA_Interim
for ((year=1980; year<=2015; year++)); 
```
Next compile and run with  ```./make_atmos.sh``` 

<h3> Create Tidal forcing - create_tides_ROMS</h3>
This toolset contains necessary files to create tidal forcing files for your ROMS grid. Start by editing the ```make_tides.sh``` file to point to the correct grid, appname, and time eperiod you wnat forcing created for.

``` bash
# Start and end period for simulation
syear=1980; eyear=2017

# Input files
gridfile=/Users/trondkr/Dropbox/NIVA/ROHO800/Grid/ROHO800_grid.nc

# Name of app
APPNAME=ROHO800
```

Make sure you download the latest <a><href="http://volkov.oce.orst.edu/tides/global.html"><TPXO9.1></a> files from Oregon State University (OSU) before you continue. Edit the path for the data input file inside the file ```Tpxo/tpxo2grid.f90``` to direct to teh path where you stored the TPXO9.1 files (grid_tpxo9.nc, h_tpxo9.v1.nc, u_tpxo9.v1.nc):
  
``` fortran
character (len=100)       :: tpxo_grd = "/Volumes/DATASETS/tpxo9_netcdf/grid_tpxo9.nc"
character (len=100)       :: tpxo_h   = "/Volumes/DATASETS/tpxo9_netcdf/h_tpxo9.v1.nc"
character (len=100)       :: tpxo_u   = "/Volumes/DATASETS/tpxo9_netcdf/u_tpxo9.v1.nc"
``` Compilation works well on local computer with gfortran compiler (can be installed using brew). 

<h4> River forcing </h4>
Move into folder ```northsea_forcing_rivers```, edit the file ```make_rivers.sh``` so that the following variables are correct according to your setup:
``` bash
years="2009 2010 2011 2012"
gridfile=/work/users/trondk/KINO/GRID/kino_norseas_800m_grid.nc
```
Next compile and run with  ```./make_atmos.sh``` or use the batch job script   ```qsub runMakeAtmos.sh```

<h3> Contact </h3>

<ul>
<li>me (at) trondkristiansen.com</li>
</ul>


<h2>License</h2>
The MIT License (MIT)

Copyright (c) <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




