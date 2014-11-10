Romstools
=========

Romstools is a toolbox that contains a variety of programs useful for Regional Ocean Model (ROMS) developers. Our goal with this and other toolboxes available on Github is to provide everything required to setup, run, and analyse the ROMS model (ROMS - http://myroms.org/). Currently, the following tools are included but split into separate packages/folders that can be used indepdentely.

<ul>
<li><strong>VolumeFlux</strong> - volume flux calculations for transects between (x,y) grid points</li>
<li><strong> Model2roms</strong> - automatically create BRY, INIT ,and CLIM files required to run ROMS using the model2roms toolbox - https://github.com/trondkr/model2roms</li>
<li><strong>Tools</strong> - a selection of useful scripts that can come in handy when working with NetCDF files and ROMS </li>
<li><strong>CreateForcing-Atmos-Tides-Rivers</strong> - programs that create atmospheric, river, and tidal forcing for your ROMS model.
</ul>


<h3> General requirements</h3>
<ul>
<li>Python installation with numpy, basemap, and matplotlib</li>
<li>Python netCDF4 interface - http://code.google.com/p/netcdf4-python/</li>
<li>In some cases a Fortran compiler combined with F2PY (part of numpy) is required to create python modules</li>
<li>A full Python distribution package such as <a href="https://store.continuum.io/cshop/anaconda/">Anaconda</a> or
<a href="https://www.enthought.com/">Enthought</a> is reccomended</li>
<li>Fortran NetCDF interface (required to compile CreateForcing-Atmos-Tides-Rivers)</li>
</ul>

<h3>Tools</h3>

The tools section contains various programs to calculate long term averages, and trends based on a number of input
files that typically cover a long time period. For example, if you want to analyse the modeled multi-year averages
of short wave fluxes for the North Sea you can do so by running the <em>createAveragesOnHexagon.py</em> script.
The script first calculates the averages on the supercomputer before downloading the resulting NetCDF4 files to
local computer. Finally, you can plot the results stored in the NetCDF4 file by running <em>createMaps.py</em>. The result
will look something like this:
![Shflux averaged 1993-2009](http://www.trondkristiansen.com/wp-content/gallery/romstools/longtermmean_shflux_time_depth_surface.jpg)

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
negative, and net transport through the transect. You can also caluclate the transport for the e.g. just the upper
500 meters of the water column by defining the minimum and maximum depths:

```Python
minDepth=0
maxDepth=500
```

<h4> Requirements </h4>
This toolbox still requires you to have a Fortran compiler to compile the Fortran programs and generate Python modules.
On Hexagon this is done by loading the gnu modules (for Fortran compiler compatible with Python and numpy). In the
terminal window type:

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

<h3> Contact </h3>

<ul>
<li>me (at) trondkristiansen.com</li>
</ul>




