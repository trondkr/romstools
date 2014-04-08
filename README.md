romstools
=========

A variety of tools for plotting and analysing output from the Regional Ocean Model (ROMS)


<h2> Volume flux calculations </h2>

To calculate the fluxes you run the script <strong>calculateFluxes.py</strong> either as a standalone python program or executed
as a job-script. If you run the program on a super-computer (e.g. Hexagon) you may want to run the script as a job
and use the script <strong>runJob.sh</strong> to queue your job (Hexagon: qsub runJob.sh).

This toolbox is a mixture of Fortrand and Python tools where the core programs are taken from the excellent pyroms
toolbox. However, to avoid having to install the entire toolbox to calculate the fluxes, this smaller toolbox was created.

This toolbox still requires you to have a Fortran compiler to compile the Fortran programs and generate Python modules.
On Hexagon this is done by loading the gnu modules (for Fortran compiler compatible with Python and numpy). In the
terminal window type:

<code>
module swap PrgEnv-pgi PrgEnv-gnu
module unload notur
f2py --verbose  -c -m iso iso.f90
f2py --verbose  -c -m obs_interp obs_interp.f
</code>

This should provide you with two python modules (obs_interp.so and iso.so) which you can try to import to python with:

<code>
Python 2.7.2 (default, Mar 22 2012, 12:32:11)
[GCC 4.6.1 20110627 (Cray Inc.)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import iso
>>> print iso.__doc__
</code>

You can define a list of transects you want the volume transport calculated for in the function
<strong>defineTransects().</strong>.  The output from running is a comma separated value file containing the positive,
negative, and net transport through the transect. You can also caluclate the transport for the e.g. just the upper
500 meters of the water column by defining the minimum and maximum depths:

<code>
minDepth=0
maxDepth=500
</code>



