#!/bin/sh

#PBS -N "FLUXES"
#PBS -A imr
#PBS -l mppwidth=1
#PBS -l walltime=90:00:00
#PBS -l mppmem=32000mb
#PBS -m abe
#PBS -M youremail@domain.com

# This script runs the volume fux calculations for MyOcean for a variety of transects
#
# Trond Kristiansen, 24.03.2014 for MYOCEAN2


module swap PrgEnv-cray PrgEnv-pgi
module load python

export drc_out=(
'/work/users/trondk/NS8km/FLUX/'
)

export MPLCONFIGDIR=$drc_out
cd $drc_out

aprun -B python ${drc_out}/calculateFluxes.py
