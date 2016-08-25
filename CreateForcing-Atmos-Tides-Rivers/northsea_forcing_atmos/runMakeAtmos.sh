#!/bin/bash
#
#  Give the job a name
#PBS -N "makeAtmosphericForcing"
#
#  Specify the project the job belongs to
#PBS -A imr
#
#PBS -q normal
#
#PBS -l mppwidth=1
#
# Request 1 hour of walltime
#PBS -l walltime=1:00:00
#
#  Send me an email on  a=abort, b=begin, e=end
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@imr.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o  makeAtmos.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e  makeAtmos.err
#

#  Make sure I am in the correct directory
cd /work/shared/imr/NS8KM/romstools/CreateForcing-Atmos-Tides-Rivers/northsea_forcing_atmos
export MPLCONFIGDIR=${pwd}

export TMP=`pwd`
module unload notur
aprun -B ./make_atmos.sh > makeAtmos_output_19072016.log
