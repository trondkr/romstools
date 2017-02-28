#!/bin/bash
#
#  Give the job a name
#PBS -N "WS4KMmakeAtmosphericForcing"
#
#  Specify the project the job belongs to
#PBS -A nn9297k
#
#PBS -q normal
#
#PBS -l mppwidth=1
#
# Request 1 hour of walltime
#PBS -l walltime=02:00:00
#
#  Send me an email on  a=abort, b=begin, e=end
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@niva.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o  makeAtmos.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e  makeAtmos.err
#

#  Make sure I am in the correct directory
cd /work/shared/nn9297k/FAABolous/atmos_forcing
export MPLCONFIGDIR=${pwd}

export TMP=`pwd`
module unload notur
aprun -B ./make_atmos.sh > output.log
