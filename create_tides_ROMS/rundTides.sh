#! /bin/sh -
#
#  Make sure I use the correct shell.
#
#PBS -S /bin/sh
#
#  Give the job a name
#
#PBS -N "KINO1600M-TIDES"
#
# Specify the project the job should be accounted on
#PBS -A nn9297k
#
#  We want 2 hours on 32 cpu's:
#
#PBS -l walltime=2:00:00
#PBS -l mppwidth=1
#PBS -l mppmem=2000MB
#PBS -l mppnppn=16
#  Send me an email on  a=abort, b=begin, e=end
#
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@niva.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o ./TIDES1600M.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e ./TIDES1600M.err
#
#  Make sure I am in the correct directory
cd /work/shared/nn9297k/A20/romstools/romstools/CreateForcing-Atmos-Tides-Rivers/northsea_forcing_tides
# For hexagon use:

module swap PrgEnv-pgi PrgEnv-gnu
module unload notur

#module unload xtpe-interlagos
aprun -B ./make_tides_step1.sh > tides.log

# Return output at end to mpiexec as exit status:
exit $?