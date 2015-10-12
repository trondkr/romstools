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
#PBS -A imr
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
#PBS -M trond.kristiansen@imr.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o ./TIDES1600M.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e ./TIDES1600M.err
#
#  Make sure I am in the correct directory
cd /work/users/trondk/KINO/romstools/CreateForcing-Atmos-Tides-Rivers/northsea_forcing_tides

# For hexagon use:
ulimit -c unlimited
aprun -B ./make_tides.sh > tides.log

# Return output at end to mpiexec as exit status:
exit $?