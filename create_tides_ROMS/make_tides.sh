
#!/bin/bash
#
#--------------------------------------
# Make file with tidal forcing for ROMS
#--------------------------------------
#

# Start and end period for simulation
syear=1980; eyear=2017
# Input files
gridfile=/Users/trondkr/Dropbox/NIVA/ROHO800/Grid/ROHO800_grid.nc

# Name of app
APPNAME=ROHO800

#
cd Tpxo
if [ -s tpxo2grid ]; then rm tpxo2grid *.o; fi
make

# I had to run this command interactively to make this work (this script does not run as expected):
./tpxo2grid ${gridfile} nc ${APPNAME}_tpxo.nc   # Extract TPXO-data to new model domain
cd ..
#
export CRAY_CPU_TARGET=x86-64
cd Tides
#module load pgi/11.10.0  # Program core-dumps with a newer pgi-compiler
if [ -s tidenc2roms ]; then rm tidenc2roms *.o; fi
make
# Substitute start simulation date and end date in sup-file, and run program
cp ./tidenc2roms.sup_ORIG tmp1.fil
perl -pe "s#YY_START#${syear}#g" < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MM_START#01#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#DD_START#01#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#HH_START#00#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MI_START#00#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#YY_END#${eyear}#g"   < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MM_END#12#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#DD_END#31#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#HH_END#00#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MI_END#00#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
mv tmp1.fil ./tidenc2roms.sup
ln -sf ../Tpxo/${APPNAME}_tpxo.nc ./tpxo.nc
./tidenc2roms
rm ./tpxo.nc
mv ./tide.nc ../${APPNAME}_tides.nc
#module unload pgi/11.10.0
#
cd ..
#
exit
#
