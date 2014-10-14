#!/bin/bash
#
#--------------------------------------
# Make file with tidal forcing for ROMS
#--------------------------------------
#
# Start and end period for simulation
syear=2009; eyear=2012
# Input files
gridfile=../nordsjoen_8km_grid.nc
###tpxodir=/work/shared/norkyst/NorKyst-800m_Forcing/Tpxo7.2
tpxodir=/work/shared/norkyst/NorKyst-800m_Forcing/Tpxo_AtlanticOcean
#
cd Tpxo
if [ -s tpxo2grid ]; then rm tpxo2grid *.o; fi
make
./tpxo2grid ${gridfile} nc nordsjoen_8km_tpxo.nc ${tpxodir}   # Extract TPXO-data to new model domain
#
cd ../Tides
module load pgi/11.10.0  # Program core-dumps with a newer pgi-compiler
if [ -s tidenc2roms ]; then rm tidenc2roms *.o; fi
make
# Substitute start simulation date and end date in sup-file, and run program
cp ./tidenc2roms.sup_ORIG tmp1.fil
perl -pe "s#YY_START#${syear}#g" < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MM_START#01#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#DD_START#01#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#HH_START#00#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MI_START#00#g"       < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#YY_END#${syear}#g"   < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MM_END#12#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#DD_END#31#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#HH_END#00#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
perl -pe "s#MI_END#00#g"         < tmp1.fil > tmp2.fil; mv tmp2.fil tmp1.fil
mv tmp1.fil ./tidenc2roms.sup
ln -sf ../Tpxo/nordsjoen_8km_tpxo.nc ./tpxo.nc
./tidenc2roms
rm ./tpxo.nc
mv ./tide.nc ../nordsjoen_8km_tides.nc
module unload pgi/11.10.0
#
cd ..
#
exit
#
