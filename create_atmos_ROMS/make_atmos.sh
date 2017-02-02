#!/bin/bash
#
#--------------------------------------------------------
# Make file(s) with atmospheric forcing
#--------------------------------------------------------
#

#
gridfile=/work/shared/nn9297k/FAABolous/WS4KM_grd.nc
atmeraintdir=/work/shared/norkyst/NorKyst-800m_Forcing/Atmos/ERA_Interim
#
cd Atmos
if [ -s get_forcing ]; then rm get_forcing *.o; fi
make
if [ -s get_forcing ]; then

for ((year=1980; year<=2015; year++)); 

	do  echo "Creating ERA INTERIM forcing file for year: $year"
    ./get_forcing << EOF
${gridfile}
${atmeraintdir}/${year}0101AN.nc
${atmeraintdir}/${year}0101FC.nc
${atmeraintdir}/ERA_fields_mask.nc
./
${year}
EOF
    for param in Uwind Vwind Pair Tair Qair cloud rain swrad lwrad; do
      mv ${param}.nc ../ws4km_${param}_${year}.nc
    done
  done  # years
else
  echo "Program Atmos/get_forcing did not compile properly!"
  exit
fi
#
cd ..
#
exit
#

