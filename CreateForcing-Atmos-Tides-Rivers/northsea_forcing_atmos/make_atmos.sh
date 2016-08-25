#!/bin/bash
#
#--------------------------------------------------------
# Make file(s) with atmospheric forcing
#--------------------------------------------------------
#
# List years (ERA-int available from 1989-today)
years="2014"
#
gridfile=/work/shared/imr/NS8KM/FORCING/GRID/nordsjoen_8km_grid_hmax20m_v4.nc

atmeraintdir=/work/shared/norkyst/NorKyst-800m_Forcing/Atmos/ERA_Interim
#
cd Atmos
if [ -s get_forcing ]; then rm get_forcing *.o; fi
make
if [ -s get_forcing ]; then
  for year in ${years}; do
    ./get_forcing << EOF
${gridfile}
${atmeraintdir}/${year}0101AN.nc
${atmeraintdir}/${year}0101FC.nc
${atmeraintdir}/ERA_fields_mask.nc
./
EOF
    for param in Uwind Vwind Pair Tair Qair cloud rain swrad lwrad; do
      mv ${param}.nc ../ns8km_${param}_${year}.nc
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

