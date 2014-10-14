#!/bin/bash
#
#--------------------------------------------------------
# Make file with river forcing dependent on vertical grid
#--------------------------------------------------------
#
# Model grid specifications
Lm=250; Mm=200; No=35
Vtransform=1; Vstretching=1; theta_s=8.0; theta_b=0.1; Tcline=10.0; dx=8000;
#
gridfile=../nordsjoen_8km_grid.nc
stationfile=RiverOutlets_Nordsjoen-8km.txt
riverfluxfile=Q_Nordsjoen-8km_daily.dat
#
cd River
# Start with editing cdl-file for nc-generation, then execute fortran-program to fill in data
perl -pe "s#NINTRHOPOINTS_I#${Lm}#g"      < frc_rivers_keyword.cdl > frc_rivers.cdl
perl -pe "s#NINTRHOPOINTS_J#${Mm}#g"      < frc_rivers.cdl > frc_rivers.tmp; mv frc_rivers.tmp frc_rivers.cdl
perl -pe "s#NUMBEROFSRHOSURFACES#${No}#g" < frc_rivers.cdl > frc_rivers.tmp; mv frc_rivers.tmp frc_rivers.cdl
#
# Find no. of rivers from reading stationfile
n_rivers=0
while read river Xpos Ypos direction signs; do
  n_rivers=`expr ${n_rivers} + 1`  # Update counter for each line
done < ${stationfile}
echo "No. of rivers that are included in the sub-domain are ${n_rivers}"
perl -pe "s#N_RIVERS#${n_rivers}#g"          < frc_rivers.cdl > frc_rivers.tmp; mv frc_rivers.tmp frc_rivers.cdl
#
# Find no. of time stamps from reading riverfluxfile
nd=0
while read yy_riv mm_riv dd_riv fluxes; do
  nd=`expr ${nd} + 1`              # Update counter for each line
done < ${riverfluxfile}
echo "Size of river time dimension is ${nd}"
perl -pe "s#N_RIVER_TIME#${nd}#g"            < frc_rivers.cdl > frc_rivers.tmp; mv frc_rivers.tmp frc_rivers.cdl
#
# Generate empty river.nc-file with cdl-information
ncgen -o river.nc frc_rivers.cdl
rm frc_rivers.cdl

# Run program to fill in data (daily climatology is calculated and filled in where runoff values are missing)
if [ -s MakeRivers_Nordsjoen ]; then rm MakeRivers_Nordsjoen; fi
make
if [ -s MakeRivers_Nordsjoen ]; then
  ./MakeRivers_Nordsjoen << EOF
${gridfile}
${stationfile}
${riverfluxfile}
${n_rivers} ${nd}
${No} ${Vtransform} ${Vstretching} ${theta_s} ${theta_b} ${Tcline} ${dx}
EOF
else
  echo "Program MakeRivers_Nordsjoen did not compile properly!"
  exit
fi
#
echo
mv river.nc ../nordsjoen_8km_river.nc
#
cd ..
#
exit
#

