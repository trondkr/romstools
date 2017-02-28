#!/bin/bash
#
if [ $# -ne 6 ]; then
  echo
  echo "Usage: ncfile_subsetting.sh <YREF> <MREF> <DREF> <name-of-nc-file> <time-var> <t-unit>"
  echo
  echo " This script separates a large ROMS input file into several subsets in time; one file for each time step"
  echo " YREF MREF DREF refer to the reference date used in ROMS"
  echo " time-var is the name of the time variable"
  echo " t-unit is either s(econds) or d(ays)"
  echo
  echo " Example:"
  echo " ncfile_subsetting.sh 1948 01 01 norkyst_800m_his.nc ocean_time s"
  exit
fi
#
BASE=`pwd`
cd ${BASE}
#
# Name of his-file
yref=${1}
mref=${2}
dref=${3}
rfil=${4}
tvar=${5}
tunit=${6}
#
# Create temporary filenames
tmpfila=tmpfil.${RANDOM}
tmpfilb=tmpfil.${RANDOM}
tmpfilc=tmpfil.${RANDOM}
tmpfild=tmpfil.${RANDOM}
tmpfile=tmpfil.${RANDOM}
tmpfilf=tmpfil.${RANDOM}
tmpfilg=tmpfil.${RANDOM}
#
# Dump all time stamps (Julian days) to ascii
ncdump -v ${tvar} -l 12 ${rfil} > ${BASE}/${tmpfila}
# Remove all except time data values
i=1
while read a b; do
  case ${a} in
  "data:" ) break ;;
  *       ) i=`expr ${i} + 1` ;;
  esac
done < ${BASE}/${tmpfila}
i=`expr ${i} + 2`
wc -l ${BASE}/${tmpfila} > ${BASE}/${tmpfilc}; read nl f < ${BASE}/${tmpfilc}; rm ${BASE}/${tmpfilc}
n=`expr ${nl} - ${i}`
tail -${n} ${BASE}/${tmpfila} > ${BASE}/${tmpfilb}
rm ${BASE}/${tmpfila}
wc -l ${BASE}/${tmpfilb} > ${BASE}/${tmpfilc}; read nl f < ${BASE}/${tmpfilc}; rm ${BASE}/${tmpfilc}
n=`expr ${nl} - 1`
head -${n} ${BASE}/${tmpfilb} > ${BASE}/${tmpfila}
rm ${BASE}/${tmpfilb}
perl -pe "s#,##g" < ${BASE}/${tmpfila} > ${BASE}/${tmpfild}
rm ${BASE}/${tmpfila}
# We now have a file with Julian days.
i=0
while read jd dummy; do
  # Update counter (time steps)
  i=`expr ${i} + 1`
  if [ ${tunit} == "d" ]; then
    jds=`echo ${jd} \* 86400.0 | bc`
  else
    jds=${jd}
  fi
  echo "date --utc --date='${yref}-${mref}-${dref} 00:00:00 UTC +${jds} seconds' '+%Y %m %d %H %M %S'" > ${BASE}/${tmpfile}
  sh ${BASE}/${tmpfile} > ${BASE}/${tmpfilf}
  read y m d h mi s < ${BASE}/${tmpfilf}
  rm ${BASE}/${tmpfilf} ${BASE}/${tmpfile}
  # Retrieve time step i
  ncea -O -F -d ${tvar},${i},${i} ${rfil} ${rfil}_${y}${m}${d}${h}${mi}
  echo "Finished time step ${i} (${y}-${m}-${d} ${h}:${mi})"
done < ${BASE}/${tmpfild}
rm ${BASE}/${tmpfild}
#
exit
#
