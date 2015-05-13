#!/bin/bash

eval `scramv1 runtime -sh`

if [ $# -ne 2 ]; then
echo "Usage: $0 [# crab jobs] [Input directory]  "
  exit;
fi


prefix=$2
echo $1 "jobs on" $prefix
#/store/user/miheejo/HiOniaAnalyzer/2011PbPb/jpsiMuMu_JpsiPt1530/

########## File 'list' should contains only .root file names WITHOUT any prefix directories.
########## Ex) (O) Jpsi_Histos_11.root,
##########     (X) /store/user/miheejo/Jpsi_Histos_11.root
#eos ls $prefix > list
#cmsLs $prefix | awk '{print $5}' | awk 'BEGIN{FS="/"}; {print $8}' > list

########## Check duplications and missing files
i=1
missings=""
while [ $i -le $1 ] #Put number of total crab jobs
do
  checkfile="Jpsi_Histos_MC_"$i"_.*root"
  result=$(grep -c "$checkfile" list)
  if [ $result -ge 2 ]; then
    files=$(grep "$checkfile" list)
    echo $(echo $files | awk -v prefix=$prefix '{for(i=2; i<=NF; ++i) print "cmsRm "prefix$i}')
  elif [ $result -eq 0 ]; then
    missings=$missings"\n"$i
  fi

  remainder=1
  let "remainder = $i % 100"
  if [ $remainder -eq 0 ]; then
    echo "Processing: " $i
  fi

  i=$(($i+1))
done

echo "Missing files"
echo -e "$missings"
