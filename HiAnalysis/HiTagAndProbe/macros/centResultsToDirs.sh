#!/bin/bash

theFiles="MuTrk_cent_20150909_Cent*_.root"
theCfg="tp_Ana_Trk_RD_isGlb_1bin.py"

for file in $theFiles; do
   tag="res`echo $file | grep -o '_Cent.*_'`"
   mkdir $tag
   root -l -b -q plot.C+'("'${file}'")'
   mv Passing* $file `echo $file | sed 's/.root/.log/g'` $tag
   cp $theCfg $tag
done
