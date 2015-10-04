#!/bin/bash

thecfg=$1
outputbase=$2

# for file in /afs/cern.ch/user/c/chflores/work/public/TnP_2015/TP_Prod_Samples/MC/tnp_Prod_MC_PbPb_Cent*_08092015.root; do
# for file in /afs/cern.ch/user/c/chflores/work/public/TnP_2015/TP_Prod_Samples/Data/tnp_Prod_Data_PbPb_Cent*_25082015.root; do
for file in /afs/cern.ch/work/e/echapon/public/TnP_2015/TP_Prod_Samples/Data/tnp_Prod_Data_PbPb_Cent*_25082015.root; do
   inputfile=$file
   outputfile="${outputbase}`echo $file | grep -o '_Cent.*_'`.root"
   outputlog="${outputbase}`echo $file | grep -o '_Cent.*_'`.log"

   echo "cmsRun $thecfg inputFiles=$inputfile outputFile=$outputfile &> $outputlog"
   nice cmsRun $thecfg inputFiles=$inputfile outputFile=$outputfile &> $outputlog
done
