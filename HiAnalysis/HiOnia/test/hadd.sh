#!/bin/sh -f

indir="/store/user/miheejo/HiOniaAnalyzer/2011PbPb/jpsiMuMu_JpsiPt1530/"
fileN="Jpsi_Histos"
final="jpsiMuMu_JpsiPt1530_Histos_cmssw445p5_RegIt_hStats.root"

#indir="/store/user/miheejo/MCsample_pythiagun_445patch1/promptReco/HiOniaAnalyzer/"
#fileN="promptReco_Tree_"
#final="promptReco_Tree.root"

#indir="/store/user/miheejo/MCsample_pythiagun_445patch1/regIt/HiOniaAnalyzer/"
#fileN="RegIt_Tree_"
#final="RegIt_Tree.root"

outdir=$(pwd)

cmsLs $indir | grep $fileN > list
list2=$(awk -v p=$indir '{print "root://eoscms//eos/cms"$5}' < list)
echo $list2

hadd /tmp/miheejo/$final $list2
#cp /tmp/miheejo/$final.root $outdir

rm -f list
