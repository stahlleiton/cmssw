#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: ./createConfigMC.sh <nEvt>"
  exit 1
fi

delete () {
    if [[ $1 == "dir" ]]; then
        if [ -f $1 ]; then
            rm -f -r $2
        fi
    fi
    if [[ $1 == "file" ]]; then
        if [ -f $1 ]; then
            rm -f $2
        fi
    fi
}

nEvt=$1
crab_GEN () {
    
    echo 'from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('\''General'\'')
config.General.requestName = '\'$1\''
config.General.workArea = '\''crab_projects'\''
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('\''JobType'\'')
config.JobType.pluginName = '\''PrivateMC'\''
config.JobType.psetName = '\'$2\''
config.JobType.maxMemoryMB = 2500

config.section_('\''Data'\'')
config.Data.outputPrimaryDataset = '\'$3\''
config.Data.splitting = '\''EventBased'\''
config.Data.unitsPerJob = '$4'
NJOBS = '$5'
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName

config.section_('\''Site'\'')
config.Site.storageSite = '\''T2_FR_GRIF_LLR'\'

}

crab_RAW () {

    echo 'from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('\''General'\'')
config.General.requestName = '\'$1\''
config.General.workArea = '\''crab_projects'\''
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('\''JobType'\'')
config.JobType.pluginName = '\''Analysis'\''
config.JobType.psetName = '\'$2\''
config.JobType.maxMemoryMB = 2500

config.section_('\''Data'\'')
config.Data.inputDataset = '\'$3\''
config.Data.inputDBS = '\''phys03'\''
config.Data.splitting = '\''FileBased'\''
config.Data.unitsPerJob = '$4'

config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName

config.section_('\''Site'\'')
config.Data.ignoreLocality = True
config.Site.whitelist = ['\''T2_FR_GRIF_*'\'']
config.Site.storageSite = '\''T2_FR_GRIF_LLR'\'

}

crab_RECO () {

    echo 'from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('\''General'\'')
config.General.requestName = '\'$1\''
config.General.workArea = '\''crab_projects'\''
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('\''JobType'\'')
config.JobType.pluginName = '\''Analysis'\''
config.JobType.psetName = '\'$2\''
config.JobType.maxMemoryMB = 2500

config.section_('\''Data'\'')
config.Data.inputDataset = '\'$3\''
config.Data.inputDBS = '\''phys03'\''
config.Data.splitting = '\''FileBased'\''
config.Data.unitsPerJob = '$4'

config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName

config.section_('\''Site'\'')
config.Data.ignoreLocality = True
config.Site.whitelist = ['\''T2_FR_GRIF_*'\'']
config.Site.storageSite = '\''T2_FR_GRIF_LLR'\'

}

crab_ANA () {

    echo 'from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('\''General'\'')
config.General.requestName = '\'$1\''
config.General.workArea = '\''crab_projects'\''
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('\''JobType'\'')
config.JobType.pluginName = '\''Analysis'\''
config.JobType.psetName = '\'$2\''
config.JobType.maxMemoryMB = 2500

config.section_('\''Data'\'')
config.Data.inputDataset ='\'$3\''
config.Data.inputDBS = '\''phys03'\''
config.Data.splitting = '\''FileBased'\''
config.Data.unitsPerJob = '$4'

config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/TTree/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

config.section_('\''Site'\'')
config.Data.ignoreLocality = True
config.Site.whitelist = ['\''T2_FR_GRIF_*'\'']
config.Site.storageSite = '\''T2_FR_GRIF_LLR'\'

}

path=$PWD

step1DATE=20160202
step2DATE=20170205
step3DATE=20170206
step4DATE=20170323
declare -a stepDates=("$step1DATE" "$step2DATE" "$step3DATE" "$step4DATE")

step1Path=$path/GEN
step2Path=$path/RAW
step3Path=$path/RECO
step4Path=$path/ANA
declare -a stepPaths=("$step1Path" "$step2Path" "$step3Path" "$step4Path")
declare -a stepTypes=("GEN" "RAW" "RECO" "ANA")

pyPath=python
rootPath=output
logPath=log
crabPath=crab
declare -a workPaths=("$pyPath" "$rootPath" "$logPath" "$crabPath")

for stepPath in "${stepPaths[@]}"; do 
    for workPath in "${workPaths[@]}"; do 
        mkdir -p $stepPath/$workPath;
    done;
done;

numJobs=4
numTotJobs=8
declare -a workJobs=(
    "Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6"
    "Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6"
    "Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola"
    "Pythia_QCD_PtHat20_MuEnrichedPt15_8160GeV_pythia8"
    "Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_reverse"
    "Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_reverse"
    "Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_reverse"
    "Pythia_QCD_PtHat20_MuEnrichedPt15_8160GeV_pythia8_reverse"
)
declare -a genOutputPDs=("DYtoMuMu" "WToMuNu" "WToTauNu" "QCD" "DYtoMuMu" "WToMuNu" "WToTauNu" "QCD")
declare -ai genNumCrabEvts='([0]="800" [1]="800" [2]="800" [3]="800" [4]="800" [5]="800" [6]="800" [7]="800" )'

declare -a rawInputPDs=(
    "/DYtoMuMu/anstahll-Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_GEN_20160202-df304e47766b2b2d38e4f1efb5cf96f9/USER"
    "/WToMuNu/anstahll-Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_GEN_20160202-17e7429fb944ce1534b03c1aa1a8be70/USER"
    "/WToTauNu/anstahll-Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_GEN_20160202-1756a8c2280a14fc5ff38f04338bfb48/USER"
    "/DYtoMuMu/anstahll-Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_reverse_GEN_20160202-232d1ce10bd15274bed0faf4e6e9f45b/USER"
    "/WToMuNu/anstahll-Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_reverse_GEN_20160202-d6fa70ff3d3710f4498eacd4411fcffa/USER"
    "/WToTauNu/anstahll-Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_reverse_GEN_20160202-30b6803042a18773b865d31a6be1313e/USER"
)
declare -ai rawNumCrabEvts='([0]="2" [1]="2" [2]="2" [3]="2" [4]="2" [5]="2")'

declare -a recoInputPDs=(
    "/DYtoMuMu/anstahll-Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_RAW_20170205-a29fe75007ca498b0c4693f898031040/USER"
    "/WToMuNu/anstahll-Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_RAW_20170205-f8522175fd15740b7e91ce3907c2d7c0/USER"
    "/WToTauNu/anstahll-Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_RAW_20170205-a29fe75007ca498b0c4693f898031040/USER"
    "/DYtoMuMu/anstahll-Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_reverse_RAW_20170205-f8522175fd15740b7e91ce3907c2d7c0/USER"
    "/WToMuNu/anstahll-Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_reverse_RAW_20170205-f8522175fd15740b7e91ce3907c2d7c0/USER"
    "/WToTauNu/anstahll-Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_reverse_RAW_20170205-a29fe75007ca498b0c4693f898031040/USER"
)
declare -ai recoNumCrabEvts='([0]="2" [1]="2" [2]="2" [3]="2" [4]="2" [5]="2")'


declare -a anaInputPDs=(
    "/DYtoMuMu/anstahll-Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_RECO_20170206-af31df56e13449bd40098ffa7c7e1173/USER"
    "/WToMuNu/anstahll-Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_RECO_20170206-af31df56e13449bd40098ffa7c7e1173/USER"
    "/WToTauNu/anstahll-Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_RECO_20170206-af31df56e13449bd40098ffa7c7e1173/USER"
    "/QCD/miheejo-PtHat20_MuEnrichedPt15_pPb_RECO_v1-af31df56e13449bd40098ffa7c7e1173/USER"
    "/DYtoMuMu/anstahll-Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_reverse_RECO_20170206-af31df56e13449bd40098ffa7c7e1173/USER"
    "/WToMuNu/anstahll-Pyquen_WToMuNu_TuneZ2_8160GeV_pythia6_reverse_RECO_20170206-af31df56e13449bd40098ffa7c7e1173/USER"
    "/WToTauNu/anstahll-Pyquen_WToTauNu_TuneZ2_8160GeV_pythia6_tauola_reverse_RECO_20170206-af31df56e13449bd40098ffa7c7e1173/USER"
    "/QCD/miheejo-PtHat20_MuEnrichedPt15_Pbp_reverse_RECO_v1-af31df56e13449bd40098ffa7c7e1173/USER"
)
declare -ai anaNumCrabEvts='([0]="2" [1]="2" [2]="2" [3]="14" [4]="2" [5]="2" [6]="2" [7]="14")'

for (( i=4; i<5; i++ )); do
    for (( j=1; j<${numTotJobs}+1; j++ )); do

        logFileName=${stepPaths[$i-1]}/$logPath/step$i\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.log
        rootFileName=step$i\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.root
        crabFileName=${stepPaths[$i-1]}/$crabPath/step$i\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.py
        crabLogFileName=${stepPaths[$i-1]}/$crabPath/step$i\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.log
        cfgFileName=${stepPaths[$i-1]}/$pyPath/step$i\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.py
        if [ $i -gt 1 ]; then
            prevRootFileName=${stepPaths[$i-2]}/$rootPath/step$(expr $i - 1)\_${workJobs[$j-1]}_${stepTypes[$i-2]}_${stepDates[$i-2]}.root
        fi
        
        if [ $i -eq 1 ]; then
            if [ ! -f $cfgFileName ]; then
                echo 'Step 1: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' cmsDriver config file'
                cd ${stepPaths[$i-1]}/$rootPath/
                if [ $j -le $numJobs ]; then
                    start=$SECONDS
                    cmsDriver.py ${workJobs[$j-1]}_cfi.py --python_filename $cfgFileName --mc --eventcontent RAWSIM --customise Configuration/StandardSequences/SimWithCastor_cff.py --datatier GEN-SIM --conditions 80X_mcRun2_pA_v4 --beamspot RealisticPPbBoost8TeV2016Collision --step GEN,SIM --scenario HeavyIons --era Run2_2016_pA --number=$nEvt --fileout=$rootFileName >& $logFileName
                    duration=$(( SECONDS - start ))
                    echo 'seconds: '$duration >> $logFileName
                    ls -lh $rootFileName >> $logFileName
                else 
                    start=$SECONDS
                    cmsDriver.py ${workJobs[$j-1]}_cfi.py --python_filename $cfgFileName --mc --eventcontent RAWSIM --customise Configuration/StandardSequences/SimWithCastor_cff.py --datatier GEN-SIM --conditions 80X_mcRun2_pA_v4 --beamspot RealisticPbPBoost8TeV2016Collision --step GEN,SIM --scenario HeavyIons --era Run2_2016_pA --number=$nEvt --fileout=$rootFileName >& $logFileName
                    duration=$(( SECONDS - start ))
                    echo 'seconds: '$duration >> $logFileName
                    ls -lh $rootFileName >> $logFileName
                fi
            fi
            if [ ! -f $crabFileName ]; then
                delete dir "${stepPaths[$i-1]}/$crabPath/crab_projects/crab_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}/"
                echo 'Step 1: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' CRAB config file'
                crab_${stepTypes[$i-1]} ${workJobs[$j-1]}\_${stepTypes[$i-1]}\_${stepDates[$i-1]} $cfgFileName ${genOutputPDs[$j-1]} ${genNumCrabEvts[$j-1]} 2500 >> $crabFileName
                cd ${stepPaths[$i-1]}/$crabPath
                crab submit -c $crabFileName --dryrun >& $crabLogFileName
            fi
        fi

        if [ $i -eq 2 ]; then
            if [ ! -f $cfgFileName ]; then
                echo 'Step 2: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' cmsDriver config file'
                cd ${stepPaths[$i-1]}/$rootPath/
                start=$SECONDS
                cmsDriver.py --python_filename $cfgFileName --filein=file:$prevRootFileName  --fileout=$rootFileName --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions 80X_mcRun2_pA_v4 --step DIGI,L1,DIGI2RAW,HLT:PIon --era Run2_2016_pA  --customise Configuration/DataProcessing/Utils.addMonitoring --number=$nEvt >& $logFileName
                duration=$(( SECONDS - start ))
                echo 'seconds: '$duration >> $logFileName
                ls -lh $rootFileName >> $logFileName
            fi
            if [ ! -f $crabFileName ]; then
                delete dir "${stepPaths[$i-1]}/$crabPath/crab_projects/crab_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}/"
                echo 'Step 2: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' CRAB config file'
                crab_${stepTypes[$i-1]} ${workJobs[$j-1]}\_${stepTypes[$i-1]}\_${stepDates[$i-1]} $cfgFileName ${rawInputPDs[$j-1]} ${rawNumCrabEvts[$j-1]} >> $crabFileName
                cd ${stepPaths[$i-1]}/$crabPath
                crab submit -c $crabFileName --dryrun >& $crabLogFileName
            fi
        fi


        if [ $i -eq 3 ]; then
            if [ ! -f $cfgFileName ]; then
                echo 'Step 3: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' cmsDriver config file'
                cd ${stepPaths[$i-1]}/$rootPath/
                start=$SECONDS
                cmsDriver.py --python_filename $cfgFileName --filein=file:$prevRootFileName  --fileout=$rootFileName --mc --eventcontent AODSIM --datatier GEN-SIM-RECO --conditions 80X_mcRun2_pA_v4 --step RAW2DIGI,L1Reco,RECO --era Run2_2016_pA  --customise_commands "process.bunchSpacingProducer.bunchSpacingOverride=cms.uint32(25)\nprocess.bunchSpacingProducer.overrideBunchSpacing=cms.bool(True)" --number=$nEvt >& $logFileName
                duration=$(( SECONDS - start ))
                echo 'seconds: '$duration >> $logFileName
                ls -lh $rootFileName >> $logFileName
            fi
            if [ ! -f $crabFileName ]; then
                delete dir "${stepPaths[$i-1]}/$crabPath/crab_projects/crab_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}/"
                echo 'Step 3: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' CRAB config file'
                crab_${stepTypes[$i-1]} ${workJobs[$j-1]}\_${stepTypes[$i-1]}\_${stepDates[$i-1]} $cfgFileName ${recoInputPDs[$j-1]} ${recoNumCrabEvts[$j-1]} >> $crabFileName
                cd ${stepPaths[$i-1]}/$crabPath
                crab submit -c $crabFileName --dryrun >& $crabLogFileName
            fi
        fi

        if [ $i -eq 4 ]; then
            cfgFileName=$path/lo.py
            if [ ! -f $crabFileName ]; then
                delete dir "${stepPaths[$i-1]}/$crabPath/crab_projects/crab_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}/"
                echo 'Step 4: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' CRAB config file'
                crab_${stepTypes[$i-1]} ${workJobs[$j-1]}\_${stepTypes[$i-1]}\_${stepDates[$i-1]} $cfgFileName ${anaInputPDs[$j-1]} ${anaNumCrabEvts[$j-1]} >> $crabFileName
                cd ${stepPaths[$i-1]}/$crabPath
                crab submit -c $crabFileName --dryrun >& $crabLogFileName
            fi
        fi

    done
done
    
