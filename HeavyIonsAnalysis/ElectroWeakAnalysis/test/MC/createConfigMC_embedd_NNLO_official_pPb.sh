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
config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/Embedded/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
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

config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/Embedded/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
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

config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/Embedded/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
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
config.Data.inputDBS = '\''global'\''
config.Data.splitting = '\''FileBased'\''
config.Data.unitsPerJob = '$4'

config.Data.outLFNDirBase = '\''/store/user/%s/EWQAnalysis2017/MC/Embedded/TTree/%s'\'' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

config.section_('\''Site'\'')
config.Site.storageSite = '\''T2_FR_GRIF_LLR'\'

}

path=$PWD

step1DATE=20160202
step2DATE=20170205
step3DATE=20170206
step4DATE=20170813
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

numJobs=7
numTotJobs=14
declare -a workJobs=(
    "POWHEG_DYtoMuMu_M_10_30_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_DYtoMuMu_M_30_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_WToMuNu_Plus_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_WToMuNu_Minus_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_WToTauNu_Plus_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_WToTauNu_Minus_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_TTall_CT14_EPPS16_8160GeV_pythia8"
    "POWHEG_DYtoMuMu_M_10_30_CT14_EPPS16_8160GeV_pythia8_reverse"
    "POWHEG_DYtoMuMu_M_30_CT14_EPPS16_8160GeV_pythia8_reverse"
    "POWHEG_WToMuNu_Plus_CT14_EPPS16_8160GeV_pythia8_reverse"
    "POWHEG_WToMuNu_Minus_CT14_EPPS16_8160GeV_pythia8_reverse"
    "POWHEG_WToTauNu_Plus_CT14_EPPS16_8160GeV_pythia8_reverse"
    "POWHEG_WToTauNu_Minus_CT14_EPPS16_8160GeV_pythia8_reverse"
    "POWHEG_TTall_CT14_EPPS16_8160GeV_pythia8_reverse"
)
declare -a genOutputPDs=("DYtoMuMu" "WToMuNu" "WToTauNu" "QCD" "DYtoMuMu" "WToMuNu" "WToTauNu" "QCD")
declare -ai genNumCrabEvts='([0]="800" [1]="800" [2]="800" [3]="800" [4]="800" [5]="800" [6]="800" [7]="800" )'

declare -a rawInputPDs=(
)
declare -ai rawNumCrabEvts='([0]="2" [1]="2" [2]="2" [3]="2" [4]="2" [5]="2")'

declare -a recoInputPDs=(
)
declare -ai recoNumCrabEvts='([0]="2" [1]="2" [2]="2" [3]="2" [4]="2" [5]="2")'


declare -a anaInputPDs=(
    "/DYtoMuMu_M-10to30_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/DYtoMuMu_M-30_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/WpToMuNu_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/WmToMuNu_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/WpToTauNu_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/WmToTauNu_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/Ttbar_pPb-EmbEPOS_8p16_Powheg/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/DYtoMuMu_M-10to30_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v1/AODSIM"
    "/DYtoMuMu_M-30_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v1/AODSIM"
    "/WpToMuNu_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v2/AODSIM"
    "/WmToMuNu_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v1/AODSIM"
    "/WpToTauNu_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v1/AODSIM"
    "/WmToTauNu_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v1/AODSIM"
    "/Ttbar_PbP-EmbEPOS_8p16_Powheg/pPb816Summer16DR-PbPEmb_80X_mcRun2_pA_v4-v1/AODSIM"
)
declare -ai anaNumCrabEvts='([0]="1" [1]="1" [2]="1" [3]="1" [4]="1" [5]="1" [6]="1" [7]="1" [8]="1" [9]="1" [10]="1" [11]="1" [12]="1" [13]="1" )'

for (( i=4; i<5; i++ )); do
    for (( j=1; j<${numTotJobs}+1; j++ )); do

        logFileName=${stepPaths[$i-1]}/$logPath/step$i\_Embedded\_Official\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.log
        rootFileName=step$i\_Embedded\_Official\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.root
        crabFileName=${stepPaths[$i-1]}/$crabPath/step$i\_Embedded\_Official\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.py
        crabLogFileName=${stepPaths[$i-1]}/$crabPath/step$i\_Embedded\_Official\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.log
        cfgFileName=${stepPaths[$i-1]}/$pyPath/step$i\_Embedded\_Official\_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}.py
        if [ $i -gt 1 ]; then
            prevRootFileName=${stepPaths[$i-2]}/$rootPath/step$(expr $i - 1)\_Embedded\_Official\_${workJobs[$j-1]}_${stepTypes[$i-2]}_${stepDates[$i-2]}.root
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
                delete dir "${stepPaths[$i-1]}/$crabPath/crab_projects/crab_Embedded_Official_${workJobs[$j-1]}_${stepTypes[$i-1]}_${stepDates[$i-1]}/"
                echo 'Step 4: Producing '${workJobs[$j-1]}' '${stepTypes[$i-1]}' CRAB config file'
                crab_${stepTypes[$i-1]} Embedded\_Official\_${workJobs[$j-1]}\_${stepTypes[$i-1]}\_${stepDates[$i-1]} $cfgFileName ${anaInputPDs[$j-1]} ${anaNumCrabEvts[$j-1]} >> $crabFileName
                cd ${stepPaths[$i-1]}/$crabPath
                crab submit -c $crabFileName --dryrun >& $crabLogFileName
            fi
        fi

    done
done
    
