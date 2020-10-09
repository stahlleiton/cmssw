# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: openHLT --python_file openHLT.py --data --era Run2_2016_pA --geometry=Extended2016,Extended2016Reco --process reHLT --conditions=auto:run2_hlt --step L1REPACK:Full,HLT:PIon --no_exec -n 5
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('reHLT',eras.Run2_2018_pp_on_AA)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('root://xrootd.cmsaf.mit.edu//store/hidata/HIRun2018A/HIZeroBiasReducedFormat/RAW/v1/000/327/466/00000/EE35AF29-16D2-F346-B579-F467BCDE2E34.root')
)

process.options = cms.untracked.PSet(

)
process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName =cms.string("openHLT.root"))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('openHLT nevts:5'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)
# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '112X_dataRun2_v4', '')

process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.RunParameters.isData = cms.untracked.bool(True)
process.hltbitanalysis.HLTProcessName = cms.string('HLT')
process.hltbitanalysis.hltresults = cms.InputTag( 'TriggerResults','','HLT' )
process.hltbitanalysis.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis")
process.hltbitanalysis.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis")
process.hltbitanalysis.gObjectMapRecord  = cms.InputTag("gtStage2ObjectMap")
process.hltbitanalysis.gmtStage2Digis    = cms.string("gmtStage2Digis")
process.hltbitanalysis.caloStage2Digis   = cms.string("caloStage2Digis")
process.hltbitanalysis.UseTFileService = cms.untracked.bool(True)
process.hltbitanalysis.getL1InfoFromEventSetup = cms.untracked.bool(False)
process.hltBitAnalysis = cms.EndPath(process.hltbitanalysis)
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("openHLT.root"))
