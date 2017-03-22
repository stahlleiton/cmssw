import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras

process = cms.Process('PAT',eras.Run2_2016_pA)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "HiEWQForest.root"
options.inputFiles =  "file:0E6009E4-65B0-E611-9C6B-02163E014503.root"
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_v19', '') #80X_dataRun2_Prompt_v15
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")


process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(
    record = cms.string('EcalLaserAlphasRcd'),
    tag = cms.string("EcalLaserAlphas_v2_prompt"),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
    ),
  cms.PSet(
    record = cms.string('EcalLaserAPDPNRatiosRcd'),
    tag = cms.string("EcalLaserAPDPNRatios_prompt_v2"),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
    )
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)#,
   # wantSummary = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)

process.load('HeavyIonsAnalysis.ElectroWeakAnalysis.metAnalyzer_cfi')
process.load('HeavyIonsAnalysis.ElectroWeakAnalysis.muonAnalyzer_cfi')
#process.load('HeavyIonsAnalysis.ElectroWeakAnalysis.pfcandAnalyzer_pA_cfi')


# Trigger Filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltEWQHI = hltHighLevel.clone(
    HLTPaths =  ['HLT_PAL3Mu12_v1'] , 
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT") , 
    throw = False , andOr = True , eventSetupPathsKey = cms.string( '' ) 
    )
process.triggerSelection = cms.Sequence( process.hltEWQHI )

process.selectMuonEvts = cms.EDFilter("PATMuonSelector", src = cms.InputTag( 'patMuons' ) , cut = cms.string("isPFMuon() && pt()>15.") )
process.selectMuonFilter = cms.EDFilter("PATCandViewCountFilter", minNumber = cms.uint32(1), maxNumber = cms.uint32(9999999), src = cms.InputTag("selectMuonEvts") )
process.muonSelection = cms.Sequence( process.selectMuonEvts + process.selectMuonFilter )

process.metAnaNoHF = process.metAna.clone(
    patMETTag      = cms.InputTag("slimmedMETsNoHF"),
    pfMETTag       = cms.InputTag("pfMetNoHF"),
    caloMETTag     = cms.InputTag(""),
    )
process.metAnaSeqNoHF = cms.Sequence( process.metAnaNoHF )

process.anaMET  = cms.EndPath( process.metAnaSeq * process.metAnaSeqNoHF )
#process.anaPath = cms.Path(process.triggerSelection + process.muonSelection + process.pfcandAnalyzer )
process.anaPath = cms.Path( process.muonAnaSeq )

#Options:
process.source    = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.inputFiles ) )
process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# Schedule definition
process.schedule = cms.Schedule(process.anaPath,process.anaMET,process.endjob_step)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
process.load('Configuration.StandardSequences.PAT_cff')
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from HeavyIonsAnalysis.ElectroWeakAnalysis.PatAlgos.slimming.miniAOD_tools import miniAOD_ForHiEWQ_customizeAllData

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_ForHiEWQ_customizeAllData(process)
