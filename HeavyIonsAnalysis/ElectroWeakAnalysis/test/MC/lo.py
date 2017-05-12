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
options.inputFiles =  "file:step3_Pyquen_DYtoMuMu_M_30_TuneZ2_8TeV16_pythia6_RECO_20170429_1.root"
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_pA_v4', '')

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)#,
    #wantSummary = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)

process.load('HeavyIonsAnalysis.ElectroWeakAnalysis.metAnalyzer_cfi')
process.load('HeavyIonsAnalysis.ElectroWeakAnalysis.muonAnalyzer_cfi')

process.metAnaNoHF = process.metAna.clone(
    patMETTag      = cms.InputTag("slimmedMETsNoHF"),
    pfMETTag       = cms.InputTag("pfMetNoHF"),
    caloMETTag     = cms.InputTag(""),
    )
process.metAnaSeqNoHF = cms.Sequence( process.metAnaNoHF )

process.anaMET  = cms.EndPath( process.metAnaSeq * process.metAnaSeqNoHF )
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
from HeavyIonsAnalysis.ElectroWeakAnalysis.PatAlgos.slimming.miniAOD_tools import miniAOD_ForHiEWQ_customizeAllMC

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_ForHiEWQ_customizeAllMC(process)
