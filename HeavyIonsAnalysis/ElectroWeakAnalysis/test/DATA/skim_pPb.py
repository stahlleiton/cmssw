import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SKIM',eras.Run2_2016_pA)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

# setup 'analysis'  options
options = VarParsing.VarParsing("Analysis")

# Input and Output File Names
options.outputFile = "HiEWQ_pPb8TeV_AOD_SKIM.root"
options.inputFiles =  "file:0E6009E4-65B0-E611-9C6B-02163E014503.root"
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v15', '')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)

# Trigger Filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltEWQHI = hltHighLevel.clone(
    HLTPaths =  ['HLT_PAL3Mu12_v1'] , 
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT") , 
    throw = False , andOr = True , eventSetupPathsKey = cms.string( '' ) 
    )
process.triggerSelection = cms.Sequence( process.hltEWQHI )

# Muon Filter
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons
process.patMuonsColl = patMuons.clone( embedCaloMETMuonCorrs = False , addGenMatch = False )
process.selectMuons = cms.EDFilter(
    "PATMuonSelector", 
    src = cms.InputTag( 'patMuonsColl' ) , 
    cut = cms.string("isPFMuon() && pt()>15.") 
    )
process.countMuons = cms.EDFilter(
    "PATCandViewCountFilter", 
    minNumber = cms.uint32(1), 
    maxNumber = cms.uint32(999999999), 
    src = cms.InputTag("selectMuons") 
    )
process.muonSelection = cms.Sequence( process.patMuonsColl + process.selectMuons + process.countMuons )

# Event Selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
process.eventSelection = cms.Sequence( process.collisionEventSelectionPA )

# Path
process.selectionPath = cms.Path( process.triggerSelection + process.muonSelection) # + process.eventSelection Remove Collision Event Selection

# Options
process.source    = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.inputFiles ) )
#process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands =  cms.untracked.vstring('keep *', 'drop *_*_*_SKIM'),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('selectionPath') )
    )
process.o = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.selectionPath,process.endjob_step,process.o)
