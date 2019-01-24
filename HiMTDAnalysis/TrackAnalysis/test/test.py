import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras

process = cms.Process('MTDAnalysis',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "HiMTDTree.root"
options.inputFiles = "file:Hydjet_RECO.root"
options.maxEvents  = 100 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)

process.load('HiMTDAnalysis.TrackAnalysis.timeAnalyzer_cfi')
process.load('RecoMTD.TrackExtender.trackExtenderWithMTD_cfi')

process.load('RecoLocalFastTime.FTLRecProducers.mtdTrackingRecHits_cfi')
process.load('RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi')

process.anaPath = cms.Path( process.mtdClusters * process.mtdTrackingRecHits * process.trackExtenderWithMTD + process.timeAnaSeq )

#Options:
process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.inputFiles ) )
process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.maxEvents    = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# Schedule definition
process.schedule = cms.Schedule(process.anaPath,process.endjob_step)
