import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023_UPC)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/hidata/HIRun2023A/HIEmptyBX/AOD/16Jan2024-v1/40000/7e3d5085-5839-4645-b4b6-9c87e6c234db.root'),
    )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v4')

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Empty BX triggers
    'HLT_HIL1NotBptxOR_v*', #0
    'HLT_HIL1UnpairedBunchBptxMinus_v*', #1
    'HLT_HIL1UnpairedBunchBptxPlus_v*', #2
]
process.eventFilterHLT = cms.Sequence(process.hltFilter)
process.eventFilterHLT_step = cms.Path(process.eventFilterHLT)

# Add the Particle tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna

process.eventAna = particleAna.clone(
  recoParticles = cms.InputTag(""),
  selectEvents = cms.string("eventFilterHLT_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_primaryVertexFilter',
      'Flag_hiClusterCompatibility',
  ),
  triggerInfo = cms.untracked.VPSet([
    # Empty BX triggers
    cms.PSet(path = cms.string('HLT_HIL1NotBptxOR_v*')), #0
    cms.PSet(path = cms.string('HLT_HIL1UnpairedBunchBptxMinus_v*')), #1
    cms.PSet(path = cms.string('HLT_HIL1UnpairedBunchBptxPlus_v*')), #2
  ]),
#   dataset = cms.untracked.string("Dataset_HIForward*"),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('ana_EmptyBX.root'))
process.p = cms.EndPath(process.eventAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilterHLT_step,
    process.p,
)

# Add the event selection filters
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.Flag_hiClusterCompatibility = cms.Path(process.hiClusterCompatibility)

eventFilterPaths = [ process.Flag_primaryVertexFilter , process.Flag_hiClusterCompatibility ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)
