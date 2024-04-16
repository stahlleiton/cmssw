import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
process = cms.Process('RECO',Run3_2023)

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_DataMapper_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("UserCode.EnergyLossPID.EnergyLossProducer_cff")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger = cms.Service("MessageLogger",
  cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    DEBUG     = cms.untracked.PSet(
      limit = cms.untracked.int32(0)
    ),
    FwkReport = cms.untracked.PSet(
      optionalPSet = cms.untracked.bool(True),
      reportEvery = cms.untracked.int32(100),
      limit = cms.untracked.int32(10000)
    )
  ),
  destinations = cms.untracked.vstring('cerr'),
)
 
# Source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
     '/store/hidata/HIRun2023A/HIForward0/RAW/v1/000/375/703/00000/0abdd956-f86e-41e0-8bb5-e7a0e0804c47.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet(
)

# Steps
process.ra2d_step = cms.Path(process.RawToDigi)
process.reco_step = cms.Path(process.reconstruction)

# Schedule
process.schedule = cms.Schedule(process.ra2d_step,
                                process.reco_step,
                                process.produceEnergyLoss)

# Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v4', '')

# Some extra
process.hcalDigis.saveQIE10DataNSamples = cms.untracked.vint32( 6) 
process.hcalDigis.saveQIE10DataTags = cms.untracked.vstring( "MYDATA" )

# Some extra
from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process, new="rawDataMapperByLabel", old="rawDataCollector")

# Early deletion
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

