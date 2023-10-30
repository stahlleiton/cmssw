### HiForest Configuration
# Input: AOD
# Type: mc

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process('HiForest',Run3_pp_on_PbPb_2023)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, AOD, 132X, mc")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/SUPERCHIC/2023_10_23/SUPERCHIC_5p36TeV_2023Run3/ChiC_SUPERCHIC_5p36TeV_2023Run3_HIFwd_RECO_AOD_2023_10_23/231022_221941/0000/SUPERCHIC_chiC_RECO_1.root'
    ),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(25)
    )

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_HI_v5', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# Define centrality binning
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db"),
        label = cms.untracked.string("HFtowers")
        ),
    ])

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestAOD_MC.root"))

###############################################################################

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.hltobject.triggerObjects = cms.InputTag("patTrigger")
process.triggerSequence = cms.Sequence(process.patTrigger + process.hltobject + process.hltanalysis + process.l1object)

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.hiEvtAnalyzer.pfCandidateSrc = cms.InputTag('particleFlow')
process.hiEvtAnalyzer.Vertex = cms.InputTag("offlinePrimaryVertices")
process.hiEvtAnalyzer.doHFfilters = cms.bool(False)
process.eventSequence = cms.Sequence(process.centralityBin * process.hiEvtAnalyzer)

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
process.hltobject.triggerNames = trigger_list_data_2023_skimmed

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.particleFlowAnalyser.pfCandidateSrc = cms.InputTag('particleFlow')
process.particleFlowAnalyser.ptMin = cms.double(0)
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.load("PhysicsTools.PatAlgos.producersLayer1.photonProducer_cfi")
process.photonSequence = cms.Sequence(process.patPhotons * process.ggHiNtuplizer)
process.ggHiNtuplizer.doMuons = cms.bool(False)
process.ggHiNtuplizer.useValMapIso = cms.bool(False) # segmentation fault
process.ggHiNtuplizer.recHitsEB = cms.untracked.InputTag("reducedEcalRecHitsEB")
process.ggHiNtuplizer.recHitsEE = cms.untracked.InputTag("reducedEcalRecHitsEE")
process.ggHiNtuplizer.pileupSrc = cms.InputTag("addPileupInfo")
process.ggHiNtuplizer.vertexSrc = cms.InputTag("offlinePrimaryVertices")
process.ggHiNtuplizer.electronSrc = cms.InputTag("gedGsfElectrons")
process.ggHiNtuplizer.photonSrc = cms.InputTag("patPhotons")
process.ggHiNtuplizer.conversionsSrc = cms.InputTag('allConversions')
process.ggHiNtuplizer.particleFlowCollection = cms.InputTag("particleFlow")
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
process.PbPbTracks.vertexSrc = cms.InputTag("offlinePrimaryVertices")
process.PbPbTracks.trackSrc = cms.InputTag("generalTracks")
process.trackSequencePbPb = cms.Sequence(process.PbPbTracks)
# muons (FTW)
process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi")
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.patMuons.addGenMatch   = cms.bool(False)
process.patMuons.embedGenMatch = cms.bool(False)
process.muonSequence = cms.Sequence(process.patMuons * process.muonAnalyzer)
process.muonAnalyzer.muonSrc = cms.InputTag("patMuons")
process.muonAnalyzer.vertexSrc = cms.InputTag("offlinePrimaryVertices")

#############################
# Gen Analyzer
#############################
process.load('PhysicsTools.PatAlgos.slimming.packedGenParticles_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.HiGenAnalyzer_cfi')
from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import prunedGenParticles
process.prunedGenParticles = prunedGenParticles.clone(select = cms.vstring("keep *"))
process.genMap = prunedGenParticles.clone(select = cms.vstring("keep *"), src =  cms.InputTag("prunedGenParticles"))
from PhysicsTools.PatAlgos.slimming.packedGenParticles_cfi import packedGenParticles
process.packedGenParticles = packedGenParticles.clone(inputCollection = cms.InputTag("prunedGenParticles"), map = cms.InputTag("genMap"))
process.HiGenParticleAna.ptMin = cms.untracked.double(0.4) # default is 5
process.HiGenParticleAna.etaMax = cms.untracked.double(5.) # default is 2.5
process.genSequence = cms.Sequence(process.prunedGenParticles * process.genMap * process.packedGenParticles * process.HiGenParticleAna)


###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.genSequence +
    process.eventSequence +
    process.triggerSequence +
    process.trackSequencePbPb +
    process.particleFlowAnalyser +
    process.photonSequence +
    process.muonSequence
    )

#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.primaryVertexFilter.src = cms.InputTag("offlinePrimaryVertices")
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.pAna = cms.EndPath(process.skimanalysis)
