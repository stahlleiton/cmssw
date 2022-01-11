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
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "HiMTDTree.root"
options.inputFiles = "/store/user/anstahll/MTD/MC/PY8EGun_TuneCP5_Deuteron_p12_eta3p2_Hydjet_PbPb_5p02TeV_RECO_20201007/MTD/PY8EGun_TuneCP5_Deuteron_p12_eta3p2_Hydjet_PbPb_5p02TeV_RECO_20201007/201010_183358/0000/step3_PY8EGun_TuneCP5_Deuteron_p12_eta3p2_Hydjet_PbPb_5p02TeV_RECO_20201007_1.root"
#options.inputFiles = "file:/afs/cern.ch/user/a/anstahll/work/MTD/Analyzer/2020/CMSSW_10_4_0_mtd5/src/RECO.root"
#options.inputFiles = "/store/user/anstahll/MTD/MC/20201103/PY8EGun_TuneCP5_3LambdaH_p12_eta3p2_To3HePi_Pythia8_pp_14TeV_RECO_20201103/MTD/PY8EGun_TuneCP5_3LambdaH_p12_eta3p2_To3HePi_Pythia8_pp_14TeV_RECO_20201103/201103_075516/0000/step3_PY8EGun_TuneCP5_3LambdaH_p12_eta3p2_To3HePi_Pythia8_pp_14TeV_RECO_20201103_1.root"
#options.inputFiles = "file:/afs/cern.ch/user/a/anstahll/work/MTD/Analyzer/2020/CMSSW_10_4_0_mtd5/src/test.root"
options.maxEvents  = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Centrality setup
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetTuneCP5MTD_v1040mtd4x1_mc"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = True
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = False
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.hiCentrality.srcTracks = cms.InputTag("generalTracks")
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices4D")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.hiCentrality.srcEBhits = cms.InputTag("HGCalRecHit","HGCHEBRecHits")
process.hiCentrality.srcEEhits = cms.InputTag("HGCalRecHit","HGCEERecHits")

process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

# Define SIM-RECO association
process.load('SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff')
process.quickTrackAssociatorByHits.ComponentName = cms.string('quickTrackAssociatorByHits')
process.simAssocSeq = cms.Sequence(process.tpClusterProducer * process.quickTrackAssociatorByHits * process.trackingParticleRecoTrackAsssociation)

# Path and EndPath definitions
process.load('HiMTDAnalysis.TrackAnalysis.timeAnalyzer_cfi')
process.timeAnaOld = process.timeAna.clone()
process.timeAnaSeq = cms.Sequence( process.timeAna * process.timeAnaOld )
process.anaPath = cms.Path( process.cent_seq * process.simAssocSeq * process.timeAnaSeq )
###############################################################################################

# MTD RE-RECO
process.reconstruction_step = cms.Path()
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.pfPileUpIso.PFCandidates = cms.InputTag("particleFlowPtrs")
process.pfNoPileUpIso.bottomCollection = cms.InputTag("particleFlowPtrs")
process.reconstruction_step += process.mtdClusters
process.reconstruction_step += process.mtdTrackingRecHits
process.trackExtenderWithMTD.UseVertex = cms.bool(True) #run trackExtender using vertex constrain
process.trackExtenderWithMTD.DZCut = 0.3
process.trackExtenderWithMTD.fixedT0Error = cms.double(0.035)
process.reconstruction_step += process.trackExtenderWithMTD
process.tofPID.vtxsSrc = cms.InputTag('offlinePrimaryVertices4D')
process.tofPID.fixedT0Error = cms.double(0.035) #put a constant 0.035 [ns] error for each track (cannot
process.reconstruction_step += process.tofPID

process.timeAna.trackBetaTag       = cms.InputTag("trackExtenderWithMTD:generalTrackBeta:MTDAnalysis")
process.timeAna.trackT0Tag         = cms.InputTag("trackExtenderWithMTD:generalTrackt0:MTDAnalysis")
process.timeAna.trackSigmaT0Tag    = cms.InputTag("trackExtenderWithMTD:generalTracksigmat0:MTDAnalysis")
process.timeAna.trackTMTDTag       = cms.InputTag("trackExtenderWithMTD:generalTracktmtd:MTDAnalysis")
process.timeAna.trackSigmaTMTDTag  = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd:MTDAnalysis")
process.timeAna.trackMomTag        = cms.InputTag("trackExtenderWithMTD:generalTrackp:MTDAnalysis")
process.timeAna.trackPathLengthTag = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength:MTDAnalysis")
process.timeAna.tofPIDT0Tag        = cms.InputTag("tofPID:t0:MTDAnalysis")
process.timeAna.tofPIDSigmaT0Tag   = cms.InputTag("tofPID:sigmat0:MTDAnalysis")
###############################################################################################

#Options:
process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.inputFiles ) )
process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.maxEvents    = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step, process.anaPath)
