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
options.inputFiles = "/store/mc/PhaseIIMTDTDRAutumn18DR/D0_PiK_prompt_pt0_y4_5p5TeV_TuneCP5_Pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/30000/FFDBA739-F97F-B347-BAB9-91EC1A7F2CE1.root"
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

# Path and EndPath definitions
process.load('HiMTDAnalysis.TrackAnalysis.timeAnalyzer_cfi')
process.anaPath = cms.Path( process.cent_seq * process.timeAnaSeq )
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
process.timeAna.doDiTracks         = cms.bool(False)
###############################################################################################

#Options:
process.source       = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.inputFiles ) )
process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.maxEvents    = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step, process.anaPath)
