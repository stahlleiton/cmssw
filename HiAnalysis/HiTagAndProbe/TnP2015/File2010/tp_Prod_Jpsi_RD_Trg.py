import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('START39_V7HI::All')

process.source = cms.Source("PoolSource", 
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
      _input_
    ),
)
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1000),
    input = cms.untracked.int32(-1),
)    

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.load("MuonAnalysis.MuonAssociators.hiRDPatMuonsWithTrigger_cff")

from MuonAnalysis.MuonAssociators.hiRDPatMuonsWithTrigger_cff import *
useL1MatchingWindowForSinglets(process)
changeTriggerProcessName(process, "HLT")   # Custom re-run HLT
#changeTriggerProcessName(process, "HLT2")   # Custom re-run HLT
#changeTriggerProcessName(process, "DATAMIX")   # Custom re-run HLT

TRACK_CUTS = ("track.numberOfValidHits > 11 && track.hitPattern.pixelLayersWithMeasurement > 1 "+
              "&& track.normalizedChi2 < 4 ") #+
              #"&& abs(track.d0) < 3 && abs(track.dz) < 30")
PT_ETA_CUTS = "(pt > 3 || (abs(eta)>1 && p > 2.6))" ## the enclosing () are very important, because there's an "||"
JPSI_ACCEPTANCE_CUT = ("(       abs(eta) <= 1.3  && pt > 3.3 || "+
                       "  1.3 < abs(eta) <= 2.2  && p >  2.9 || "+
                       "  2.2 < abs(eta) <= 2.4  && pt > 0.8 )  ")

PASSING_GLB_CUT = ("isGlobalMuon && globalTrack.normalizedChi2 < 20  && "+
                   "globalTrack.hitPattern.numberOfValidMuonHits > 0 && "+
                   "muonID('TrackerMuonArbitrated') &&  muonID('TMLastStationAngTight')")

PASS_HLT_L2DoubleMu3 = "(!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty())"
PASS_HLT_L1DoubleMuOpen = "(!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty())"
PASS_TAG_MU = "(!triggerObjectMatchesByFilter('hltL1sHIL1SingleMu3').empty())"
TAG_CUTS_MU = "isGlobalMuon && " + TRACK_CUTS +' && '+ PASS_TAG_MU

process.recoMuFilter = cms.EDFilter("CandViewCountFilter", 
        src = cms.InputTag("muons"),    
        minNumber = cms.uint32(1)
)

process.tagMuons = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("patMuonsWithTrigger"),
        cut = cms.string(TAG_CUTS_MU), 
        #cut = cms.string("isGlobalMuon && !triggerObjectMatchesByFilter('hltL1sHIL1SingleMu3').empty()"), 
)

process.probeMuons = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("patMuonsWithTrigger"),
        cut = cms.string("isGlobalMuon" + TRACK_CUTS + " && " + PT_ETA_CUTS),
        #cut = cms.string("isGlobalMuon"),
)

# Passing Probe collection and tpPairs
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.6'),
        decay = cms.string('tagMuons@+ probeMuons@-')
)

# Make the fit tree and save it in the "MuonID" directory
process.TrgEff = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      p   = cms.string("p"),
      phi = cms.string("phi"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    flags = cms.PSet(
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(False),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuons"),
    #addRunLumiInfo = cms.bool(True),
)

process.tnpSimpleSequence = cms.Sequence(
        process.tagMuons   * 
        process.probeMuons * 
        process.tpPairs    *
        process.TrgEff
)

process.tagAndProbe = cms.Path(
#        process.recoMuFilter +
#        process.patMuonsWithTriggerSequence *
        process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("_output_file_"))
