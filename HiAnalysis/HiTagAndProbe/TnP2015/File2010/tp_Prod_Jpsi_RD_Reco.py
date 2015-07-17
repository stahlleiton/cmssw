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
process.GlobalTag.globaltag = 'START39_V7HI::All'

process.source = cms.Source("PoolSource", 
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
      _input_
    ),
)
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(100),
    input = cms.untracked.int32(-1),
)    

#process.load("MuonAnalysis.MuonAssociators.hiPatMuonsWithTrigger_cff")
process.load("MuonAnalysis.MuonAssociators.hiRDPatMuonsWithTrigger_cff")

from MuonAnalysis.MuonAssociators.hiRDPatMuonsWithTrigger_cff import *
useL1MatchingWindowForSinglets(process)
changeTriggerProcessName(process, "HLT")   # Custom re-run HLT

### Torsten's cut
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
PASS_TAG_1MU = "(!triggerObjectMatchesByFilter('hltL1sHIL1SingleMu3').empty())"
PASS_TAG_2MU = "(!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty())"
TAG_CUTS_MU = "isGlobalMuon && " + TRACK_CUTS +' && '+ PASS_TAG_2MU

ptMinCut = 'pt > 2 || (abs(eta) > 1 && p > 2)';

process.tagMuons = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("patMuonsWithTrigger"),
        cut = cms.string(TAG_CUTS_MU), 
)

process.goodTracks = cms.EDFilter("TrackSelector",
    src  = cms.InputTag("hiGlobalPrimTracks"),     # L4 
    cut = cms.string(ptMinCut),
)

process.betterTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("goodTracks"),
    cut = cms.string(TRACK_CUTS.replace("track.","")), # this is a Track object, so I have to remove the 'track.'
)
process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
        src  = cms.InputTag("betterTracks"),     # L4 
        particleType = cms.string("mu+"),
) 
process.tkProbes = cms.EDFilter("CandViewRefSelector",
        src = cms.InputTag("tkTracks"),
        cut = cms.string(PT_ETA_CUTS),
)

# Passing Probe collection
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.6'),
        decay = cms.string('tagMuons@+ tkProbes@-')
)

### Match to MC truth
process.tagMuonsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        src = cms.InputTag("tagMuons"),
        matched = cms.InputTag("genParticles"),
        pdgId = cms.vint32(13),
        distMin = cms.double(0.3),
)
process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "tkProbes")

### Passing Probe ###
process.glbMuons = cms.EDFilter("PATMuonRefSelector",
        src = cms.InputTag("patMuonsWithTrigger"),
        cut = cms.string(PASSING_GLB_CUT), 
)
### Sta Muon
### option 3
process.probeMuonsSta3 = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("isStandAloneMuon && " + PT_ETA_CUTS),
)

### opion 1
process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)
from PhysicsTools.PatAlgos.tools.helpers import *
process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(process, process.patMuonsWithTriggerSequence, "Sta")

process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerSta"),
    cut = cms.string("outerTrack.isNonnull && " + PT_ETA_CUTS), # no real cut now
)

### option 2
process.staTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    particleType = cms.string("mu+"),
)
process.probeMuonsSta2 = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("staTracks"),
    cut = cms.string(PT_ETA_CUTS),
)

process.probeMuonsMCMatchSta = process.tagMuonsMCMatch.clone(src = "probeMuonsSta2")
process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta2@-", cut = "2 < mass < 5")
#process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = "2 < mass < 5")

process.tkToGlbMatch = cms.EDProducer("MatcherUsingTracks",
        src     = cms.InputTag("tkTracks"), # all tracks are available for matching
        matched = cms.InputTag("glbMuons"), # to all global muons
        algorithm = cms.string("byDirectComparison"), # check that they
        srcTrack = cms.string("tracker"),             # have the same 
        srcState = cms.string("atVertex"),            # tracker track
        matchedTrack = cms.string("tracker"),         # can't check ref
        matchedState = cms.string("atVertex"),        # because of the
        maxDeltaR        = cms.double(0.01),          # embedding.
        maxDeltaLocalPos = cms.double(0.01),
        maxDeltaPtRel    = cms.double(0.01),
        sortBy           = cms.string("deltaR"),
)
process.tkPassingGlb = cms.EDProducer("MatchedCandidateSelector",
        src   = cms.InputTag("tkProbes"),
        match = cms.InputTag("tkToGlbMatch"),
)

process.staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("staTracks"),
    matched = cms.InputTag("tkTracks"),  
    algorithm = cms.string("byDirectComparison"), 
    srcTrack     = cms.string("muon"),    
    srcState = cms.string("atVertex"), 
    matchedTrack = cms.string("tracker"), 
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(1.),   # large range in DR (we can tighten it later)
    maxDeltaEta      = cms.double(0.4),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(5),   # |pt(sta) - pt(tk)|/pt(tk)
    sortBy           = cms.string("deltaR"),
)
    
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("probeMuonsSta2"),
    match = cms.InputTag("staToTkMatch"),
)

# Make the fit tree and save it in the "MuonID" directory
process.MuonID = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      HLTL1Mu3 = cms.string("!triggerObjectMatchesByFilter('hltL1sHIL1SingleMu3').empty()"),
      HLTL2Mu3 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    #isMC = cms.bool(True),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("tkProbes"),
)

process.MuonIDsta = process.MuonID.clone(
    tagProbePairs = "tpPairsSta",
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        pt  = cms.string("pt"),
        p   = cms.string("p"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        PassingSta = cms.InputTag("staPassingTk"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      HLTL1Mu3 = cms.string("!triggerObjectMatchesByFilter('hltL1sHIL1SingleMu3').empty()"),
      HLTL2Mu3 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      Acc_JPsi   = cms.string(JPSI_ACCEPTANCE_CUT),
    ),
    #isMC = cms.bool(True),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = "probeMuonsSta2",
)
    
process.tnpSimpleSequence = cms.Sequence(
        process.tagMuons * #process.tagMuonsMCMatch *
        process.glbMuons * 
        process.goodTracks * process.betterTracks *
        process.tkTracks * process.tkProbes * #process.probeMuonsMCMatch *
        process.staTracks * process.probeMuonsSta2 * #process.probeMuonsMCMatchSta +
        process.tkToGlbMatch * process.tkPassingGlb *
        process.staToTkMatch * process.staPassingTk *
        process.tpPairs * 
        process.tpPairsSta *
        process.MuonID * 
        process.MuonIDsta
)

process.tagAndProbe = cms.Path(
        process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("_output_file_"))
