import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)

inputFiles='file:F00776F3-3075-E211-AD4E-0025901D5E10.root'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

#process.GlobalTag.globaltag= 'STARTHI53_V28::All'
process.GlobalTag.globaltag= 'GR_P_V43D::All' #Data

IN_ACCEPTANCE = '( (abs(eta)<1.0 && pt>=3.4) || (1.0<=abs(eta)<1.5 && pt>=5.8-2.4*abs(eta)) || (1.5<=abs(eta)<2.4 && pt>=3.3667-7.0/9.0*abs(eta)) )'

# several selection cuts
TRACK_CUTS    = "track.numberOfValidHits > 10 && track.normalizedChi2 < 4 && track.hitPattern.pixelLayersWithMeasurement > 0"
GLB_CUTS      = "isGlobalMuon && globalTrack.normalizedChi2 < 20"
QUALITY_CUTS  =  "(" + GLB_CUTS + ' && ' + TRACK_CUTS + ")"
DXYZ_CUTS = "abs(dB) < 3 && abs(track.dz) < 15"
TAG_CUTS = "isTrackerMuon && muonID('TrackerMuonArbitrated')"

staOnlyVariables = cms.PSet(
staQoverP      = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverp"),
staQoverPerror = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverpError"),
staValidStations = cms.string("? outerTrack.isNull() ? -1 : outerTrack.hitPattern.muonStationsWithValidHits()"),
staNumValidHits = cms.string("? outerTrack.isNull() ? -1 : outerTrack.hitPattern.numberOfValidMuonHits()"),
)

TrigTagFlags = cms.PSet(
    HLTL2Mu3 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()"),
    HLTL2Mu7 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()"),
    HLTL2Mu15= cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty()"),
)

TrackQualityVariables = cms.PSet(
    dB          = cms.string("dB"),
    tkValidHits = cms.string("? track.isNull ? 0 : track.numberOfValidHits"),
    tkValidPixelHits = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidPixelHits"),
    tkPixelLay  = cms.string("? track.isNull ? 0 : track.hitPattern.pixelLayersWithMeasurement"),
    tkExpHitIn  = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsInner.numberOfLostHits"),
    tkExpHitOut = cms.string("? track.isNull ? 0 : track.trackerExpectedHitsOuter.numberOfLostHits"),
    tkHitFract  = cms.string("? track.isNull ? 0 : track.hitPattern.numberOfValidHits/(track.hitPattern.numberOfValidHits+track.hitPattern.numberOfLostHits+track.trackerExpectedHitsInner.numberOfLostHits+track.trackerExpectedHitsOuter.numberOfLostHits)"),
    tkChi2 = cms.string("? track.isNull ? -1 : track.normalizedChi2"),
    tkPtError = cms.string("? track.isNull ? -1 : track.ptError"),
    tkSigmaPtOverPt = cms.string("? track.isNull ? -1 : track.ptError/track.pt"),
)
GlobalTrackQualityVariables = cms.PSet(
    glbChi2 = cms.string("? globalTrack.isNull ? -1 : globalTrack.normalizedChi2"),
    glbPtError = cms.string("? globalTrack.isNull ? -1 : globalTrack.ptError"),
    glbSigmaPtOverPt = cms.string("? globalTrack.isNull ? -1 : globalTrack.ptError/globalTrack.pt"),
)




   
process.source.fileNames = cms.untracked.vstring(inputFiles)


process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_PAMu3_v*','HLT_PAMu7_v*','HLT_PAMu12_v*' )
process.triggerResultsFilter.l1tResults = ''
process.triggerResultsFilter.throw = True
process.triggerResultsFilter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
process.HLTMu   = process.triggerResultsFilter.clone(triggerConditions = ['HLT_PAMu3_v*','HLT_PAMu7_v*','HLT_PAMu12_v*'])


## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    mergeCaloMuons = cms.bool(False), # AOD
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("track.isNonnull && pt > 1.5"),
    caloMuonsCut = cms.string("pt > 1.5"),
    tracksCut    = cms.string("pt > 1.5"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
useExtendedL1Match(process)
addHLTL1Passthrough(process)

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuonsSglTrg = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(QUALITY_CUTS + ' && ' + IN_ACCEPTANCE + ' && ' + DXYZ_CUTS + ' && ' + TAG_CUTS + " && (!triggerObjectMatchesByPath('HLT_PAMu3_v*').empty() && !triggerObjectMatchesByFilter('hltL3fL2sMu3L3Filtered3').empty()) || (!triggerObjectMatchesByPath('HLT_PAMu7_v*').empty() && !triggerObjectMatchesByFilter('hltL3fL2sMu7L3Filtered7').empty()) || (!triggerObjectMatchesByPath('HLT_PAMu12_v*').empty() && !triggerObjectMatchesByFilter('hltL3fL2sMu12L3Filtered12').empty())"),
)

process.probeMuonsGenTrk = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(TRACK_CUTS + ' && ' + IN_ACCEPTANCE),
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('2. < mass < 5.'),
    decay = cms.string('tagMuonsSglTrg@+ probeMuonsGenTrk@-')
)


# check if there is at least one Tag and Probe pair
process.tpPairsFilter = cms.EDFilter("CandViewCountFilter",
     src = cms.InputTag('tpPairs'),
     minNumber = cms.uint32(1),
)


process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
     # choice of tag and probe pairs, and arbitration
     tagProbePairs = cms.InputTag("tpPairs"),
     arbitration   = cms.string("OneProbe"),
     # probe variables: all useful ones
     variables = cms.PSet(
     KinematicVariables,
     staOnlyVariables,
     ),
     flags = cms.PSet(
     isSTA = cms.string("isStandAloneMuon"),
     outerValidHits  = cms.string("? outerTrack.isNull() ? 0 : outerTrack.numberOfValidHits > 0"),
     ),
     tagVariables = cms.PSet(
     TrackQualityVariables,
     GlobalTrackQualityVariables,
     pt  = cms.string("pt"),
     eta = cms.string("eta"),
     abseta = cms.string("abs(eta)"),
     l2dr  = cms.string("? triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() ? 999 : "+
                            " deltaR( eta, phi, " +
                            "         triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).eta, "+
                            "         triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).phi ) "),

     ),
     tagFlags     = cms.PSet(
     TrigTagFlags,
     ),
     pairVariables = cms.PSet(
     pt  = cms.string("pt"),
     y = cms.string("rapidity"),
     absy = cms.string("abs(rapidity)"),
     ),
     pairFlags = cms.PSet(),
     isMC           = cms.bool(False),
     #addRunLumiInfo = cms.bool(True),
     allProbes     = cms.InputTag("probeMuonsGenTrk"),
     # addCentralityInfo = cms.bool(False) 
)



process.tnpSimpleSequence = cms.Sequence(
    process.tagMuonsSglTrg +
    process.probeMuonsGenTrk +
    process.tpPairs    +
    process.tpPairsFilter +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.PAcollisionEventSelection *
    process.HLTMu *
    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence *
    process.tnpSimpleSequence
)



process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_pp_GenTrkSTA_Data.root"))
