import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(),
)


inputFiles='file:onia2MuMuPAT_regit_1000_1_yLu.root'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

#process.GlobalTag.globaltag = 'GR_P_V27A::All'    # prompt
process.GlobalTag.globaltag = 'STARTHI44_V12::All'

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
   
process.source.fileNames = cms.untracked.vstring(inputFiles)

process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')

process.HeavyIonGlobalParameters = cms.PSet(
         centralitySrc = cms.InputTag("hiCentrality"),
         centralityVariable = cms.string("HFtowers"),
         #centralityVariable = cms.string("HFhits"),
         nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
         #nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
      )

process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
process.centralityFilter.selectedBins = [2,3] # 0510 

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuonsSglTrg = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(QUALITY_CUTS + ' && ' + IN_ACCEPTANCE + ' && ' + DXYZ_CUTS + ' && ' + TAG_CUTS + ' && ' + "(!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()) || (!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()) || (!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty())"),
)


#--------------------------------------------------------------------
##    ____                   _____               _      ____            _               
##   | __ )  __ _ _ __ ___  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___ 
##   |  _ \ / _` | '__/ _ \   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |_) | (_| | | |  __/   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__,_|_|  \___|   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##                                                                                      
##   
process.goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("hiSelectedTracks"),
    cut = cms.string(TRACK_CUTS.replace("track.","")),
)

process.tkTracks  = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodTracks"),      
    particleType = cms.string("mu+"),
)
##### CHANGE CUTS : acceptance and quality of tracks
process.tkProbes = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("tkTracks"),
    cut = cms.string(IN_ACCEPTANCE),
)


##    ____  _                  _    _    _                    ____            _               
##   / ___|| |_ __ _ _ __   __| |  / \  | | ___  _ __   ___  |  _ \ _ __ ___ | |__   ___  ___ 
##   \___ \| __/ _` | '_ \ / _` | / _ \ | |/ _ \| '_ \ / _ \ | |_) | '__/ _ \| '_ \ / _ \/ __|
##    ___) | || (_| | | | | (_| |/ ___ \| | (_) | | | |  __/ |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__\__,_|_| |_|\__,_/_/   \_\_|\___/|_| |_|\___| |_|   |_|  \___/|_.__/ \___||___/
##                                                                                            
##   

process.goodSta = cms.EDFilter("TrackSelector",
    src = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    cut = cms.string("numberOfValidHits>0"),
)

process.staTracks = cms.EDProducer("ConcreteChargedCandidateProducer", 
    src  = cms.InputTag("goodSta"), 
    particleType = cms.string("mu+"),
)
process.staProbes = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("staTracks"),
    cut = cms.string("pt > 0"),
)



process.allTagsAndProbes = cms.Sequence( process.tagMuonsSglTrg +
                    process.goodTracks * process.tkTracks * process.tkProbes +
                    process.goodSta * process.staTracks * process.staProbes
)

##    ____               _               ____            _                   _____               _    _             
##   |  _ \ __ _ ___ ___(_)_ __   __ _  |  _ \ _ __ ___ | |__   ___  ___ _  |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##   | |_) / _` / __/ __| | '_ \ / _` | | |_) | '__/ _ \| '_ \ / _ \/ __(_)   | || '__/ _` |/ __| |/ / | '_ \ / _` |
##   |  __/ (_| \__ \__ \ | | | | (_| | |  __/| | | (_) | |_) |  __/\__ \_    | || | | (_| | (__|   <| | | | | (_| |
##   |_|   \__,_|___/___/_|_| |_|\__, | |_|   |_|  \___/|_.__/ \___||___(_)   |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                               |___/                                                                        |___/ 
##

process.tkToStaMatchAtVtx = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("tkTracks"), # all tk tracks
    matched = cms.InputTag("staTracks"),  # to all sta tracks
    algorithm = cms.string("byDirectComparison"), # using parameters at PCA
    srcTrack = cms.string("tracker"),  # trkTrack is a 'RecoChargedCandidate'
    srcState = cms.string("atVertex"), # matching done at the outermost layer of the tracker (?)
    matchedTrack = cms.string("tracker"),
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(1.),   # large range in DR
    maxDeltaEta      = cms.double(0.2),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(3),
    sortBy           = cms.string("deltaR"),
)

process.tkToStaMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("tkTracks"), # all tk tracks
    matched = cms.InputTag("staTracks"),  # to all sta tracks
    algorithm = cms.string("byPropagatingSrc"), # using parameters at PCA
    srcTrack = cms.string("tracker"),  # trkTrack is a 'RecoChargedCandidate'
    srcState = cms.string("outermost"), # matching done at the outermost layer of the tracker (?)
    matchedTrack = cms.string("muon"),
    matchedState = cms.string("innermost"),
    maxDeltaR        = cms.double(1.),   # large range in DR
    maxDeltaEta      = cms.double(0.2),  # small in eta, which is more precise
    maxDeltaLocalPos = cms.double(100),
    maxDeltaPtRel    = cms.double(3),
    sortBy           = cms.string("deltaR"),
)
process.tkPassingSta = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("tkProbes"),
    match = cms.InputTag("tkToStaMatch"),
)

process.tkPassingStaAtVtx = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("tkProbes"),
    match = cms.InputTag("tkToStaMatchAtVtx"),
)


process.allPassingProbes = cms.Sequence(
    process.tkToStaMatch * process.tkPassingSta +
    process.tkToStaMatchAtVtx * process.tkPassingStaAtVtx
)


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
process.tpTkSta = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuonsSglTrg@+ tkProbes@-"), # charge conjugate states are implied
    cut   = cms.string("2 < mass < 5"),
)


#####
## Tk from Sta
process.fitTkFromSta = cms.EDAnalyzer("TagProbeFitTreeProducer",
     tagProbePairs = cms.InputTag("tpTkSta"),
     arbitration   = cms.string("OneProbe"),
     # probe variables: all useful ones
     variables = cms.PSet(
     KinematicVariables,
     #staOnlyVariables,
     ),
     # choice of tag and probe pairs, and arbitration
     # choice of what defines a 'passing' probe
     flags = cms.PSet(
        passing = cms.InputTag("tkPassingSta"),
        passingAtVtx = cms.InputTag("tkPassingStaAtVtx"),
       # outerValidHits  = cms.string("? outerTrack.isNull() ? 0 : outerTrack.numberOfValidHits > 0"), 
     ),
     tagVariables = cms.PSet(
     pt  = cms.string("pt"),
     eta = cms.string("eta"),
     abseta = cms.string("abs(eta)"),
     ),
     tagFlags     = cms.PSet(
     ),
     pairVariables = cms.PSet(
     pt  = cms.string("pt"),
     y = cms.string("rapidity"),
     absy = cms.string("abs(rapidity)"),
     ),
     pairFlags = cms.PSet(),
     isMC           = cms.bool(False),
     #addRunLumiInfo = cms.bool(True),
     allProbes     = cms.InputTag("tkProbes"),
     eventWeight = cms.double(1198.723),
)


##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuonsSglTrg +
    process.tpTkSta +
    process.fitTkFromSta
)

process.tagAndProbe = cms.Path(
    process.centralityFilter *
    process.allTagsAndProbes *
    process.allPassingProbes *
    process.tnpSimpleSequence
)



process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_PbPb_GenTrkSTA_MC.root"))


