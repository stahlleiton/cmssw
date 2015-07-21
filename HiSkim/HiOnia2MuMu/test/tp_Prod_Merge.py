import FWCore.ParameterSet.Config as cms

process = cms.Process("TnPTrg")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V27A::All'    # prompt 

# centrality stuff
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')

process.HeavyIonGlobalParameters = cms.PSet(
      centralityVariable = cms.string("HFtowers"),
      nonDefaultGlauberModel = cms.string(""),
      centralitySrc = cms.InputTag("hiCentrality")
      )


process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
#process.centralityFilter.selectedBins = [0,1,2,3] # 0 - 10 %
#process.centralityFilter.selectedBins = [4,5,6,7] # 10 - 20 %
#process.centralityFilter.selectedBins = [8,9,10,11,12,13,14,15,16,17,18,19] # 20 - 50 %
#process.centralityFilter.selectedBins = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] # 50 - 100 %
process.centralityFilter.selectedBins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] # MinBias

process.source = cms.Source("PoolSource", 
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
       "file:tnp_regit.root",
    ),
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
)    

# Selection : acceptance and muId, for innerTrack and globalTrack
IN_ACCEPTANCE = '( (abs(eta)<1.0 && pt>=3.4) || (1.0<=abs(eta)<1.5 && pt>=5.8-2.4*abs(eta)) || (1.5<=abs(eta)<2.4 && pt>=3.3667-7.0/9.0*abs(eta)) )'
# tracke muon cuts 
TRACK_CUTS    = "isTrackerMuon && innerTrack.numberOfValidHits > 10 && innerTrack.normalizedChi2 < 4 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && abs(innerTrack.dB(\"PV3D\")) < 3 && abs(innerTrack.dz(\"PV3D\")) < 15  && muonID('TrackerMuonArbitrated')"
# global muon cuts
GLB_CUTS      = "isGlobalMuon && globalTrack.normalizedChi2 < 20"
QUALITY_CUTS  =  "(" + GLB_CUTS + ' && ' + TRACK_CUTS + ")"

MuonIDFlags = cms.PSet(
    GlobalMu  = cms.string("isGlobalMuon"),
    TrackerMu = cms.string("isTrackerMuon"),
    TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
    TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
    TrackCuts	= cms.string(TRACK_CUTS),
    GlobalCuts	= cms.string(GLB_CUTS),
    QualityMu	= cms.string(QUALITY_CUTS)
  
)

########## TAG & Pair DEFINITIONS!
# from the input collection, make sure they pass also the MachByFilter
process.tagMuonsSglTrgNew = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("tagMuonsSglTrg"),
        cut = cms.string("(!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()) || (!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()) || (!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty())"),
)

# pairs for trigger efficiency (probe muons are created in the skim step, and pass all the analysis cuts except the trigger matching)
process.tpPairsTrigNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5'),
        decay = cms.string('tagMuonsSglTrgNew@+ probeMuons@-')
)

# pairs for muId (and muId*trigger) efficiency (probe muons are created in the skim step, and the only conditions: isglobal and inACCeptance)
process.tpPairsMuIdNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5'),
        decay = cms.string('tagMuonsSglTrgNew@+ probeMuonsTrk@-')
)


# pairs for tracking efficiency (probe muons are created in the skim step, and the only condition: isSta)
process.tpPairsStaNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 1 < mass < 5 '),
        decay = cms.string('tagMuonsSglTrgNew@+ probeMuonsSta@-')
)


################################################################################
# a bunch of predefined things that you want to keep in (prety awesome!)
# https://github.com/CMS-HIN-dilepton/cmssw/blob/CMSSW_4_4_X_Lxyz/MuonAnalysis/TagAndProbe/python/common_variables_cff.py
from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

##############
# Trigger efficiency by itself: make the fit tree and save it in the "Trigger" directory
process.MuonTrg = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsTrigNew"),
    arbitration   = cms.string("None"), # have unique tag-probe for each event
    variables = cms.PSet( # probe variables that will be stored in the output tree
    KinematicVariables,
	L1Variables
    ),
    flags = cms.PSet( # passing probe condition (the efficiency we want to probe); decision 1 or 0 wil be stored in the output tree
      HLTL1v0 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
      HLTL1v1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL1v2 = cms.string("(!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty())"),
    ),
    tagVariables = cms.PSet( # tag variables that will be stored in the output tree
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      l2dr  = cms.string("? triggerObjectMatchesByCollection('hltL2MuonCandidates').empty() ? 999 : "+
                        " deltaR( eta, phi, " +
                        "         triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).eta, "+
                        "         triggerObjectMatchesByCollection('hltL2MuonCandidates').at(0).phi ) ")
        ),
    tagFlags = cms.PSet(# tag cut; decision 1 or 0 wil be stored in the output tree
      HLTL2Mu3 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()"),
      HLTL2Mu7 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()"),
      HLTL2Mu15= cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty()"),
    ),
    pairVariables = cms.PSet( #pair variables
      pt  = cms.string("pt"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
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

##############
# Tracking efficiency by itself: make the fit tree and save it in the "Tracking" directory
process.MuonTrk = cms.EDAnalyzer("TagProbeFitTreeProducer",
        tagProbePairs = cms.InputTag("tpPairsStaNew"),
        arbitration   = cms.string("OneProbe"),
        variables = cms.PSet(
        KinematicVariables, 
        StaOnlyVariables, # undeclared!
            ## track matching variables
            tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
            tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
            tk_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ","deltaR"),
        tk_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ","deltaEta"),
    ),
    flags = cms.PSet(
        outerValidHits  = cms.string("outerTrack.numberOfValidHits > 0"),
        Tk              = cms.string("track.isNonnull"),
        StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
        GlobalMu        = cms.string("isGlobalMuon"),
        TrackerMu       = cms.string("isTrackerMuon"),
    	TrackCuts       = cms.string(TRACK_CUTS),
        GlobalCuts	    = cms.string(GLB_CUTS),
       	PassingSta = cms.string("isGlobalMuon && " + GLB_CUTS),
    ),
    tagVariables = cms.PSet( # tag variables that will be stored in the output tree
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)")
       ),
    tagFlags = cms.PSet(# tag cut; decision 1 or 0 wil be stored in the output tree
      HLTL2Mu3 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()"),
      HLTL2Mu7 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()"),
      HLTL2Mu15= cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty()"),
        ),
        pairVariables = cms.PSet(
        pt  = cms.string("pt"),
        y = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
        ),  
    pairFlags = cms.PSet(),
    isMC = cms.bool(False),
        #tagMatches = cms.InputTag("tagMuonsMCMatch"),
        #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
        #motherPdgId = cms.int32(23),
        #makeMCUnbiasTree = cms.bool(True),
        #checkMotherInUnbiasEff = cms.bool(True),
        allProbes     = cms.InputTag("probeMuonsSta"),
    )


##############
# Muon ID efficiency (or muId*trigger): make the fit tree and save it in the "ID" directory
process.MuonID = cms.EDAnalyzer("TagProbeFitTreeProducer",
            tagProbePairs = cms.InputTag("tpPairsMuIdNew"),
            arbitration   = cms.string("OneProbe"),
            variables = cms.PSet(KinematicVariables,
                                TrackQualityVariables,
                                GlobalTrackQualityVariables,
                                L1Variables,
                                IP = cms.string('abs(dB(\"PV3D\"))'),
                                ),
            flags = cms.PSet(MuonIDFlags,
                            HLTL1v0 = cms.string("!triggerObjectMatchesByPath(\'HLT_HIL1DoubleMu0_HighQ_v*\').empty()"),
                            HLTL1v1 = cms.string("!triggerObjectMatchesByFilter(\'hltHIDoubleMuLevel1PathL1HighQFiltered\').empty()"),
                            HLTL1v2 = cms.string("(!triggerObjectMatchesByPath(\'HLT_HIL1DoubleMu0_HighQ_v*\').empty() && !triggerObjectMatchesByFilter(\'hltHIDoubleMuLevel1PathL1HighQFiltered\').empty())"),
                            PassingGlb = cms.string("isTrackerMuon && (!triggerObjectMatchesByPath(\'HLT_HIL1DoubleMu0_HighQ_v*\').empty() && !triggerObjectMatchesByFilter(\'hltHIDoubleMuLevel1PathL1HighQFiltered\').empty()) && " + QUALITY_CUTS),
        ),
        tagVariables = cms.PSet( # tag variables that will be stored in the output tree
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)")
       ),
    tagFlags = cms.PSet(# tag cut; decision 1 or 0 wil be stored in the output tree
      HLTL2Mu3 = cms.string("!triggerObjectMatchesByPath(\'HLT_HIL2Mu3_NHitQ_v*\').empty() && !triggerObjectMatchesByFilter(\'hltHIL2Mu3NHitL2Filtered\').empty()"),
      HLTL2Mu7 = cms.string("!triggerObjectMatchesByPath(\'HLT_HIL2Mu7_v*\').empty() && !triggerObjectMatchesByFilter(\'hltHIL2Mu7L2Filtered\').empty()"),
      HLTL2Mu15= cms.string("!triggerObjectMatchesByPath(\'HLT_HIL2Mu15_v*\').empty() && !triggerObjectMatchesByFilter(\'hltHIL2Mu15L2Filtered\').empty()"),
        ),
        
        pairVariables = cms.PSet(
            pt  = cms.string("pt"),
            y = cms.string("rapidity"),
            absy = cms.string("abs(rapidity)"),
            ),  
        pairFlags = cms.PSet(),
        isMC = cms.bool(False),
        #tagMatches = cms.InputTag("tagMuonsMCMatch"),
        #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
        #motherPdgId = cms.int32(23),
        #makeMCUnbiasTree = cms.bool(True),
        #checkMotherInUnbiasEff = cms.bool(True),
        allProbes     = cms.InputTag("probeMuonsTrk"),
    )




process.tnpSimpleSequence = cms.Sequence(
   process.tagMuonsSglTrgNew *
   process.tpPairsTrigNew *
   process.tpPairsMuIdNew *
   process.tpPairsStaNew *
   process.MuonTrg *
   process.MuonTrk *
   process.MuonID
)

process.tagAndProbe = cms.Path(
    process.centralityFilter *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_Regit_PbPb_All.root"))

