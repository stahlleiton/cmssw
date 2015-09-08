import FWCore.ParameterSet.Config as cms

process = cms.Process("TnPTrg")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'STARTHI44_V12::All'


# centrality stuff
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')

process.HeavyIonGlobalParameters = cms.PSet(
        centralitySrc = cms.InputTag("hiCentrality"),
        centralityVariable = cms.string("HFtowers"),
        #centralityVariable = cms.string("HFhits"),
        nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
        #nonDefaultGlauberModel = cms.string("Hydjet_Bass"), 
    )


process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
process.centralityFilter.selectedBins = [16,17,18,19] # 40 - 50 %


process.source = cms.Source("PoolSource", 
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
       "file:tnp_regit_100_1_8G0.root",
    ),
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
)    

# Selection : acceptance and muId, for innerTrack and globalTrack
IN_ACCEPTANCE = '( (abs(eta)<1.0 && pt>=3.4) || (1.0<=abs(eta)<1.5 && pt>=5.8-2.4*abs(eta)) || (1.5<=abs(eta)<2.4 && pt>=3.3667-7.0/9.0*abs(eta)) )'
# tracke muon cuts 
TRACK_CUTS    = "isTrackerMuon && innerTrack.numberOfValidHits > 10 && innerTrack.normalizedChi2 < 4 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && muonID('TrackerMuonArbitrated')"
# global muon cuts
GLB_CUTS      = "isGlobalMuon && globalTrack.normalizedChi2 < 20"
QUALITY_CUTS  =  "(" + GLB_CUTS + ' && ' + TRACK_CUTS + ")"

# Old cuts Pre-HP 2015
TRACK_CUTS2 = "isTrackerMuon && track.numberOfValidHits > 10 && track.normalizedChi2 < 4 && track.hitPattern.pixelLayersWithMeasurement > 0"
GLB_CUTS2 = "isGlobalMuon && globalTrack.normalizedChi2 < 20  && abs(dB) < 3 && abs(track.dz) < 15 && muonID('TrackerMuonArbitrated')"#move id cut to tracking efficiency
QUALITY_CUTS2 =  GLB_CUTS2 + ' && ' + TRACK_CUTS2
DXYZ_CUTS = "abs(dB) < 3 && abs(track.dz) < 15"



muonIDFlags = cms.PSet(
    GlobalMu  = cms.string("isGlobalMuon"),
    TrackerMu = cms.string("isTrackerMuon"),
    TMA    = cms.string("muonID('TrackerMuonArbitrated')"),
    TMLSAT = cms.string("muonID('TMLastStationAngTight')"),
    TrackCuts	= cms.string(TRACK_CUTS),
    GlobalCuts	= cms.string(GLB_CUTS),
    QualityCuts	= cms.string(QUALITY_CUTS),
    TrackCuts2	= cms.string(TRACK_CUTS2),
    GlobalCuts2	= cms.string(GLB_CUTS2),
    QualityCuts2 = cms.string(QUALITY_CUTS2),
    dxyzCuts = cms.string(DXYZ_CUTS),
)

staOnlyVariables = cms.PSet(
    staQoverP      = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverp"),
    staQoverPerror = cms.string("? outerTrack.isNull() ? 0 : outerTrack.qoverpError"),
    staValidStations = cms.string("? outerTrack.isNull() ? -1 : outerTrack.hitPattern.muonStationsWithValidHits()"),
    staNumValidHits = cms.string("? outerTrack.isNull() ? -1 : outerTrack.hitPattern.numberOfValidMuonHits()"),
)

muonIDVariables = cms.PSet(
    globalMuVar  = cms.string("isGlobalMuon"),
    trackerMuVar = cms.string("isTrackerMuon"),
    tmaVar    = cms.string("muonID('TrackerMuonArbitrated')"),
    tmlsatVar = cms.string("muonID('TMLastStationAngTight')"),
)

TrigTagFlags = cms.PSet(
    HLTL2Mu3 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()"),
    HLTL2Mu7 = cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()"),
    HLTL2Mu15= cms.string("!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty()"),
)


TrigProbeFlags = cms.PSet(
       HLTL1v0 = cms.string("!triggerObjectMatchesByPath(\'HLT_HIL1DoubleMu0_HighQ_v*\').empty()"),
       HLTL1v1 = cms.string("!triggerObjectMatchesByFilter(\'hltHIDoubleMuLevel1PathL1HighQFiltered\').empty()"),
       HLTL1v2 = cms.string("(!triggerObjectMatchesByPath(\'HLT_HIL1DoubleMu0_HighQ_v*\').empty() && !triggerObjectMatchesByFilter(\'hltHIDoubleMuLevel1PathL1HighQFiltered\').empty())"),
)

####new probes
process.probeMuonsTrg = cms.EDFilter("PATMuonSelector",
                src = cms.InputTag("probeMuonsTrk"),
                cut = cms.string(QUALITY_CUTS),
)

process.muonDxyPVdzMinTrg = cms.EDProducer("MuonDxyPVdzmin",
    probes = cms.InputTag("probeMuonsTrg"),
)
process.muonDxyPVdzMinIDTrg = cms.EDProducer("MuonDxyPVdzmin",
    probes = cms.InputTag("probeMuonsTrk"),
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
        decay = cms.string('tagMuonsSglTrgNew@+ probeMuonsTrg@-')
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
    arbitration   = cms.string("OneProbe"), # have unique tag-probe for each event
    variables = cms.PSet( # probe variables that will be stored in the output tree
    KinematicVariables,
    L1Variables,
    dxyPVdzmin       = cms.InputTag("muonDxyPVdzMinTrg","dxyPVdzmin"),
    dzPV       = cms.InputTag("muonDxyPVdzMinTrg","dzPV"),
    absdB = cms.string("abs(dB)"),
    absdz = cms.string("abs(track.dz)"),
 
    ),
    flags = cms.PSet(TrigProbeFlags, # passing probe condition (the efficiency we want to probe); decision 1 or 0 wil be stored in the output tree
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
    tagFlags = cms.PSet(TrigTagFlags,# tag cut; decision 1 or 0 wil be stored in the output tree
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
    allProbes     = cms.InputTag("probeMuonsTrg"),
    eventWeight = cms.double(160.944),
    #addRunLumiInfo = cms.bool(True),
)





##############
# Tracking efficiency by itself: make the fit tree and save it in the "Tracking" directory
process.MuonTrk = cms.EDAnalyzer("TagProbeFitTreeProducer",
        tagProbePairs = cms.InputTag("tpPairsStaNew"),
        arbitration   = cms.string("OneProbe"),
        variables = cms.PSet(
        KinematicVariables, 
        staOnlyVariables,
        absdB = cms.string("? innerTrack.isNull() ? -1 : abs(dB)"),
        absdz = cms.string("? innerTrack.isNull() ? -1 : abs(track.dz)"),

    ),
    flags = cms.PSet(
        muonIDFlags,
        outerValidHits  = cms.string("outerTrack.numberOfValidHits > 0"),
        Tk              = cms.string("track.isNonnull"),
        StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
       	PassingSta = cms.string("isGlobalMuon &&" + GLB_CUTS),
    ),
    tagVariables = cms.PSet( # tag variables that will be stored in the output tree
        pt  = cms.string("pt"),
        eta = cms.string("eta"),
        abseta = cms.string("abs(eta)")
       ),
    tagFlags = cms.PSet(# tag cut; decision 1 or 0 wil be stored in the output tree
        TrigTagFlags,
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
	eventWeight = cms.double(160.944),
    )


##############
# Muon ID efficiency (or  muId*trigger): make the fit tree and save it in the "ID" directory
process.MuonIDTrg = cms.EDAnalyzer("TagProbeFitTreeProducer",
            tagProbePairs = cms.InputTag("tpPairsMuIdNew"),
            arbitration   = cms.string("OneProbe"),
            variables = cms.PSet(KinematicVariables,
                                TrackQualityVariables,
                                GlobalTrackQualityVariables,
                                L1Variables,
                                muonIDVariables,			
                                dxyPVdzmin       = cms.InputTag("muonDxyPVdzMinIDTrg","dxyPVdzmin"),
                                dzPV       = cms.InputTag("muonDxyPVdzMinIDTrg","dzPV"),
                                absdB = cms.string("abs(dB)"),
                                absdz = cms.string("abs(track.dz)"),

                                
                                ),
            flags = cms.PSet(muonIDFlags,
                             TrigProbeFlags,
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
                       TrigTagFlags,
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
	eventWeight = cms.double(160.944),
    )




process.tnpSimpleSequence = cms.Sequence(process.probeMuonsTrg *
                                         process.muonDxyPVdzMinTrg *
                                         process.muonDxyPVdzMinIDTrg *
                                         process.tagMuonsSglTrgNew *
	                                 process.tpPairsTrigNew *
                                         process.tpPairsMuIdNew *
                                         process.tpPairsStaNew *
                                         process.MuonTrk *
                                         process.MuonIDTrg *
                                         process.MuonTrg
)

process.tagAndProbe = cms.Path(
    process.centralityFilter *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_Regit_PbPb_All.root"))

