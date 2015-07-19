import FWCore.ParameterSet.Config as cms

process = cms.Process("TnPTrg")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_44_V10::50100'

process.HeavyIonGlobalParameters = cms.PSet(
        centralityVariable = cms.string("HFtowers"),
        nonDefaultGlauberModel = cms.string(""),
        centralitySrc = cms.InputTag("hiCentrality")
)

process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
#process.centralityFilter.selectedBins = [0,1,2,3] # 0 - 10 %
#process.centralityFilter.selectedBins = [4,5,6,7] # 10 - 20 %
#process.centralityFilter.selectedBins = [8,9,10,11,12,13,14,15,16,17,18,19] # 20 - 50 %
process.centralityFilter.selectedBins = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] # 50 - 100 %
#process.centralityFilter.selectedBins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] # MinBias

process.source = cms.Source("PoolSource", 
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
        "tnp_Regit_Skim.root",
    ),
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
)    

TRACK_CUTS = "track.numberOfValidHits > 10 && track.normalizedChi2 < 4 && track.hitPattern.pixelLayersWithMeasurement > 0"
GLB_CUTS = "isGlobalMuon && globalTrack.normalizedChi2 < 20  && abs(dB) < 3 && abs(track.dz) < 15 && muonID('TrackerMuonArbitrated')"#move id cut to tracking efficiency
QUALITY_CUTS =  GLB_CUTS + ' && ' + TRACK_CUTS
IN_ACCEPTANCE = '( (abs(eta)<1.0 && pt>=3.4) || (1.0<=abs(eta)<1.5 && pt>=5.8-2.4*abs(eta)) || (1.5<=abs(eta)<2.4 && pt>=3.3667-7.0/9.0*abs(eta)) )'

process.tagMuonsSglTrg2 = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("tagMuonsSglTrg"),
        cut = cms.string("(!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu3NHitL2Filtered').empty()) || (!triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu7L2Filtered').empty()) || (!triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty() && !triggerObjectMatchesByFilter('hltHIL2Mu15L2Filtered').empty())"),
)

process.tagMuonsDblTrg2 = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("tagMuonsDblTrg"),
        cut = cms.string("(!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty())"),
)

process.tpPairsTrigNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuons@-')
)

process.tpPairsTrigNew2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuons@-')
)

process.tpPairsTrigNew3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 3.0 && pt < 30.0 && abs(y) > 1.6 && abs(y) < 2.4'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuons@-')
)

process.tpPairsMuIdNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0'),
        decay = cms.string('tagMuonsDblTrg2@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNew2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsDblTrg2@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNew3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 3.0 && pt < 30.0 && abs(y) > 1.6 && abs(y) < 2.4'),
        decay = cms.string('tagMuonsDblTrg2@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNewSgl = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNewSgl2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNewSgl3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 3.0 && pt < 30.0 && abs(y) > 1.6 && abs(y) < 2.4'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsTrk@-')
)


process.tpPairsStaNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 '),
        decay = cms.string('tagMuonsDblTrg2@+ probeMuonsSta@-')
)

process.tpPairsStaNewSgl = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 '),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsSta@-')
)

process.tpPairsStaNewSgl2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 && pt > 6.5 && pt < 30.0 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsSta@-')
)

process.tpPairsStaNewSgl3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 && pt > 3.0 && pt < 30.0 && abs(y) > 1.6 && abs(y) < 2.4'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsSta@-')
)

process.probeMuonsSta2 = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("probeMuonsSta"),
        cut = cms.string(IN_ACCEPTANCE),
)

process.tpPairsStaNew2Sgl = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 '),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsSta2@-')
)

process.tpPairsStaNew2Sgl2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 && pt > 6.5 && pt < 30.0 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsSta2@-')
)

process.tpPairsStaNew2Sgl3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 && pt > 3.0 && pt < 30.0 && abs(y) > 1.6 && abs(y) < 2.4'),
        decay = cms.string('tagMuonsSglTrg2@+ probeMuonsSta2@-')
)



# Make the fit tree and save it in the "MuonID" directory
process.MuonTrg = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsTrigNew"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      HLTL1v0 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
      HLTL1v1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL1v2 = cms.string("(!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty())"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      HLTL1 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      GLB = cms.string("isGlobalMuon"),
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

process.MuonTrg2 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsTrigNew2"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      HLTL1v0 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
      HLTL1v1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL1v2 = cms.string("(!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty())"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      HLTL1 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      GLB = cms.string("isGlobalMuon"),
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

process.MuonTrg3 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsTrigNew3"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      HLTL1v0 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
      HLTL1v1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL1v2 = cms.string("(!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty() && !triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty())"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      HLTL1 = cms.string("!triggerObjectMatchesByPath('HLT_HIL1DoubleMu0_HighQ_v*').empty()"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      GLB = cms.string("isGlobalMuon"),
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


process.MuonTrk = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsSta"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
        isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta"),
)

process.MuonTrkSgl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNewSgl"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta"),
)

process.MuonTrkSgl2 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNewSgl2"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta"),
)

process.MuonTrkSgl3 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNewSgl3"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta"),
)

process.MuonTrk2Sgl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNew2Sgl"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta2"),
)

process.MuonTrk2Sgl2 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNew2Sgl2"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta2"),
)

process.MuonTrk2Sgl3 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNew2Sgl3"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta2"),
)




process.MuonTrkDbl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsStaNew"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      isGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatchSta"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsSta"),
)

process.MuonID = cms.EDAnalyzer("TagProbeFitTreeProducer",
    #tagProbePairs = cms.InputTag("tpPairsMuId"),
    tagProbePairs = cms.InputTag("tpPairsTracks"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.MuonIDNew = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIdNew"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.MuonIDNew2 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIdNew2"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.MuonIDNew3 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIdNew3"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.MuonIDSgl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIdNewSgl"),
    #tagProbePairs = cms.InputTag("tpPairsTracks"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.MuonIDSgl2 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIdNewSgl2"),
    #tagProbePairs = cms.InputTag("tpPairsTracks"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.MuonIDSgl3 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairsMuIdNewSgl3"),
    #tagProbePairs = cms.InputTag("tpPairsTracks"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      phi = cms.string("phi"),
    ),
    flags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
    ),
    tagFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
      #PassingGlb = cms.InputTag("tkPassingGlb"),
    ),
    pairVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
    ),  
    pairFlags = cms.PSet(
      PassingGlb = cms.string(QUALITY_CUTS),
    ),
    isMC = cms.bool(False),
    #tagMatches = cms.InputTag("tagMuonsMCMatch"),
    #probeMatches  = cms.InputTag("probeMuonsMCMatch"),
    #motherPdgId = cms.int32(23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    allProbes     = cms.InputTag("probeMuonsTrk"),
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuonsSglTrg2 *
    process.tagMuonsDblTrg2 *
    process.probeMuonsSta2 *
    process.tpPairsTrigNew * # tagMuonsSglTrg2@+ probeMuons@-, 6.5 pt cut
    process.tpPairsTrigNew2 * # tagMuonsSglTrg2@+ probeMuons@-, 6.5 pt cut, y < 1.6
    process.tpPairsTrigNew3 * # tagMuonsSglTrg2@+ probeMuons@-, 6.5 pt cut, 1.6 < y < 2.4
    process.tpPairsMuIdNew * # tagMuonsDblTrg2@+ probeMuonsTrk@-
    process.tpPairsMuIdNew2 * # tagMuonsDblTrg2@+ probeMuonsTrk@-, pt > 6.5 && abs(y) < 1.6
    process.tpPairsMuIdNew3 * # tagMuonsDblTrg2@+ probeMuonsTrk@-, pt > 6.5 && abs(y) > 1.6 && abs(y) < 2.4
    process.tpPairsMuIdNewSgl * # tagMuonsSglTrg2@+ probeMuonsTrk@-, pt > 6.5
    process.tpPairsMuIdNewSgl2 * # tagMuonsSglTrg2@+ probeMuonsTrk@-, pt > 6.5 && abs(y) < 1.6
    process.tpPairsMuIdNewSgl3 * # tagMuonsSglTrg2@+ probeMuonsTrk@-, pt > 6.5 && abs(y) > 1.6 && abs(y) < 2.4
    process.tpPairsStaNew * # tagMuonsDblTrg2@+ probeMuonsSta@-
    process.tpPairsStaNewSgl * # tagMuonsSglTrg2@+ probeMuonsSta@-
    process.tpPairsStaNewSgl2 * # tagMuonsSglTrg2@+ probeMuonsSta@-
    process.tpPairsStaNewSgl3 * # tagMuonsSglTrg2@+ probeMuonsSta@-
    process.tpPairsStaNew2Sgl *
    process.tpPairsStaNew2Sgl2 *
    process.tpPairsStaNew2Sgl3 *
    process.MuonTrg *
    process.MuonTrg2 *
    process.MuonTrg3 *
    process.MuonTrk *
    process.MuonTrkSgl *
    process.MuonTrkSgl2 *
    process.MuonTrkSgl3 *
    process.MuonTrk2Sgl *
    process.MuonTrk2Sgl2 *
    process.MuonTrk2Sgl3 *
    process.MuonTrkDbl *
    process.MuonID *
    process.MuonIDNew *
    process.MuonIDNew2 *
    process.MuonIDNew3 *
    process.MuonIDSgl * 
    process.MuonIDSgl2 * 
    process.MuonIDSgl3 
)

process.tagAndProbe = cms.Path(
    process.centralityFilter *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_Regit_PbPb_50100.root"))

