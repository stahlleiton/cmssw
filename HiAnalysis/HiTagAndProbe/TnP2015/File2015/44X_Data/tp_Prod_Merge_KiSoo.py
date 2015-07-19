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
process.GlobalTag.globaltag = 'GR_R_44_V10::All'

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

process.tpPairsTrigNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0'),
        decay = cms.string('tagMuonsSglTrg@+ probeMuons@-')
)

process.tpPairsTrigNew2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsSglTrg@+ probeMuons@-')
)

process.tpPairsTrigNew3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 3.0 && 1.6 < abs(y) < 2.4'),
        decay = cms.string('tagMuonsSglTrg@+ probeMuons@-')
)

process.tpPairsMuIdNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && pt < 30.0'),
        decay = cms.string('tagMuonsDblTrg@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNew2 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 6.5 && abs(y) < 1.6'),
        decay = cms.string('tagMuonsDblTrg@+ probeMuonsTrk@-')
)

process.tpPairsMuIdNew3 = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5 && pt > 3.0 && 1.6 < abs(y) < 2.4'),
        decay = cms.string('tagMuonsDblTrg@+ probeMuonsTrk@-')
)

process.tpPairsStaNew = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string(' 2 < mass < 5 '),
        decay = cms.string('tagMuonsDblTrg@+ probeMuonsSta@-')
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
      #HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      HLTL3 = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()"),
      HLTL3NoCow = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      #HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      HLTL3 = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()"),
      HLTL3NoCow = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
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
##### KiSoo start #####
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
      #HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      HLTL3 = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()"),
      HLTL3NoCow = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      #HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      HLTL3 = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()"),
      HLTL3NoCow = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
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
      #HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      HLTL3 = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()"),
      HLTL3NoCow = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
    ),
    tagVariables = cms.PSet(
      pt  = cms.string("pt"),
      eta = cms.string("eta"),
      abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
      #HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1OpenFiltered').empty()"),
      HLTL1 = cms.string("!triggerObjectMatchesByFilter('hltHIDoubleMuLevel1PathL1HighQFiltered').empty()"),
      HLTL2 = cms.string("!triggerObjectMatchesByFilter('hltHIL2DoubleMu3L2Filtered').empty()"),
      HLTL3 = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredOpen').empty()"),
      HLTL3NoCow = cms.string("!triggerObjectMatchesByFilter('hltHIDimuonL3FilteredMg2OSnoCowboy').empty()"),
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
##### KiSoo end #####
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
##### KiSoo start #####
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
##### KiSoo end #####
process.tnpSimpleSequence = cms.Sequence(
    process.tpPairsTrigNew * # tagMuonsSglTrg@+ probeMuons@-, 6.5 pt cut
##### KiSoo start #####
    process.tpPairsTrigNew2 * # tagMuonsSglTrg@+ probeMuons@-, 6.5 pt cut
    process.tpPairsTrigNew3 * # tagMuonsSglTrg@+ probeMuons@-, 6.5 pt cut
##### KiSoo end #####
    process.tpPairsMuIdNew * # tagMuonsDblTrg@+ probeMuonsTrk@-
##### KiSoo start #####
    process.tpPairsMuIdNew2 * # tagMuonsDblTrg@+ probeMuonsTrk@-
    process.tpPairsMuIdNew3 * # tagMuonsDblTrg@+ probeMuonsTrk@-
##### KiSoo end #####
    process.tpPairsStaNew * # tagMuonsDblTrg@+ probeMuonsSta@-
    process.MuonTrg *
##### KiSoo start #####
    process.MuonTrg2 *
    process.MuonTrg3 *
##### KiSoo end #####
    process.MuonTrk *
    process.MuonIDNew *
##### KiSoo start #####
    process.MuonIDNew2 *
    process.MuonIDNew3 *
##### KiSoo start #####
)

process.tagAndProbe = cms.Path(
    process.centralityFilter *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp_Regit_PbPb_All.root"))

