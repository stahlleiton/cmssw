### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2025_OXY_cff import Run3_2025_OXY
process = cms.Process('HiForest', Run3_2025_OXY)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 150X, data")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/hidata/OORun2025/IonPhysics10/MINIAOD/PromptReco-v1/000/394/217/00000/e1cb21d9-45cf-46f8-b4f6-b1f9e27bc0da.root'),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_Prompt_v1', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# Define centrality binning
#process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin = cms.EDProducer('HICentralityBinProducer')
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("PFhf")
process.centralityBin.table = cms.vdouble(
    0,0.91197,1.82394,2.73591,3.64788,4.55985,5.47182,6.38379,7.29576,8.20773,9.1197,10.0317,10.6667,10.7669,11.1084,11.4511,11.7967,12.1453,12.4985,12.8586,13.2279,13.6047,13.9872,14.378,14.7734,15.1775,15.5875,16.0048,16.4344,16.8729,17.3249,17.7878,18.2572,18.7342,19.2173,19.7139,20.221,20.7368,21.2654,21.7955,22.3333,22.8819,23.446,24.0097,24.5881,25.179,25.7856,26.3914,27.0037,27.6247,28.2619,28.9063,29.5473,30.2005,30.8754,31.5574,32.2518,32.9558,33.6616,34.3925,35.1279,35.8711,36.618,37.3905,38.1709,38.9568,39.76,40.5757,41.3886,42.2135,43.0613,43.921,44.8002,45.6699,46.5569,47.4536,48.3673,49.2971,50.2437,51.1962,52.1656,53.1554,54.1508,55.1587,56.1668,57.1942,58.2508,59.3168,60.406,61.5058,62.6151,63.7513,64.8938,66.0488,67.2111,68.4189,69.6168,70.8422,72.0903,73.3433,74.6371,75.9344,77.2465,78.5776,79.9278,81.3019,82.7002,84.105,85.5389,86.9959,88.4836,89.9969,91.5244,93.0648,94.6373,96.2243,97.8198,99.454,101.105,102.798,104.466,106.186,107.925,109.692,111.487,113.265,115.092,116.936,118.856,120.819,122.768,124.744,126.732,128.739,130.809,132.907,135.033,137.173,139.34,141.543,143.777,146.028,148.298,150.674,153.07,155.464,157.912,160.369,162.846,165.358,167.907,170.503,173.094,175.76,178.495,181.231,183.991,186.784,189.674,192.575,195.507,198.523,201.592,204.653,207.796,210.983,214.245,217.554,220.881,224.256,227.741,231.249,234.786,238.444,242.096,245.842,249.627,253.523,257.577,261.697,265.88,270.16,274.601,279.162,283.888,288.745,293.767,298.968,304.526,310.308,316.408,322.893,329.956,337.457,345.746,354.995,365.685,378.512,395.16,420.836,526.238
)
from RecoHI.HiCentralityAlgos.CentralityBin_cfi import centralityBin
process.centralityBin2024PbPb = centralityBin.clone()
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(
        record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
    ),
])

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

###############################################################################

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.hiEvtAnalyzer.CentralityBinSrc = "centralityBin:PFhf"
process.hiEvtAnalyzer.doHFfilters = False
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')
process.metFilters = process.skimanalysis.clone(hltresults = "TriggerResults::RECO")

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2025_OO_skimmed
process.hltobject.triggerNames = trigger_list_data_2025_OO_skimmed
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import dummy_branches_for_OO_2025_HLT
process.hltanalysis.hltdummybranches = dummy_branches_for_OO_2025_HLT

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.load('HeavyIonsAnalysis.EGMAnalysis.hiElectrons_cfi')
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedPatElectronProducer_cfi')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJRhoAnalyzer_cff')
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedPatPhotonProducer_cfi')
process.correctedPhotons = process.correctedPatPhotonProducer.clone(src = "slimmedPhotons", centrality = "centralityBin:PFhf", minPt = 15)
process.correctedPhotons.correctionFile = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2025_pO/SuperClusterSS_pO2025_DATA.dat"
from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier
process.correctedElectrons = process.correctedPatElectronProducer.clone(src = "slimmedElectrons", centrality = "centralityBin:PFhf", epCombConfig = regressionModifier.eleRegs.epComb, minPt = 15, calibrateSuperCluster = False)
process.correctedElectrons.correctionFile = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2025_pO/ElectronSS_pO2025_DATA.dat"
process.hiElectrons.centrality = "centralityBin2024PbPb:HFtowers"
process.hiElectrons.electrons = "correctedElectrons"
process.hiElectrons.file_idModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleid_BDT.ubj"
process.hiElectrons.file_isoModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleiso_BDT.ubj"
process.hiElectrons.file_corr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2024_PbPb.json.gz"
process.hiElectrons.era = "Run3_2024_PbPb"
process.ggHiNtuplizer.muonSrc = "slimmedMuons"
process.ggHiNtuplizer.useValMapIso = False
process.ggHiNtuplizer.photonSrc = "correctedPhotons"
process.ggHiNtuplizer.electronSrc = "hiElectrons"
process.egammaSequence = cms.Sequence(process.centralityBin2024PbPb * process.rhoSequence * process.correctedElectrons * process.correctedPhotons * process.hiElectrons * process.ggHiNtuplizer)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
# muons
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.muonAnalyzer.muonSrc = "slimmedMuons"
###############################################################################

#########################
# ZDC RecHit Producer && Analyzer
#########################
# to prevent crash related to HcalSeverityLevelComputerRcd record
process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
process.load('HeavyIonsAnalysis.ZDCAnalysis.ZDCAnalyzersPbPb_cff')
process.zdcanalyzer.doZdcDigis = False

###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.centralityBin +
    process.hiEvtAnalyzer +
    process.hltanalysis +
    process.hltobject +
    process.l1object +
    process.unpackedTracksAndVertices +
    process.particleFlowAnalyser +
    process.egammaSequence +
    process.metFilters +
    process.zdcSequencePbPb
    )

#customisation
process.particleFlowAnalyser.ptMin = 0.0
process.ggHiNtuplizer.muonPtMin = 0.0

# Select the types of jets filled
jetPtMin = 15
jetAbsEtaMax = 2.5

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = True        # Add jet phi and eta for WTA axis

# add candidate tagging
for jetR, doFlow in zip([0.4], [False]):
    R = str(int(jetR*10))
    from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettings_cff import candidateBtaggingMiniAOD
    candidateBtaggingMiniAOD(process, isMC = False, jetPtMin = jetPtMin, jetR = jetR, jetCorrLevels = ['L2Relative', 'L2L3Residual'], doFlow = doFlow, addNegTag = True, era = "Run3_2025_PbPb")

    # setup jet analyzer
    jL = f"Cs{R}Flow" if doFlow else f"Cs{R}"
    setattr(process,f'ak{jL}PFJetAnalyzer', process.akCs4PFJetAnalyzer.clone())
    getattr(process,f'ak{jL}PFJetAnalyzer').jetTag = f'selectedUpdatedPatJetsAK{jL}DeepFlavour'
    getattr(process,f'ak{jL}PFJetAnalyzer').jetName = f'ak{jL}PF'
    getattr(process,f'ak{jL}PFJetAnalyzer').rParam = jetR
    getattr(process,f'ak{jL}PFJetAnalyzer').doHiJetID = doHIJetID
    getattr(process,f'ak{jL}PFJetAnalyzer').doWTARecluster = doWTARecluster
    getattr(process,f'ak{jL}PFJetAnalyzer').jetPtMin = jetPtMin
    getattr(process,f'ak{jL}PFJetAnalyzer').useRawPt = True
    getattr(process,f'ak{jL}PFJetAnalyzer').jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
    getattr(process,f'ak{jL}PFJetAnalyzer').pfJetProbabilityBJetTag = cms.untracked.string(f"pfJetProbabilityBJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfNegativeOnlyJetProbabilityBJetTag = cms.untracked.string(f"pfNegativeOnlyJetProbabilityBJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfDeepCSVJetTags = cms.untracked.string(f"pfDeepCSVJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfDeepFlavourJetTags = cms.untracked.string(f"pfDeepFlavourJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfParticleTransformerAK4JetTags = cms.untracked.string(f"pfParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string(f"pfUnifiedParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfNegativeUnifiedParticleTransformerAK4JetTags = cms.untracked.string(f"pfNegativeUnifiedParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    process.forest += getattr(process,f'ak{jL}PFJetAnalyzer')


#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilterPF_cfi')
process.pAna = cms.EndPath(process.skimanalysis)

process.goodMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons::@skipCurrentProcess"),
    cut = cms.string("pt >= 15.0 && passed('CutBasedIdLoose')")
)
process.goodElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons::@skipCurrentProcess"),
    cut = cms.string("pt >= 15.0")
)
process.goodPhotons = cms.EDFilter("PATPhotonSelector",
    src = cms.InputTag("slimmedPhotons::@skipCurrentProcess"),
    cut = cms.string("pt >= 25.0")
)
process.oneLepton = cms.EDFilter("PATCountFilter",
    electronSource = cms.InputTag("goodElectrons"),
    muonSource     = cms.InputTag("goodMuons"),
    tauSource      = cms.InputTag(""),
    photonSource   = cms.untracked.InputTag("goodPhotons"),
    countElectrons = cms.bool(True),
    countMuons     = cms.bool(True),
    countTaus      = cms.bool(False),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000000),
)
process.leptonSelection = cms.Sequence(process.goodElectrons * process.goodMuons * process.goodPhotons * process.oneLepton)
process.filterSequence = cms.Sequence(
    process.primaryVertexFilter *
    process.leptonSelection
)

process.superFilterPath = cms.Path(process.filterSequence)
process.skimanalysis.superFilters = cms.vstring("superFilterPath")

for path in process.paths:
    if path != "superFilterPath":
        getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq
