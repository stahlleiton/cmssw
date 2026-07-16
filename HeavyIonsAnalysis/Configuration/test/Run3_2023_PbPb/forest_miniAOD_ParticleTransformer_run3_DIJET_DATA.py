### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process('HiForest', Run3_pp_on_PbPb_2023)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/375/790/00000/56ad580f-b228-4f3c-b8e3-17f9d95c7654.root'),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# Define centrality binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

###############################################################################

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')
process.metFilters = process.skimanalysis.clone(hltresults = "TriggerResults::RECO")

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
process.hltobject.triggerNames = trigger_list_data_2023_skimmed
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import dummy_branches_for_PbPb_2023_HLT
process.hltanalysis.hltdummybranches = dummy_branches_for_PbPb_2023_HLT

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.egammaSequence = cms.Sequence(process.ggHiNtuplizer)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJRhoAnalyzer_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJSoftKillerAnalyzer_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFlowRhoAnalyzer_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
# muons
process.load("HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi")
process.load('HeavyIonsAnalysis.MuonAnalysis.hiMuons_cfi')
process.hiMuons.muon_minPt = 10
process.hiMuons.file_isoModel = "HeavyIonsAnalysis/MuonAnalysis/data/muiso_BDT.root"
process.hiMuons.file_isoCorr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2023_PbPb.json.gz"
process.hiMuons.era = "Run3_2023_PbPb"
process.unpackedMuons.muons = "hiMuons"
process.muonSequence = cms.Sequence(process.hiMuons * process.unpackedMuons)
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
###############################################################################

#########################
# ZDC RecHit Producer && Analyzer
#########################
process.load('HeavyIonsAnalysis.ZDCAnalysis.ZDCAnalyzersHC2023_cff')
process.zdcSequencePbPb = process.zdcSequence.copy()

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
    process.rhoSequence +
    process.muonSequence +
    process.egammaSequence +
    process.hiFJSoftKillerAnalyzer +
    process.rhoFlowDataSequence +
    process.metFilters +
    process.zdcSequencePbPb
    )

#customisation
process.particleFlowAnalyser.ptMin = 0.0
process.ggHiNtuplizer.muonPtMin = 0.0

#########################
# Apply egamma regression
#########################
from HeavyIonsAnalysis.EGMAnalysis.applyEgammaRegression_cfi import applyEgammaRegression
process = applyEgammaRegression(process, era = "Run3_2023_PbPb")

# Select the types of jets filled
jetPtMin = 15
jetAbsEtaMax = 2.5

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = True        # Add jet phi and eta for WTA axis

# add candidate tagging
for jetR, doFlow in zip([0.3], [False]):
    R = str(int(jetR*10))
    from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettings_cff import candidateBtaggingMiniAOD
    candidateBtaggingMiniAOD(process, isMC = False, jetPtMin = jetPtMin, jetR = jetR, jetCorrLevels = ['L2Relative', 'L2L3Residual'], doFlow = doFlow, addNegTag = False, era = "Run3_2023_PbPb")

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
    getattr(process,f'ak{jL}PFJetAnalyzer').pfDeepCSVJetTags = cms.untracked.string(f"pfDeepCSVJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfDeepFlavourJetTags = cms.untracked.string(f"pfDeepFlavourJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfParticleTransformerAK4JetTags = cms.untracked.string(f"pfParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string(f"pfUnifiedParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    process.forest += getattr(process,f'ak{jL}PFJetAnalyzer')


#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.HiZDCfilter_cfi')
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
process.oneLeptonOrTwoJets1 = cms.EDFilter("PATCountFilter",
    electronSource = cms.InputTag("goodElectrons"),
    muonSource     = cms.InputTag("goodMuons"),
    tauSource      = cms.InputTag(""),
    photonSource   = cms.untracked.InputTag("goodPhotons"),
    jetSource      = cms.untracked.InputTag("akCs3PFJets"),
    jetTagSource   = cms.untracked.InputTag(""),
    jetMinPt       = cms.untracked.double(15),
    countElectrons = cms.bool(True),
    countMuons     = cms.bool(True),
    countTaus      = cms.bool(False),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000000),
    minJets = cms.untracked.uint32(2),
)
process.oneLeptonOrTwoJets2 = cms.EDFilter("PATCountFilter",
    electronSource = cms.InputTag("goodElectrons"),
    muonSource     = cms.InputTag("goodMuons"),
    tauSource      = cms.InputTag(""),
    photonSource   = cms.untracked.InputTag("goodPhotons"),
    jetSource      = cms.untracked.InputTag("akCs3PFJets"),
    jetTagSource   = cms.untracked.InputTag("pfUnifiedParticleTransformerAK4JetTagsAKCs3DeepFlavour:ptcorr"),
    jetMinPt       = cms.untracked.double(30),
    countElectrons = cms.bool(True),
    countMuons     = cms.bool(True),
    countTaus      = cms.bool(False),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000000),
    minJets = cms.untracked.uint32(2),
)
process.leptonJetSelection = cms.Sequence(process.goodElectrons * process.goodMuons * process.goodPhotons * process.PackedPFTowers * process.hiPuRho * process.oneLeptonOrTwoJets1 * process.oneLeptonOrTwoJets2)
process.filterSequence = cms.Sequence(
    process.primaryVertexFilter *
    process.leptonJetSelection
)

process.superFilterPath = cms.Path(process.filterSequence)
process.skimanalysis.superFilters = cms.vstring("superFilterPath")

for path in process.paths:
    if path != "superFilterPath":
        getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq
