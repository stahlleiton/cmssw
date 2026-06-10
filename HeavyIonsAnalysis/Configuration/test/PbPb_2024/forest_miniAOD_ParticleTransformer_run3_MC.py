### HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2024_cff import Run3_pp_on_PbPb_2024
process = cms.Process('HiForest', Run3_pp_on_PbPb_2024)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 141X, mc")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('/store/mc/HINPbPbWinter24MiniAOD/DYto2E_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/MINIAODSIM/141X_mcRun3_2024_realistic_HI_v14-v2/120000/0015e615-1187-4d36-b2dc-af894515d657.root'),
    #fileNames = cms.untracked.vstring('file:006a8402-fe3b-4f12-93dd-a837d92f4deb.root'),
    #fileNames = cms.untracked.vstring('/store/mc/HINPbPbWinter24MiniAOD/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/MINIAODSIM/141X_mcRun3_2024_realistic_HI_v14-v2/110000/006a8402-fe3b-4f12-93dd-a837d92f4deb.root'),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun3_2024_realistic_HI_v17', '')
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

#############################
# Gen Analyzer
#############################
process.load('HeavyIonsAnalysis.EventAnalysis.HiGenAnalyzer_cfi')

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')
process.metFilters = process.skimanalysis.clone(hltresults = "TriggerResults::PAT")

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2024_skimmed
process.hltobject.triggerNames = trigger_list_data_2024_skimmed

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.load('HeavyIonsAnalysis.EGMAnalysis.hiElectrons_cfi')
process.hiElectrons.file_idModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleid_BDT.ubj"
process.hiElectrons.file_isoModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleiso_BDT.ubj"
process.hiElectrons.file_corr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2024_PbPb.json.gz"
process.hiElectrons.era = "Run3_2024_PbPb"
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.ggHiNtuplizer.genParticleSrc = "prunedGenParticles"
process.ggHiNtuplizer.doPackedGenParticle = False
process.ggHiNtuplizer.electronSrc = "hiElectrons"
process.egammaSequence = cms.Sequence(process.hiElectrons * process.ggHiNtuplizer)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_mc_cff')
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
process.hiMuons.file_isoModel = "HeavyIonsAnalysis/MuonAnalysis/data/Run3_2024_PbPb/muiso_BDT.ubj"
process.hiMuons.file_isoCorr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2024_PbPb.json.gz"
process.hiMuons.era = "Run3_2024_PbPb"
process.unpackedMuons.muons = "hiMuons"
process.muonSequence = cms.Sequence(process.hiMuons * process.unpackedMuons)
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.muonAnalyzer.doGen = cms.bool(True)
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
    process.HiGenParticleAna +
    process.rhoSequence +
    process.muonSequence +
    process.egammaSequence +
    process.hiFJSoftKillerAnalyzer +
    process.rhoFlowMCSequence +
    process.metFilters +
    process.zdcSequencePbPb
    )

#customisation
process.particleFlowAnalyser.ptMin = 0.0
process.ggHiNtuplizer.muonPtMin = 0.0

# Select the types of jets filled
matchJets = True             # Enables q/g and heavy flavor jet identification in MC
jetPtMin = 15
jetAbsEtaMax = 2.5

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = False        # Add jet phi and eta for WTA axis

# add candidate tagging
for jetR, doFlow in zip([0.3, 0.3, 0.4], [False, True, True]):
    R = str(int(jetR*10))
    from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettings_cff import candidateBtaggingMiniAOD
    candidateBtaggingMiniAOD(process, isMC = True, jetPtMin = jetPtMin, jetR = jetR, jetCorrLevels = ['L2Relative', 'L3Absolute'], doFlow = doFlow)

    # setup jet analyzer
    jL = f"Cs{R}Flow" if doFlow else f"Cs{R}"
    setattr(process,f'ak{jL}PFJetAnalyzer', process.akCs4PFJetAnalyzer.clone())
    getattr(process,f'ak{jL}PFJetAnalyzer').genjetTag = f'ak{R}GenJetsRecluster'
    getattr(process,f'ak{jL}PFJetAnalyzer').jetTag = f'selectedUpdatedPatJetsAK{jL}DeepFlavour'
    getattr(process,f'ak{jL}PFJetAnalyzer').jetName = f'ak{jL}PF'
    getattr(process,f'ak{jL}PFJetAnalyzer').rParam = jetR
    getattr(process,f'ak{jL}PFJetAnalyzer').matchJets = matchJets
    getattr(process,f'ak{jL}PFJetAnalyzer').matchTag = f'patJetsAK{R}PFUnsubJets'
    getattr(process,f'ak{jL}PFJetAnalyzer').unsubjet_map = cms.untracked.InputTag(f"unsubAK{jL}JetMap")
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
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.pphfCoincFilter4Th2 = cms.Path(process.phfCoincFilter4Th2)
process.pphfCoincFilter1Th3 = cms.Path(process.phfCoincFilter1Th3)
process.pphfCoincFilter2Th3 = cms.Path(process.phfCoincFilter2Th3)
process.pphfCoincFilter3Th3 = cms.Path(process.phfCoincFilter3Th3)
process.pphfCoincFilter4Th3 = cms.Path(process.phfCoincFilter4Th3)
process.pphfCoincFilter5Th3 = cms.Path(process.phfCoincFilter5Th3)
process.pphfCoincFilter1Th4 = cms.Path(process.phfCoincFilter1Th4)
process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
process.pphfCoincFilter3Th4 = cms.Path(process.phfCoincFilter3Th4)
process.pphfCoincFilter4Th4 = cms.Path(process.phfCoincFilter4Th4)
process.pphfCoincFilter5Th4 = cms.Path(process.phfCoincFilter5Th4)
process.pphfCoincFilter1Th5 = cms.Path(process.phfCoincFilter1Th5)
process.pphfCoincFilter2Th5 = cms.Path(process.phfCoincFilter2Th5)
process.pphfCoincFilter3Th5 = cms.Path(process.phfCoincFilter3Th5)
process.pphfCoincFilter4Th5 = cms.Path(process.phfCoincFilter4Th5)
process.pphfCoincFilter5Th5 = cms.Path(process.phfCoincFilter5Th5)
process.pphfCoincFilter1Th6 = cms.Path(process.phfCoincFilter1Th6)
process.pphfCoincFilter2Th6 = cms.Path(process.phfCoincFilter2Th6)
process.pphfCoincFilter3Th6 = cms.Path(process.phfCoincFilter3Th6)
process.pphfCoincFilter4Th6 = cms.Path(process.phfCoincFilter4Th6)
process.pphfCoincFilter5Th6 = cms.Path(process.phfCoincFilter5Th6)
process.pAna = cms.EndPath(process.skimanalysis)
