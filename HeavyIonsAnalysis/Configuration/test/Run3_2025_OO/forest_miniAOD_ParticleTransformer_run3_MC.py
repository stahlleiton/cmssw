### HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2025_OXY_cff import Run3_2025_OXY
process = cms.Process('HiForest', Run3_2025_OXY)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 150X, mc")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/mc/HINOOSpring25MiniAOD/DYto2E_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2025_forOO_realistic_v9-v3/100000/031e2d39-1e57-41c0-b0fe-f619cec09668.root'),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun3_2025_forOO_realistic_v9', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# Define centrality binning
#process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin = cms.EDProducer('HICentralityBinProducer')
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("PFhf")
process.centralityBin.table = cms.vdouble(
    0, 4.10985, 4.77625, 5.26157, 5.68513, 6.07021, 6.42387, 6.75839, 7.07004, 7.37778, 7.68672, 7.98486, 8.28043, 8.57545, 8.86176, 9.14369, 9.42725, 9.71304, 10.0147, 10.3105, 10.6073, 10.9038, 11.2042, 11.5095, 11.8193, 12.1275, 12.446, 12.7678, 13.095, 13.4164, 13.7444, 14.0892, 14.4311, 14.784, 15.1463, 15.5155, 15.8793, 16.2598, 16.6363, 17.0308, 17.4318, 17.8363, 18.2507, 18.6695, 19.1158, 19.5598, 20.0042, 20.4468, 20.9093, 21.3787, 21.8542, 22.3545, 22.8449, 23.3585, 23.8655,
    24.3914, 24.928, 25.4831, 26.0269, 26.6125, 27.1776, 27.7521, 28.3604, 28.9708, 29.5734, 30.2092, 30.8553, 31.5125, 32.1607, 32.8557, 33.5435, 34.248, 34.9785, 35.7214, 36.4689, 37.2423, 38.0336, 38.8336, 39.651, 40.464, 41.3046, 42.1642, 43.0373, 43.9275, 44.8478, 45.784, 46.6997, 47.6568, 48.6178, 49.6245, 50.6359, 51.6727, 52.7001, 53.7518, 54.8302, 55.9117, 57.0425, 58.1931, 59.375, 60.5682, 61.7735, 62.9594, 64.1904, 65.4397, 66.7374, 68.0604, 69.3719, 70.7033, 72.0891, 73.4636,
    74.8995, 76.3443, 77.8348, 79.3369, 80.8528, 82.3963, 83.9709, 85.5451, 87.1865, 88.8456, 90.4798, 92.1981, 93.9828, 95.7577, 97.5642, 99.4073, 101.278, 103.22, 105.137, 107.099, 109.138, 111.174, 113.263, 115.385, 117.522, 119.691, 121.983, 124.244, 126.547, 128.864, 131.327, 133.831, 136.317, 138.873, 141.404, 144.05, 146.784, 149.55, 152.275, 155.033, 157.909, 160.833, 163.726, 166.732, 169.794, 172.864, 176.065, 179.353, 182.606, 185.968, 189.345, 192.925, 196.385, 199.963, 203.548,
    207.266, 211.009, 214.894, 218.785, 222.834, 226.924, 231.082, 235.36, 239.681, 244.136, 248.668, 253.374, 258.162, 263.07, 268.116, 273.287, 278.555, 284.024, 289.552, 295.42, 301.42, 307.779, 314.293, 321.209, 328.506, 336.377, 344.46, 353.27, 362.966, 373.42, 385.572, 399.896, 416.711, 439.198, 473.479, 752.978
)
from RecoHI.HiCentralityAlgos.CentralityBin_cfi import centralityBin
process.centralityBin2024PbPb = centralityBin.clone()
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(
        record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetCello_v1401x0_official_MC2024"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
    ),
])

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
process.hiEvtAnalyzer.CentralityBinSrc = "centralityBin:PFhf"
process.hiEvtAnalyzer.doHFfilters = False
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')
process.metFilters = process.skimanalysis.clone(hltresults = "TriggerResults::PAT")

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2025_OO_skimmed
process.hltobject.triggerNames = trigger_list_data_2025_OO_skimmed
from HeavyIonsAnalysis.EventAnalysis.dummybranches_cff import dummy_branches_for_OO_2025_HLT
process.hltanalysis.hltdummybranches = dummy_branches_for_OO_2025_HLT

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.particleFlowAnalyser.addInfo = True
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.load('HeavyIonsAnalysis.EGMAnalysis.hiElectrons_cfi')
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedPatElectronProducer_cfi')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJRhoAnalyzer_cff')
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedPatPhotonProducer_cfi')
process.correctedPhotons = process.correctedPatPhotonProducer.clone(src = "slimmedPhotons", centrality = "centralityBin:PFhf", minPt = 15)
process.correctedPhotons.correctionFile = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2025_pO/SuperClusterSS_pO2025_MC.dat"
from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier
process.correctedElectrons = process.correctedPatElectronProducer.clone(src = "slimmedElectrons", centrality = "centralityBin:PFhf", epCombConfig = regressionModifier.eleRegs.epComb, minPt = 15, calibrateSuperCluster = False)
process.correctedElectrons.correctionFile = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2025_pO/ElectronSS_pO2025_MC.dat"
process.hiElectrons.centrality = "centralityBin2024PbPb:HFtowers"
process.hiElectrons.electrons = "correctedElectrons"
process.hiElectrons.file_idModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleid_BDT.ubj"
process.hiElectrons.file_isoModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleiso_BDT.ubj"
process.hiElectrons.file_corr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2024_PbPb.json.gz"
process.hiElectrons.era = "Run3_2024_PbPb"
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.ggHiNtuplizer.genParticleSrc = "prunedGenParticles"
process.ggHiNtuplizer.doPackedGenParticle = False
process.ggHiNtuplizer.muonSrc = "slimmedMuons"
process.ggHiNtuplizer.useValMapIso = False
process.ggHiNtuplizer.photonSrc = "correctedPhotons"
process.ggHiNtuplizer.electronSrc = "hiElectrons"
process.egammaSequence = cms.Sequence(process.centralityBin2024PbPb * process.rhoSequence * process.correctedElectrons * process.correctedPhotons * process.hiElectrons * process.ggHiNtuplizer)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_mc_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
# muons
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.muonAnalyzer.doGen = cms.bool(True)
process.muonAnalyzer.muonSrc = "slimmedMuons"
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
    process.egammaSequence +
    process.metFilters
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
for jetR, doFlow in zip([0.4], [False]):
    R = str(int(jetR*10))
    from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettings_cff import candidateBtaggingMiniAOD
    candidateBtaggingMiniAOD(process, isMC = True, jetPtMin = jetPtMin, jetR = jetR, jetCorrLevels = ['L2Relative', 'L3Absolute'], doFlow = doFlow, addNegTag = True, era = "Run3_2025_PbPb")

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
    #mute JetFlavourClustering error
    getattr(process,f'patJetFlavourAssociationAK{jL}PF').relPtTolerance = cms.double(1.0)
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
