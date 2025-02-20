### HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process('HiForest', Run3_pp_on_PbPb_2023)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, mc")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        '/store/mc/HINPbPbSpring23MiniAOD/TT_TuneCP5_5p36TeV_powheg-pythia8/MINIAODSIM/132X_mcRun3_2023_realistic_HI_v9-v3/2810000/003e1086-4b39-4210-aecf-2335e661fe92.root'
    ),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_HI_v10', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_MC103X_2018PbPb_v4"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
         )
])

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

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
process.hltobject.triggerNames = trigger_list_data_2023_skimmed

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.ggHiNtuplizer.genParticleSrc = "prunedGenParticles"
process.ggHiNtuplizer.doPackedGenParticle = False
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
process.load('HeavyIonsAnalysis.MuonAnalysis.hiIsoMuons_cfi')
process.hiIsoMuons.muon_minPt = 10
process.unpackedMuons.muons = "hiIsoMuons"
process.muonSequence = cms.Sequence(process.hiIsoMuons * process.unpackedMuons)
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.muonAnalyzer.doGen = cms.bool(True)
###############################################################################

# ZDC RecHit Producer
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')

process.zdcdigi.SOI = cms.untracked.int32(2)
process.zdcanalyzer.doZDCRecHit = False
process.zdcanalyzer.doZDCDigi = True
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.zdcDigiSrc = cms.InputTag("hcalDigis", "ZDC")
process.zdcanalyzer.calZDCDigi = False
process.zdcanalyzer.verbose = False
process.zdcanalyzer.nZdcTs = cms.int32(6)

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
    process.ggHiNtuplizer +
    process.hiFJSoftKillerAnalyzer +
    process.rhoFlowMCSequence +
    process.metFilters +
    process.zdcanalyzer
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
for jetR in [0.3, 0.4]:
    R = str(int(jetR*10))
    from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettings_cff import candidateBtaggingMiniAOD
    candidateBtaggingMiniAOD(process, isMC = True, jetPtMin = jetPtMin, jetR = jetR, jetCorrLevels = ['L2Relative', 'L3Absolute'])

    # setup jet analyzer
    setattr(process,f'akCs{R}PFJetAnalyzer', process.akCs4PFJetAnalyzer.clone())
    getattr(process,f'akCs{R}PFJetAnalyzer').genjetTag = f'ak{R}GenJetsRecluster'
    getattr(process,f'akCs{R}PFJetAnalyzer').jetTag = f'selectedUpdatedPatJetsAKCs{R}DeepFlavour'
    getattr(process,f'akCs{R}PFJetAnalyzer').jetName = f'akCs{R}PF'
    getattr(process,f'akCs{R}PFJetAnalyzer').rParam = jetR
    getattr(process,f'akCs{R}PFJetAnalyzer').matchJets = matchJets
    getattr(process,f'akCs{R}PFJetAnalyzer').matchTag = f'patJetsAK{R}PFUnsubJets'
    getattr(process,f'akCs{R}PFJetAnalyzer').unsubjet_map = cms.untracked.InputTag(f"unsubAK{R}JetMap")
    getattr(process,f'akCs{R}PFJetAnalyzer').doHiJetID = doHIJetID
    getattr(process,f'akCs{R}PFJetAnalyzer').doWTARecluster = doWTARecluster
    getattr(process,f'akCs{R}PFJetAnalyzer').useNewBtaggers = True
    getattr(process,f'akCs{R}PFJetAnalyzer').jetPtMin = jetPtMin
    getattr(process,f'akCs{R}PFJetAnalyzer').useRawPt = True
    getattr(process,f'akCs{R}PFJetAnalyzer').jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
    getattr(process,f'akCs{R}PFJetAnalyzer').pfJetProbabilityBJetTag = cms.untracked.string(f"pfJetProbabilityBJetTagsAKCs{R}DeepFlavour")
    getattr(process,f'akCs{R}PFJetAnalyzer').pfDeepCSVJetTags = cms.untracked.string(f"pfDeepCSVJetTagsAKCs{R}DeepFlavour")
    getattr(process,f'akCs{R}PFJetAnalyzer').pfDeepFlavourJetTags = cms.untracked.string(f"pfDeepFlavourJetTagsAKCs{R}DeepFlavour")
    getattr(process,f'akCs{R}PFJetAnalyzer').pfParticleTransformerAK4JetTags = cms.untracked.string(f"pfParticleTransformerAK4JetTagsAKCs{R}DeepFlavour")
    getattr(process,f'akCs{R}PFJetAnalyzer').pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string(f"pfUnifiedParticleTransformerAK4JetTagsAKCs{R}DeepFlavour")
    process.forest += getattr(process,f'akCs{R}PFJetAnalyzer')


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


###############################################################################
# quench particles
from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process, new="hiQuenchedParticles:packedPFCandidates", old="packedPFCandidates")
MassReplaceInputTag(process, new="hiQuenchedParticles:lostTracks", old="lostTracks")
MassReplaceInputTag(process, new="hiQuenchedParticles:lostTrackseleTracks", old="lostTracks:eleTracks")

process.unpackedTracksAndVertices.packedCandidates = cms.VInputTag("packedPFCandidates", "lostTracks", "lostTracks:eleTracks")
process.load('HeavyIonsAnalysis.EventAnalysis.hiQuenchedParticles_cfi')
'''
process.PackedPFTowersUnquenched = process.PackedPFTowers.clone(src = cms.InputTag("packedPFCandidates"))
process.hiPuRhoUnquenched = process.hiPuRho.clone(src = cms.InputTag("PackedPFTowersUnquenched"))
process.ak4PFJetsForFlowUnquenched = process.ak4PFJetsForFlow.clone(src = cms.InputTag("PackedPFTowersUnquenched"))
process.hiFJRhoFlowModulationUnquenched = process.hiFJRhoFlowModulation.clone(
    jetTag = cms.InputTag("ak4PFJetsForFlowUnquenched"),
    pfCandSource = cms.InputTag("packedPFCandidates")
)
process.akCs4PFUnquenchedJets = process.akCs4PFJets.clone(
    etaMap = cms.InputTag("hiPuRhoUnquenched","mapEtaEdges"),
    rhom = cms.InputTag("hiPuRhoUnquenched","mapToRhoM"),
    rho = cms.InputTag("hiPuRhoUnquenched","mapToRho"),
    rhoFlowFitParams = cms.InputTag("hiFJRhoFlowModulationUnquenched","rhoFlowFitParams"),
    src = cms.InputTag("packedPFCandidates"),
    jetPtMin = 15.
)
process.quenchTask = cms.Task(process.centralityBin, process.PackedPFTowersUnquenched, process.hiPuRhoUnquenched, process.ak4PFJetsForFlowUnquenched, process.hiFJRhoFlowModulationUnquenched, process.akCs4PFUnquenchedJets, process.hiQuenchedParticles)
'''
process.ak4GenJetsWithNuUnquenched = process.ak4GenJetsWithNu.clone(jetPtMin = 15.)
process.hiQuenchedParticles.jets = "ak4GenJetsWithNuUnquenched"
process.quenchTask = cms.Task(process.centralityBin, process.ak4GenJetsWithNuUnquenched, process.hiQuenchedParticles)
process.forest.associate(process.quenchTask)
