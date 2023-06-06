### HiForest Configuration
# Collisions: PbPb
# Type: Data
# Input: AOD

import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('HiForest',eras.Run2_2018_pp_on_AA)

###############################################################################
# HiForest labelling info
###############################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest 103X")
import subprocess, os
version = subprocess.check_output(['git',
    '-C', os.path.expandvars('$CMSSW_BASE/src'), 'describe', '--tags'])
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

###############################################################################
# Input source
###############################################################################

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        "file:28E4F661-1674-2E42-B8BF-6C9F643E185A.root"
),
    )

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

###############################################################################
# Load Global Tag, Geometry, etc.
###############################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_v6', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

print('\n\033[31m~*~ USING CENTRALITY TABLE FOR PbPb 2018 DATA ~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = True
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcZDChits = cms.InputTag("QWzdcreco")
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","reRECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")

###############################################################################
# Define tree output
###############################################################################

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestAOD.root"))

###############################################################################
# Additional Reconstruction and Analysis: Main Body
###############################################################################

# Make PAT Muons
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")

# Add trigger matching
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import switchOffAmbiguityResolution, addHLTL1Passthrough, useL1Stage2Candidates
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
# Use L1 Stage2 (Run2 2016-onwards)
addHLTL1Passthrough(process)
useL1Stage2Candidates(process)
if "hltL1extraParticles" in process.patTrigger.collections:
    process.patTrigger.collections.remove("hltL1extraParticles")
if "hltGtStage2Digis:Muon" not in process.patTrigger.collections:
    process.patTrigger.collections.append("hltGtStage2Digis:Muon")
process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGtStage2Digis:Muon")')
process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
process.muonMatchHLTL1.preselection = cms.string("")
# Modifications for L2 and L3 HLT muons in PbPb 2018s
if "hltL2MuonCandidates" in process.patTrigger.collections:
    process.patTrigger.collections.remove("hltL2MuonCandidates")
if "hltL2MuonCandidatesPPOnAA" not in process.patTrigger.collections:
    process.patTrigger.collections.append("hltL2MuonCandidatesPPOnAA")
process.muonMatchHLTL2.matchedCuts = cms.string('coll("hltL2MuonCandidatesPPOnAA")')
if "hltIterL3MuonCandidates" in process.patTrigger.collections:
    process.patTrigger.collections.remove("hltIterL3MuonCandidates")
if "hltIterL3MuonCandidatesPPOnAA" not in process.patTrigger.collections:
    process.patTrigger.collections.append("hltIterL3MuonCandidatesPPOnAA")
process.muonMatchHLTL3.matchedCuts = cms.string('coll("hltIterL3MuonCandidatesPPOnAA")')
# Selection used to match online-offline muons
process.muonL1Info.maxDeltaR = cms.double(0.3)
process.muonL1Info.maxDeltaEta   = cms.double(0.2)
process.muonL1Info.fallbackToME1 = cms.bool(True)
process.muonMatchHLTL1.maxDeltaR = cms.double(0.3)
process.muonMatchHLTL1.maxDeltaEta   = cms.double(0.2)
process.muonMatchHLTL1.fallbackToME1 = cms.bool(True)
process.muonMatchHLTL2.maxDeltaR = cms.double(0.3)
process.muonMatchHLTL2.maxDPtRel = cms.double(10.0)
process.muonMatchHLTL3.maxDeltaR = cms.double(0.1)
process.muonMatchHLTL3.maxDPtRel = cms.double(10.0)

process.muonTree = cms.EDAnalyzer("MuonTree",
    muons = cms.InputTag("patMuonsWithTrigger"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    genparticles = cms.InputTag("genParticles"),
)

###############################################################################

############################
# Event Analysis
############################
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')

###############################################################################
#Recover peripheral primary vertices
#https://twiki.cern.ch/twiki/bin/view/CMS/HITracking2018PbPb#Peripheral%20Vertex%20Recovery
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")

# ZDC RecHit Producer
process.load('RecoHI.ZDCRecHit.QWZDC2018Producer_cfi')
process.load('RecoHI.ZDCRecHit.QWZDC2018RecHit_cfi')

###############################################################################

############################
# Event Selection
############################

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Double muon triggers
    'HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v*', # Peripheral OS dimuons
    'HLT_HIL1DoubleMuOpen_Centrality_50_100_v*', # Peripheral dimuons
    'HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v*', # Bottomonia
    'HLT_HIL1DoubleMu10_v*', # Z bosons
    'HLT_HIUPC_DoubleMu0_NotMBHF2AND_v*', # UPC dimuons
    # Single muon triggers
    'HLT_HIL1MuOpen_Centrality_80_100_v*', # Peripheral muons
    'HLT_HIL3Mu12_v*', # Electroweak bosons
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*', # UPC muons
]

# Add muon selection
process.twoMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"), minNumber = cms.uint32(2))

process.evtSel = cms.Sequence(process.hltFilter * process.twoMuons)
process.eventSelection = cms.Path( process.evtSel )

###############################################################################

#########################
# Main analysis list
#########################

process.ana_step = cms.Path(
    process.evtSel +
    process.offlinePrimaryVerticesRecovery +
    process.hltanalysis +
    process.zdcdigi +
    process.QWzdcreco +
    process.hiCentrality +
    process.centralityBin +
    process.hiEvtAnalyzer +
    process.patMuonsWithTriggerSequence +
    process.muonTree
    )

###############################################################################

#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.evtSel * process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.evtSel * process.primaryVertexFilter)
process.pBeamScrapingFilter = cms.Path(process.evtSel * process.beamScrapingFilter)
process.collisionEventSelectionAOD = cms.Path(process.evtSel * process.collisionEventSelectionAOD)
process.collisionEventSelectionAODv2 = cms.Path(process.evtSel * process.collisionEventSelectionAODv2)

process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
process.phfCoincFilter2Th4 = cms.Path(process.evtSel * process.hfCoincFilter2Th4)
process.phfPosFilterTh3 = cms.Path(process.evtSel * process.hfPosFilterTh3_seq)
process.phfPosFilterTh4 = cms.Path(process.evtSel * process.hfPosFilterTh4_seq)
process.phfPosFilterTh5 = cms.Path(process.evtSel * process.hfPosFilterTh5_seq)
process.phfPosFilterTh6 = cms.Path(process.evtSel * process.hfPosFilterTh6_seq)
process.phfPosFilterTh7 = cms.Path(process.evtSel * process.hfPosFilterTh7_seq)
process.phfPosFilterTh8 = cms.Path(process.evtSel * process.hfPosFilterTh8_seq)
process.phfPosFilterTh7p3 = cms.Path(process.evtSel * process.hfPosFilterTh7p3_seq)
process.phfNegFilterTh3 = cms.Path(process.evtSel * process.hfNegFilterTh3_seq)
process.phfNegFilterTh4 = cms.Path(process.evtSel * process.hfNegFilterTh4_seq)
process.phfNegFilterTh5 = cms.Path(process.evtSel * process.hfNegFilterTh5_seq)
process.phfNegFilterTh6 = cms.Path(process.evtSel * process.hfNegFilterTh6_seq)
process.phfNegFilterTh7 = cms.Path(process.evtSel * process.hfNegFilterTh7_seq)
process.phfNegFilterTh8 = cms.Path(process.evtSel * process.hfNegFilterTh8_seq)
process.phfNegFilterTh7p6 = cms.Path(process.evtSel * process.hfNegFilterTh7p6_seq)

process.skimanalysis.superFilters = cms.vstring("eventSelection")
process.pAna = cms.EndPath(process.skimanalysis)

from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"

###############################################################################

# Customization
