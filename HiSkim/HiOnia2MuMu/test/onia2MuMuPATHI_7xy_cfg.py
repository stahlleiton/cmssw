# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# set up process
process = cms.Process("Onia2MuMuPAT")

# Conditions
isPbPb = True;          
isMC = False;
keepGeneralTracks = True;
keepEventPlane = False;
muonSelection = "GlbTrk" # Single muon selection: Glb, GlbTrk, Trk are availale

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.inputFiles = 'file:/afs/cern.ch/user/t/tuos/work/public/reco2AOD/round2April28/DIMUON/step2_RAW2DIGI_L1Reco_DIMUONskim_AOD.root'#step2_reRECO_740_100_1_lRV.root'
#/afs/cern.ch/user/t/tuos/work/public/reco2AOD/round2April28/DIMUON/step2_RAW2DIGI_L1Reco_DIMUONskim_AOD.root'#step2_reRECO_740_100_1_lRV.root'
options.outputFile = 'onia2MuMuPAT_740.root'

options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")


# Event plane (Not working currently)
#process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
#process.ProdEvtPlane = cms.Path(process.hiEvtPlane)

# tag for running on 2011 data in 7xy
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if isMC:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_HIon', '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')

'''
# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('*','L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
    )
'''

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths = [
    "*",
    "HLT_HIL1DoubleMu0_v*",
    "HLT_HIL1DoubleMu0_2HF_v*",
    "HLT_HIL1DoubleMu0_2HF0_v*",
    "HLT_HIL1DoubleMu10_v*",
    "HLT_HIL2DoubleMu0_NHitQ_v*",
    "HLT_HIL2DoubleMu0_NHitQ_2HF_v*",
    "HLT_HIL2DoubleMu0_NHitQ_2HF0_v*",
    "HLT_HIL3DoubleMu0_2HF_v*",
    "HLT_HIL3DoubleMu0_OS_m2p5to4p5_v*",
    "HLT_HIL3DoubleMu0_OS_m7to14_v*",
    "HLT_HIL1DoubleMu0_2HF_Cent30100_v*",
    "HLT_HIL1DoubleMu0_2HF0_Cent30100_v*",
    "HLT_HIL2DoubleMu0_2HF_Cent30100_NHitQ_v*",
    "HLT_HIL2DoubleMu0_2HF0_Cent30100_NHitQ_v*",
    "HLT_HIL1DoubleMu0_Cent30_v*",
    "HLT_HIL2DoubleMu0_Cent30_NHitQ_v*",
    "HLT_HIL3DoubleMu0_Cent30_v*",
    "HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v*",
    "HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v*",
    "HLT_HIL2Mu3_NHitQ10_v*",
    "HLT_HIL2Mu3_NHitQ10_2HF_v*",
    "HLT_HIL2Mu3_NHitQ10_2HF0_v*",
    "HLT_HIL3Mu3_NHitQ15_2HF_v*",
    "HLT_HIL3Mu3_NHitQ15_2HF0_v*",
    "HLT_HIL2Mu15_NHitQ10_v*",
    "HLT_HIL2Mu15_NHitQ10_2HF_v*",
    "HLT_HIL2Mu15_NHitQ10_2HF0_v*",
    "HLT_HIL3Mu15_NHitQ10_v*",
    "HLT_HIL3Mu15_NHitQ15_2HF_v*",
    "HLT_HIL3Mu15_NHitQ15_2HF0_v*"
                             ]

process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True
process.hltOniaHI.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=isMC, HLT="HLT", Filter=False)

##### Onia2MuMuPAT input collections/options
process.onia2MuMuPatGlbGlb.dimuonSelection          = cms.string("mass > 0")
process.onia2MuMuPatGlbGlb.resolvePileUpAmbiguity   = False
if isPbPb:
  process.onia2MuMuPatGlbGlb.primaryVertexTag         = cms.InputTag("hiSelectedVertex")
  process.patMuonsWithoutTrigger.pvSrc                = cms.InputTag("hiSelectedVertex")
  process.onia2MuMuPatGlbGlb.addMuonlessPrimaryVertex = False
  if isMC:
    process.genMuons.src = "genParticles"
    process.onia2MuMuPatGlbGlb.genParticles = "genParticles"
else: # ispp
  process.onia2MuMuPatGlbGlb.primaryVertexTag         = cms.InputTag("offlinePrimaryVertices")
  process.patMuonsWithoutTrigger.pvSrc                = cms.InputTag("offlinePrimaryVertices")
  # Adding muonLessPV gives you lifetime values wrt. muonLessPV only
  process.onia2MuMuPatGlbGlb.addMuonlessPrimaryVertex = True
  if isMC:
    process.genMuons.src = "genParticles"
    process.onia2MuMuPatGlbGlb.genParticles = "genParticles"

##### Remove few paths for MC
if isMC:
  process.patMuonSequence.remove(process.hltOniaHI)
if not isMC:
  process.patMuonSequence.remove(process.hltOniaHI)

##### Dimuon pair selection
commonP1 = "|| (innerTrack.isNonnull && genParticleRef(0).isNonnull)"
commonP2 = " && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35"
if muonSelection == "Glb":
  highP = "isGlobalMuon"; # At least one muon must pass this selection
  process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("("+highP+commonP1+")"+commonP2)
  lowP = "isGlobalMuon"; # BOTH muons must pass this selection
  process.onia2MuMuPatGlbGlb.lowerPuritySelection = cms.string("("+lowP+commonP1+")"+commonP2)
elif muonSelection == "GlbTrk":
  highP = "isGlobalMuon";
  process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("("+highP+commonP1+")"+commonP2)
  lowP = "(isGlobalMuon && isTrackerMuon)";
  process.onia2MuMuPatGlbGlb.lowerPuritySelection = cms.string("("+lowP+commonP1+")"+commonP2)
elif muonSelection == "Trk":
  highP = "isTrackerMuon";
  process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("("+highP+commonP1+")"+commonP2)
  lowP = "isTrackerMuon";
  process.onia2MuMuPatGlbGlb.lowerPuritySelection = cms.string("("+lowP+commonP1+")"+commonP2)
else:
  print "Using default settings in HiSkim/HiOnia2MuMu/python/onia2MuMuPAT_cff.py file!"

##### If single track collection has to be kept
if keepGeneralTracks:
  process.outOnia2MuMu.outputCommands.append("keep *_hiGeneralTracks_*_*")

##### If event plane collection has to be kept
if keepEventPlane:
  process.outOnia2MuMu.outputCommands.extend(('keep *_hiEvtPlane_*_*','keep *_hiEvtPlaneFlat_*_*'))



process.source.fileNames      = cms.untracked.vstring(options.inputFiles)
process.maxEvents             = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.outOnia2MuMu.fileName = cms.untracked.string( options.outputFile )
process.e                     = cms.EndPath(process.outOnia2MuMu)
process.schedule              = cms.Schedule(process.Onia2MuMuPAT,process.e)

from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process)
