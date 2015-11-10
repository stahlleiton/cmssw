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
options.outputFile = "Jpsi_TEST.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = 'file:/home/llr/cms/chapon/data_CMS/promptskims2015/CMSSW_7_5_4/test/step2_reRECO_740_9_1_e5r.root'
#step2_reRECO_740_100_1_lRV.root'
#/afs/cern.ch/user/t/tuos/work/public/reco2AOD/round2April28/DIMUON/step2_RAW2DIGI_L1Reco_DIMUONskim_AOD.root'#step2_reRECO_740_100_1_lRV.root'
#options.maxEvents = 10 # -1 means all events

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


if isPbPb:
  process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
  process.centralityBin.Centrality = cms.InputTag("hiCentrality")
  if isMC:
    process.centralityBin.centralityVariable = cms.string("HFtowersHydjetDrum5")
    process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")  
    process.GlobalTag.toGet.extend([cms.PSet(
          record  = cms.string("HeavyIonRcd"),
          tag     = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v750x02_mc"),
          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
          label   = cms.untracked.string("HFtowersHydjetDrum5")
          )])
  else:
    process.centralityBin.centralityVariable = cms.string("HFtowers")
    process.centralityBin.nonDefaultGlauberModel = cms.string("")    
    process.GlobalTag.toGet.extend([cms.PSet(
          record  = cms.string("HeavyIonRcd"),
          tag     = cms.string("CentralityTable_HFtowers200_Glauber2010A_eff99_run1v750x01_offline"),
          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
          label   = cms.untracked.string("HFtowers")
          )])
'''
# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths = [
    "*",
    "HLT_HIL1DoubleMu0_v*",
    "HLT_HIL1DoubleMu0_2HF_v*",
    "HLT_HIL1DoubleMu0_2HF0_v*",
                             ]
'''

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


process.hionia = cms.EDAnalyzer('HiOniaAnalyzer',
                                #-- Collections
                                srcMuon             = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig       = cms.InputTag("patMuonsWithoutTrigger"),
                                src                 = cms.InputTag("onia2MuMuPatGlbGlb"),
                                srcTracks           = cms.InputTag("hiGeneralTracks"),
                                EvtPlane            = cms.InputTag("hiEvtPlane","recoLevel"),
                                EvtPlaneFlat        = cms.InputTag("hiEvtPlaneFlat",""),

                                triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                                useRapidity = cms.bool(True),
                                
                                #--
                                maxAbsZ = cms.double(24.0),
                                
                                pTBinRanges      = cms.vdouble(0.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 40.0),
                                etaBinRanges     = cms.vdouble(0.0, 2.5),
                                centralityRanges = cms.vdouble(20,40,100),

                                onlyTheBest     = cms.bool(False),		
                                applyCuts       = cms.bool(True),
                                storeEfficiency = cms.bool(False),
                      
                                removeSignalEvents = cms.untracked.bool(False),
                                removeTrueMuons    = cms.untracked.bool(False),
                                storeSameSign      = cms.untracked.bool(True),
                                
                                #-- Gen Details
                                oniaPDG = cms.int32(443),
                                muonSel = cms.string(muonSelection),
                                isHI = cms.untracked.bool(isPbPb),
                                isPA = cms.untracked.bool(False),
                                isMC = cms.untracked.bool(isMC),
                                isPromptMC = cms.untracked.bool(False),
                                useEvtPlane = cms.untracked.bool(keepEventPlane),
                                useGeTracks = cms.untracked.bool(keepGeneralTracks),
                                runVersionChange = cms.untracked.uint32(182133),

                                #-- Histogram configuration
                                combineCategories = cms.bool(False),
                                fillRooDataSet    = cms.bool(False),
                                fillTree          = cms.bool(True),
                                fillHistos        = cms.bool(False),
                                minimumFlag       = cms.bool(True),
                                fillSingleMuons   = cms.bool(True),
                                fillRecoTracks    = cms.bool(False),
                                histFileName      = cms.string(options.outputFile),		
                                dataSetName       = cms.string(options.secondaryOutputFile),
                                
                                #--
                                # NumberOfTriggers = cms.uint32(8),
                                dblTriggerPathNames = cms.vstring("HLT_HIL1DoubleMu0_HighQ_v2",
                                                                  "HLT_HIL2DoubleMu3_v2",
                                                                  "HLT_HIL3DoubleMuOpen_v2",
                                                                  "HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2"),
                                dblTriggerFilterNames = cms.vstring("hltHIDoubleMuLevel1PathL1HighQFiltered",
                                                                    "hltHIL2DoubleMu3L2Filtered",
                                                                    "hltHIDimuonL3FilteredOpen",
                                                                    "hltHIDimuonL3FilteredMg2OSnoCowboy"),
                                sglTriggerPathNames = cms.vstring("HLT_HIL2Mu3_NHitQ_v2",
                                                                  "HLT_HIL2Mu7_v2",
                                                                  "HLT_HIL2Mu15_v2",
                                                                  "HLT_HIL3Mu3_v2"),
                                sglTriggerFilterNames = cms.vstring("hltHIL2Mu3NHitL2Filtered",
                                                                    "hltHIL2Mu7L2Filtered",
                                                                    "hltHIL2Mu15L2Filtered",
                                                                    "hltHISingleMu3L3Filtered")
                                )

if isPbPb:
  process.hionia.primaryVertexTag = cms.InputTag("hiSelectedVertex")
  process.hionia.genParticles     = cms.InputTag("genParticles")
  process.hionia.muonLessPV       = cms.bool(False)
  process.hionia.CentralitySrc    = cms.InputTag("hiCentrality")
  if isMC:
    process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowersHydjetDrum5")  
  else:
    process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
  
else:    
  process.hionia.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")
  process.hionia.genParticles     = cms.InputTag("genParticles")
  process.hionia.muonLessPV       = cms.bool(True)
  process.hionia.CentralitySrc    = cms.InputTag("")
  process.hionia.CentralityBinSrc = cms.InputTag("")


##### Remove few paths for MC
if isMC:
  if isPbPb:
    process.Onia2MuMuPAT = cms.Path(
       process.genMuons *
       process.patMuonsWithTriggerSequence *
       process.patMuonSequence *
       process.onia2MuMuPatGlbGlb *
       process.onia2MuMuPatGlbGlbFilter *
       process.centralityBin *
       process.hionia
       )
  else:
    process.Onia2MuMuPAT = cms.Path(
       process.genMuons *
       process.patMuonsWithTriggerSequence *
       process.patMuonSequence *
       process.onia2MuMuPatGlbGlb *
       process.onia2MuMuPatGlbGlbFilter *
       process.hionia
       )
    
if not isMC:
  if isPbPb:
    process.Onia2MuMuPAT = cms.Path(
       process.patMuonsWithTriggerSequence *
       process.patMuonSequence *
       process.onia2MuMuPatGlbGlb *
       process.onia2MuMuPatGlbGlbFilter *
       process.centralityBin *
       process.hionia
       )
  else:
    process.Onia2MuMuPAT = cms.Path(
       process.patMuonsWithTriggerSequence *
       process.patMuonSequence *
       process.onia2MuMuPatGlbGlb *
       process.onia2MuMuPatGlbGlbFilter *
       process.hionia
       )

process.schedule              = cms.Schedule(process.Onia2MuMuPAT)

from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process)
