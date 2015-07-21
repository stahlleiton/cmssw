# configuration file to run in 44x on the RegIt OniaSkim, to re-make T&P skim out of the existing patMuonWithrigger matching collection

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# set up process
process = cms.Process("Onia2MuMuPATtnp")

# setup 'analysis'  options
# options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
# secondaryInputFiles = 'file:onia2MuMuPAT_regit_1000_1_yLu.root'
# inputFiles = 'file:onia2MuMuPAT_regit_1000_1_yLu.root'
# inputFiles = '/store/user/tdahms/PyquenEvtGen_jpsiMuMu_JPsiPt03/Onia2MuMu_JpsiMuMu_JpsiPt03_HiFall11_STARTHI44_V12-v0/1f8a659433db9516a729e51c7b2867f9/onia2MuMuPAT_MC_regit_100_1_K0z.root'
inputFiles = '/store/user/tdahms/PyquenEvtGen_jpsiMuMu_JPsiPt1530/Onia2MuMu_JpsiMuMu_JpsiPt1530_HiFall11_STARTHI44_V12-v0/1f8a659433db9516a729e51c7b2867f9/onia2MuMuPAT_MC_regit_100_1_bMb.root'
# inputFiles = 'file:/tmp/onia2MuMuPAT_MC_regit_100_1_K0z.root'
outputFile = 'tnp_regit.root'

maxEvents = -1 # -1 means all events

# skip events when an object is missing
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.GlobalTag.globaltag = 'GR_P_V27A::All'    # prompt reco data
process.GlobalTag.globaltag = 'STARTHI44_V12::All'   # MC
    
  # Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

    # Drop stuff on input
process.source = cms.Source("PoolSource",
      inputCommands = cms.untracked.vstring("keep *", 
         'drop *_tagMuonsDblTrgMCMatch__Onia2MuMuPAT',                # tagMuons MC matches for efficiency
         'drop *_tagMuonsSglTrgMCMatch__Onia2MuMuPAT',                # tagMuons MC matches for efficiency
         'drop patMuons_tagMuonsDblTrg__Onia2MuMuPAT',                # tagMuons for efficiency
         'drop patMuons_tagMuonsSglTrg__Onia2MuMuPAT',                # tagMuons for efficiency
         'drop patMuons_probeMuonsSta__Onia2MuMuPAT',           # probeMuons for efficiency
         'drop patMuons_probeMuonsTrk__Onia2MuMuPAT',                    # probeTracks for efficiency
         'drop patMuons_probeMuons__Onia2MuMuPAT',              # probeMuons for efficiency
         'drop *_hiEvtPlane_*_*'),
      fileNames = cms.untracked.vstring()
      )

# centrality stuff
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')

process.HeavyIonGlobalParameters = cms.PSet(
      centralityVariable = cms.string("HFtowers"),
      nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
      centralitySrc = cms.InputTag("hiCentrality")
      )

#### call the onia2MuMuPAT
from HiSkim.HiOnia2MuMu.onia2MuMuPAT_TnPonOniaSkim_cff import *
tnpOnOniaSkim(process, GlobalTag=process.GlobalTag.globaltag, MC=True, HLT="HLT", Filter=True)

process.source.fileNames = cms.untracked.vstring(inputFiles)
# process.source.secondaryFileNames = cms.untracked.vstring(secondaryInputFiles)
process.maxEvents           = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
# process.source.eventsToSkip = cms.untracked.VEventRange('182124:1559639-182124:1562415')
process.outTnP.fileName  = cms.untracked.string( outputFile )

# produce patMuons that use the STA momentum information
process.patMuonsWithTriggerSta = cms.EDProducer("RedefineMuonP4FromTrackPAT",
      src   = cms.InputTag("patMuonsWithTrigger"),
      track = cms.string("outer")
      )

process.thePatMuonsWithTriggerSta = cms.Path(process.patMuonsWithTriggerSta)

process.e             = cms.EndPath(process.outTnP)
process.schedule = cms.Schedule(
      process.thePatMuonsWithTriggerSta, 
      process.TagAndProbeTrig,
      process.TagAndProbeSta, 
      process.TagAndProbeMuID,
      process.e)

