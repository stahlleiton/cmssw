# configuration file to run in 44x on the RegIt OniaSkim, to re-make T&P skim out of the existing patMuonWithrigger matching collection

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# set up process
process = cms.Process("Onia2MuMuPATtnp")

# setup 'analysis'  options
# options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
# secondaryInputFiles = 'file:onia2MuMuPAT_regit_1000_1_yLu.root'
inputFiles = 'file:onia2MuMuPAT_regit_1000_1_yLu.root'
# inputFiles = '/store/user/tdahms/HIDiMuon/Onia2MuMu_RegIT-Skim_v3/16e111b93439e581c6bdad4acd2262ef/onia2MuMuPAT_regit_1000_1_yLu.root'
# inputFiles = 'file:tnp_regit_1000_1_fEL.root'
outputFile = 'tnp_regit.root'

maxEvents = -1 # -1 means all events

# get and parse the command line arguments

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V27A::All'    # prompt 
    
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
      nonDefaultGlauberModel = cms.string(""),
      centralitySrc = cms.InputTag("hiCentrality")
      )

#### call the onia2MuMuPAT
from HiSkim.HiOnia2MuMu.onia2MuMuPAT_TnPonOniaSkim_cff import *
tnpOnOniaSkim(process, GlobalTag=process.GlobalTag.globaltag, MC=False, HLT="HLT", Filter=True)

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

