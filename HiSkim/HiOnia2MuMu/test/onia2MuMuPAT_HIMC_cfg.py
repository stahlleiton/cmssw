# for the list of used tags please see:
# https://twiki.cern.ch/twiki/bin/view/CMS/Onia2MuMuSamples

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# set up process
process = cms.Process("Onia2MuMuPAT")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.inputFiles = 'file:/tmp/tdahms/6AEAB3A5-1A21-E211-A508-D4AE5264CC75.root'
options.outputFile = 'onia2MuMuPAT_MC.root'

options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'STARTHI44_V12::All'


# event plane stuff
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hievtplaneflatproducer_cfi")
process.load('RecoHI.HiCentralityAlgos.CentralityBin_cfi')

# produce missing l1extraParticles
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco_step = cms.Path(process.l1extraParticles)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# BSC or HF coincidence (masked unprescaled L1 bits)
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.bscOrHfCoinc = process.l1Filter.clone(
    algorithms = cms.vstring('L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1', 'L1_NotBsc2_BscMinBiasOR', 'L1_HcalHfCoincidencePm')
    )

# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltOniaHI = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltOniaHI.HLTPaths = ["HLT_HIL1DoubleMu0_HighQ_v*",
                              "HLT_HIL2Mu3_NHitQ_v*",
                              "HLT_HIL2Mu7_v*","HLT_HIL2Mu15_v*",
                              "HLT_HIL2DoubleMu0_NHitQ_v*",
                              "HLT_HIL2DoubleMu3_v*",
                              "HLT_HIL3Mu3_v*",
                              "HLT_HIL3DoubleMuOpen_v*","HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v*"
                              ]
process.hltOniaHI.throw = False
process.hltOniaHI.andOr = True
process.hltOniaHI.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")


from HiSkim.HiOnia2MuMu.onia2MuMuPAT_cff import *
onia2MuMuPAT(process, GlobalTag=process.GlobalTag.globaltag, MC=True, HLT="HLT", Filter=False)

process.onia2MuMuPatGlbGlb.dimuonSelection  = cms.string("mass > 0")
process.onia2MuMuPatGlbGlb.addMuonlessPrimaryVertex = False
process.onia2MuMuPatGlbGlb.resolvePileUpAmbiguity = False

process.source.fileNames = cms.untracked.vstring(
    options.inputFiles
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.outOnia2MuMu.fileName = cms.untracked.string( options.outputFile )
#process.outTnP.fileName = cms.untracked.string( options.outputFile )

# add event plane flattening
'''
process.load("CondCore.DBCommon.CondDBCommon_cfi");
process.CondDBCommon.connect = "sqlite_file:flatparms_PbPb_2011.db"
process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
                                      process.CondDBCommon,
                                      toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                 tag = cms.string('flatParamtest')
                                                                 )
                                                        )
                                      )
'''
# add event plane information
#process.ProdEvtPlane = cms.Path(process.centralityBin * process.hiEvtPlane * process.hiEvtPlaneFlat)
process.ProdEvtPlane = cms.Path(process.centralityBin)

process.outOnia2MuMu.outputCommands.extend(cms.untracked.vstring('keep *_generator_*_*'))
process.outTnP.outputCommands.extend(cms.untracked.vstring('keep *_generator_*_*'))

process.e = cms.EndPath(process.outOnia2MuMu)# + process.outTnP)


process.patMuonSequence = cms.Sequence(
    #process.bscOrHfCoinc *
    process.hltOniaHI *
    process.collisionEventSelection *
    process.genMuons *
    process.patMuonsWithTriggerSequence
    )

process.schedule = cms.Schedule(process.L1Reco_step,process.ProdEvtPlane,process.Onia2MuMuPAT,
                                #process.TagAndProbeSta, process.TagAndProbeMuID, process.TagAndProbeTrig,
                                process.e)
