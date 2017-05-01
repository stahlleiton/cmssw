import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras

process = cms.Process('PAT',eras.Run2_2016_pA)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "HiChiForest.root"
options.inputFiles =  "/store/hidata/PARun2016C/PADoubleMuon/AOD/PromptReco-v1/000/286/496/00000/FE9E703F-8ABD-E611-A367-FA163EB9FCBE.root"
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_v19', '') #80X_dataRun2_Prompt_v15


process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(
    record = cms.string('EcalLaserAlphasRcd'),
    tag = cms.string("EcalLaserAlphas_v2_prompt"),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
    ),
  cms.PSet(
    record = cms.string('EcalLaserAPDPNRatiosRcd'),
    tag = cms.string("EcalLaserAPDPNRatios_prompt_v2"),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
    )
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)#,
   # wantSummary = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)

process.load('HeavyIonsAnalysis.ConversionAnalysis.conversionAnalyzer_cfi')
process.load('HeavyIonsAnalysis.ConversionAnalysis.muonAnalyzer_cfi')

process.muonAna.pfCandidatesTag = cms.InputTag("")


# Trigger Filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltL1DoubleMu = hltHighLevel.clone(
    HLTPaths =  ['HLT_PAL1DoubleMuOpen_v1'] , 
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT") , 
    throw = False , andOr = True , eventSetupPathsKey = cms.string( '' ) 
    )
process.triggerSelection = cms.Sequence( process.hltL1DoubleMu )


# HLT dimuon trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltL1DoubleMuon = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltL1DoubleMuon.HLTPaths = ["HLT_PAL1DoubleMuOpen_v1"]
process.hltL1DoubleMuon.throw = False
process.hltL1DoubleMuon.andOr = True

# selection of valid vertex
process.primaryVertexFilterForZMM = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2"), 
    filter = cms.bool(True),   # otherwise it won't filter the events
    )

# selection of dimuons (at least TRK+TRK) with mass in Onia range
process.muonSelector = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isTrackerMuon && pt > 1.0 && abs(eta) < 2.4"),
    filter = cms.bool(True)
    )

process.muonFilter = cms.EDFilter("MuonCountFilter",
    src = cms.InputTag("muonSelector"),
    minNumber = cms.uint32(2)
    )

process.dimuonMassCut = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('( 2.6 < mass < 3.6 ) || ( 9.0  < mass < 10.0 ) '),
    decay = cms.string("muonSelector@+ muonSelector@-")
    )

process.dimuonMassCutFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("dimuonMassCut"),
    minNumber = cms.uint32(1)
    )

# selection of conversions
process.conversionSelector = cms.EDFilter("ConversionSelector",
    src = cms.InputTag("allConversions"),
    cut = cms.string("(tracks().size() == 2) && conversionVertex().isValid()"),
    filter = cms.bool(True)
    )

process.conversionFilter = cms.EDFilter("ConversionCountFilter",
    src = cms.InputTag("conversionSelector"),
    minNumber = cms.uint32(1)
    )

# Gamma+Onia->mumu skim sequence
process.oniaGMMSkimSequence = cms.Sequence(
    process.hltL1DoubleMuon *
    process.primaryVertexFilterForZMM *
    process.muonSelector *
    process.muonFilter *
    process.dimuonMassCut *
    process.dimuonMassCutFilter *
    process.conversionSelector *
    process.conversionFilter 
    )

process.selectMuonEvts = cms.EDFilter("PATMuonSelector", src = cms.InputTag( 'patMuons' ) , cut = cms.string("isTrackerMuon()") )
process.selectMuonFilter = cms.EDFilter("PATCandViewCountFilter", minNumber = cms.uint32(2), maxNumber = cms.uint32(9999999), src = cms.InputTag("selectMuonEvts") )
process.muonSelection = cms.Sequence( process.selectMuonEvts + process.selectMuonFilter )

process.anaPath = cms.Path(process.oniaGMMSkimSequence + process.muonAnaSeq + process.convAnaSeq )

#Options:
process.source    = cms.Source("PoolSource", fileNames = cms.untracked.vstring( options.inputFiles ) )
process.TFileService = cms.Service("TFileService", fileName = cms.string( options.outputFile ) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# Schedule definition
process.schedule = cms.Schedule(process.anaPath,process.endjob_step)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
process.load('Configuration.StandardSequences.PAT_cff')
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)

# customisation of the process.
from PhysicsTools.PatAlgos.tools.coreTools import runOnData
runOnData( process, outputModules = [] )
