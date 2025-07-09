# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --python_filename Run3_2025_OXY_RECO.py --conditions 150X_dataRun3_Prompt_v3 -s RAW2DIGI,L1Reco,RECO --datatier AOD --eventcontent AOD --data --process RECO --scenario pp --era Run3_2025_OXY --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3_2025_OXY --nThreads 8 -n -1 --filein /store/data/pORun2025/IonPhysics0/RAW/v1/000/393/974/00000/00027414-031a-473e-8182-24cc3c7b4c21.root --fileout file:output.root --customise_commands process.dedxHitInfo = process.dedxAllHitInfo.clone()
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2025_UPC_OXY_cff import Run3_2025_UPC_OXY

process = cms.Process('RECO',Run3_2025_UPC_OXY)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/pORun2025/IonPhysics0/RAW/v1/000/393/974/00000/00027414-031a-473e-8182-24cc3c7b4c21.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--python_filename nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.load('RecoTracker.DeDxCalibration.collisionEventSelection_cff')
process.hpTracks = cms.EDFilter("TrackSelector", src = cms.InputTag("generalTracks"), cut = cms.string("quality(\"highPurity\") && abs(eta)<3.0"))
process.twoHPTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("hpTracks"), minNumber = cms.uint32(2))
process.eventSel = cms.Sequence(process.clusterCompatibilityFilter * process.hpTracks * process.twoHPTracks * process.primaryVertexFilter)

from SimGeneral.MixingModule.SiStripSimParameters_cfi import SiStripSimBlock as _SiStripSimBlock
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import siPixelClusters as _siPixelClusters

process.microDstProducer = cms.EDAnalyzer("MicroDstProducer",
    MeVPerElectron = cms.double(1000*_SiStripSimBlock.GevPerElectron.value()),
    VCaltoElectronGain = _siPixelClusters.VCaltoElectronGain,
    VCaltoElectronGain_L1 = _siPixelClusters.VCaltoElectronGain_L1,
    VCaltoElectronOffset = _siPixelClusters.VCaltoElectronOffset,
    VCaltoElectronOffset_L1 = _siPixelClusters.VCaltoElectronOffset_L1,
    pixelSaturationThr = cms.int32(254),
    trackProducer = cms.InputTag('generalTracks'),
    dedxHitInfo   = cms.InputTag('dedxHitInfo'),
    dedxMomentum = cms.InputTag('dedxHitInfo:momentumAtHit')
)

process.trackAnalyzer = cms.EDAnalyzer('TrackAnalyzer',
    doTrack = cms.untracked.bool(True),
    trackPtMin = cms.untracked.double(0.01),
    vertexSrc = cms.InputTag("offlinePrimaryVertices"),
    trackSrc = cms.InputTag("generalTracks"),
    beamSpotSrc = cms.untracked.InputTag('offlineBeamSpot'),
    dedxEstimators = cms.VInputTag(["dedxAllLikelihood", "dedxPixelLikelihood", "dedxStripLikelihood", "dedxHarmonic2", "dedxPixelHarmonic2"])
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hadronTree.root')
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_Prompt_v3', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.eventSel * process.L1Reco)
process.reconstruction_step = cms.Path(process.eventSel * process.reconstruction)
process.produceMicroDst = cms.Path(process.eventSel * process.microDstProducer * process.trackAnalyzer)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.produceMicroDst)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 4
process.options.numberOfStreams = 0

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run3_2025_OXY 

#call to customisation function customisePostEra_Run3_2025_OXY imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run3_2025_OXY(process)

# End of customisation functions


# Customisation from command line

process.dedxHitInfo = process.dedxAllHitInfo.clone(storeMomentumAtHit = True)
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
