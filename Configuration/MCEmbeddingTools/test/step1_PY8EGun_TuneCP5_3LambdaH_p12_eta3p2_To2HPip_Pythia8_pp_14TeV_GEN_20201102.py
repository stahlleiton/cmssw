# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/HypertritonP12Eta3p2_To2HPip_cfi.py --fileout=step1_PY8EGun_TuneCP5_3LambdaH_p12_eta3p2_To2HPip_Pythia8_pp_14TeV_GEN_20201102.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --step GEN,SIM --nThreads 8 --geometry Extended2023D35 --era Phase2C4_timing_layer_new --python_filename /afs/cern.ch/user/a/anstahll/work/MTD/MC/2020/CMSSW_10_4_0_mtd1/src/GEN/python/step1_PY8EGun_TuneCP5_3LambdaH_p12_eta3p2_To2HPip_Pythia8_pp_14TeV_GEN_20201102.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Phase2C4_timing_layer_new)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/HypertritonP12Eta3p2_To2HPip_cfi.py nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('step1_PY8EGun_TuneCP5_3LambdaH_p12_eta3p2_To2HPip_Pythia8_pp_14TeV_GEN_20201102.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')

process.generator = cms.EDFilter("Pythia8EGun",
    PGunParameters = cms.PSet(
        AddAntiParticle = cms.bool(True),
        MaxE = cms.double(12.37),
        MaxEta = cms.double(3.2),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(3.0),
        MinEta = cms.double(-3.2),
        MinPhi = cms.double(-3.14159265359),
        ParticleID = cms.vint32(1010010030)
    ),
    PythiaParameters = cms.PSet(
        hyperNuclei = cms.vstring(
            'ParticleDecays:limitTau0 = off', 
            '1010010030:new = 3LambdaH 3LambdaHbar 2 3 0 2.99116 0 0 0 6.2e+01', 
            '1010010030:oneChannel = 1 1.0 0 1000010020 -211 2212'
        ),
        lightNuclei = cms.vstring(
            '1000010020:new = 2H 2Hbar 1 3 0 1.87561', 
            '1000010030:new = 3H 3Hbar 2 3 0 2.80925', 
            '1000020030:new = 3He 3Hebar 2 6 0 2.80923', 
            '1000020040:new = 4He 4Hebar 1 6 0 3.72742'
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CP5Settings', 
            'lightNuclei', 
            'hyperNuclei'
        ),
        pythia8CP5Settings = cms.vstring(
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:ecmPow=0.03344', 
            'PDF:pSet=20', 
            'MultipartonInteractions:bProfile=2', 
            'MultipartonInteractions:pT0Ref=1.41', 
            'MultipartonInteractions:coreRadius=0.7634', 
            'MultipartonInteractions:coreFraction=0.63', 
            'ColourReconnection:range=5.176', 
            'SigmaTotal:zeroAXB=off', 
            'SpaceShower:alphaSorder=2', 
            'SpaceShower:alphaSvalue=0.118', 
            'SigmaProcess:alphaSvalue=0.118', 
            'SigmaProcess:alphaSorder=2', 
            'MultipartonInteractions:alphaSvalue=0.118', 
            'MultipartonInteractions:alphaSorder=2', 
            'TimeShower:alphaSorder=2', 
            'TimeShower:alphaSvalue=0.118'
        ),
        pythia8CommonSettings = cms.vstring(
            'Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'
        )
    ),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

process.signal = process.generator.clone()
from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import pythia8CP5SettingsBlock
process.background = cms.EDFilter("Pythia8GeneratorFilter",
    crossSection = cms.untracked.double(71.39e+09),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(14000.0),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring("SoftQCD:inelastic = on"),
        parameterSets = cms.vstring("pythia8CommonSettings",
                                    "pythia8CP5Settings",
                                    "processParameters",
                                    )
        )
)
process.RandomNumberGeneratorService.signal = process.RandomNumberGeneratorService.generator.clone()
process.RandomNumberGeneratorService.background = process.RandomNumberGeneratorService.generator.clone()
delattr(process.RandomNumberGeneratorService, "generator")
process.generator = cms.EDProducer("HepMCMerger", signal = cms.string("signal"), background = cms.string("background"))
for path in process.paths:
    getattr(process,path).replace(process.generator, process.signal+process.background+process.generator)
