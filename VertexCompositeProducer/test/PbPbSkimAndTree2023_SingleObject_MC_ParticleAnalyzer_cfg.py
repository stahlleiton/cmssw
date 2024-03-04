import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023_UPC)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/a/anstahll/work/Run3_2023/Simulation/Production/Reconstruction/CMSSW_13_2_10/src/UPCReco/STARLIGHT/STARLIGHT_double_diff_RECO.root"),
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/STARLIGHT/2024_01_06/STARLIGHT_5p36TeV_2023Run3/double_diff_STARLIGHT_5p36TeV_2023Run3_UPCRECO_2024_01_06/240106_193637/0000/STARLIGHT_double_diff_RECO_1.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_mcRun3_2023_realistic_HI_v7')

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

generalTrackParticles = generalParticles.clone(
    recoToSimTrackMap = cms.InputTag('trackingParticleRecoTrackAsssociation')
)

process.muons = generalTrackParticles.clone(
    pdgId = cms.uint32(13),
    muons = cms.InputTag('patMuons'),
)

process.electrons = generalTrackParticles.clone(
    pdgId = cms.uint32(11),
    electrons = cms.InputTag('patElectrons')
)

process.lowPtElectrons = generalTrackParticles.clone(
    pdgId = cms.uint32(11),
    electrons = cms.InputTag('patLowPtElectrons')
)

process.photons = generalParticles.clone(
    pdgId = cms.uint32(22),
    photons = cms.InputTag('patPhotons')
)

process.convertedPhotons = generalParticles.clone(
    pdgId = cms.uint32(22),
    conversions = cms.InputTag('allConversions')
)

process.tracks = generalTrackParticles.clone(
    tracks = cms.InputTag('generalTracks'),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
)

process.pixelTracks = generalParticles.clone(
    tracks = cms.InputTag('hiConformalPixelTracks'),
    recoToSimTrackMap = cms.InputTag('trackingParticlePixelTrackAsssociation')
)

process.pfCandidates = generalTrackParticles.clone(
    pfParticles = cms.InputTag('particleFlow'),
    tracks = cms.InputTag('')
)

# Add PAT objects
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons, doPATElectrons, doPATPhotons
doPATMuons(process)
doPATElectrons(process)
doPATPhotons(process)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.colEvtSel = cms.Sequence(process.hiClusterCompatibility * process.primaryVertexFilter)

# Add sim-reco matching
process.load('SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff')
process.load('SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi')
process.load('SimTracker.TrackerHitAssociation.tpClusterProducer_cfi')
process.trackingParticlePixelTrackAsssociation = process.trackingParticleRecoTrackAsssociation.clone(label_tr = cms.InputTag("hiConformalPixelTracks"))
process.simRecoTrackAssocSeq = cms.Sequence(process.tpClusterProducer * process.quickTrackAssociatorByHitsTrackerHitAssociator * process.trackingParticleRecoTrackAsssociation * process.trackingParticlePixelTrackAsssociation)

# Define the analysis steps
process.muon_step = cms.Path(process.simRecoTrackAssocSeq * process.patMuonSequence * process.muons)
process.electron_step = cms.Path(process.simRecoTrackAssocSeq * process.patElectronSequence * process.electrons * process.lowPtElectrons)
process.photon_step = cms.Path(process.patPhotonSequence * process.photons * process.convertedPhotons)
process.track_step = cms.Path(process.simRecoTrackAssocSeq * process.tracks * process.pixelTracks * process.pfCandidates)

# Add the Particle tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc


process.muonAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("muons"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

process.elecAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("electrons"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.lowPtElecAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("lowPtElectrons"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.phoAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("photons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.convAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("convertedPhotons"),
  maxGenDeltaR = cms.untracked.double(0.3),
  maxGenDeltaPtRel = cms.untracked.double(1.0),
)

process.trackAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("tracks"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

process.pixelTrackAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("pixelTracks"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

process.pfAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("pfCandidates"),
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('obj_ana_mc.root'))
process.p = cms.EndPath(process.muonAna * process.elecAna * process.lowPtElecAna * process.phoAna * process.convAna * process.trackAna * process.pixelTrackAna * process.pfAna)

# Add fix for SuperChic
#process.genParticles = cms.EDProducer('GenParticleSuperChicFixer', genPars = cms.InputTag('genParticles::SIM'))
#process.gen_step = cms.Path(process.genParticles)

# Define the process schedule
process.schedule = cms.Schedule(
    #process.gen_step,
    process.muon_step,
    process.electron_step,
    process.photon_step,
    process.track_step,
    process.p
)
