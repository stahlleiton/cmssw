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
process.options.numberOfThreads=cms.untracked.uint32(1)

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/jiazhao/STARlight/2023Run3/Reco/STARlight_CohPhiToKK_Reco_132X_240125_044529/STARlight/CohPhiToKK/240125_034539/0000/step3_STARlight_Reco_10.root"),
    # fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/MC/STARLIGHT/2024_01_19/STARLIGHT_5p36TeV_2023Run3/double_diff_STARLIGHT_5p36TeV_2023Run3_UPCRECO_2024_01_19/240119_020712/0000/STARLIGHT_double_diff_RECO_10.root"),

)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_mcRun3_2023_realistic_HI_v9')


## ##############################################################################################################################
## Variables Production #########################################################################################################

#* Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

generalTrackParticles = generalParticles.clone(
    # recoToSimTrackMap = cms.InputTag('trackingParticleRecoTrackAsssociation')
)

process.tracks = generalTrackParticles.clone(
    tracks = cms.InputTag('generalTracks'),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
)

process.pixelTracks = generalParticles.clone(
    tracks = cms.InputTag('hiConformalPixelTracks'),
    # recoToSimTrackMap = cms.InputTag('trackingParticlePixelTrackAsssociation')
)

process.pfCandidates = generalTrackParticles.clone(
    pfParticles = cms.InputTag('particleFlow'),
    tracks = cms.InputTag('')
)


# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hiClusterCompatibility * process.primaryVertexFilter)

# Add sim-reco matching
# process.load('SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff')
# process.load('SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi')
# process.load('SimTracker.TrackerHitAssociation.tpClusterProducer_cfi')
# process.trackingParticlePixelTrackAsssociation = process.trackingParticleRecoTrackAsssociation.clone(label_tr = cms.InputTag("hiConformalPixelTracks"))
# process.simRecoTrackAssocSeq = cms.Sequence(process.tpClusterProducer * process.quickTrackAssociatorByHitsTrackerHitAssociator * process.trackingParticleRecoTrackAsssociation * process.trackingParticlePixelTrackAsssociation)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.colEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# process.track_step = cms.Path(process.simRecoTrackAssocSeq * process.tracks * process.pixelTracks * process.pfCandidates)
process.track_step = cms.Path(process.tracks * process.pixelTracks * process.pfCandidates)

## Adding the VertexComposite tree ################################################################################################

event_filter = cms.untracked.vstring(
        "Flag_colEvtSel",
        "Flag_clusterCompatibilityFilter",
        "Flag_primaryVertexFilter",
        "Flag_hfPosFilterNTh7",
        "Flag_hfPosFilterNTh7p3",
        "Flag_hfPosFilterNTh8",
        "Flag_hfPosFilterNTh10",
        "Flag_hfNegFilterNTh7",
        "Flag_hfNegFilterNTh7p6",
        "Flag_hfNegFilterNTh8",
        "Flag_hfNegFilterNTh10",
    )

trig_info = cms.untracked.VPSet([
    # zero bias triggers
    cms.PSet(path = cms.string('HLT_HIZeroBias_v*')),
    cms.PSet(path = cms.string('HLT_HIZeroBias_HighRate_v*')),
    # UPC low pT triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # UPC ZDC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*')),
  ])

from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc

process.trackAna = particleAna_mc.clone(
    recoParticles = cms.InputTag("tracks"),
    selectEvents = cms.string(""),
    maxGenDeltaR = cms.untracked.double(0.3),
    maxGenDeltaPtRel = cms.untracked.double(0.5),
    eventFilterNames = event_filter,
    triggerInfo = trig_info,
)

process.pixelTrackAna = particleAna_mc.clone(
	recoParticles = cms.InputTag("pixelTracks"),
    selectEvents = cms.string(""),
    maxGenDeltaR = cms.untracked.double(0.3),
    maxGenDeltaPtRel = cms.untracked.double(0.5),
    eventFilterNames = event_filter,
    triggerInfo = trig_info,
)

process.pfCandidatesAna = particleAna_mc.clone(
	recoParticles = cms.InputTag("pfCandidates"),
	selectEvents = cms.string(""),
	maxGenDeltaR = cms.untracked.double(0.3),
	maxGenDeltaPtRel = cms.untracked.double(0.5),
	eventFilterNames = event_filter,
	triggerInfo = trig_info,
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('track_ana_mc.root'))
process.p = cms.EndPath(process.trackAna * process.pixelTrackAna * process.pfCandidatesAna)

#! Define the process schedule !!!!!!!!!!!!!!!!!!
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.track_step,
    process.p
)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Add the event selection filters ###############################################################################################
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_clusterCompatibilityFilter = cms.Path(process.hiClusterCompatibility)
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.Flag_hfPosFilterNTh7 = cms.Path(process.hfPosFilterNTh7_seq)
process.Flag_hfPosFilterNTh7p3 = cms.Path(process.hfPosFilterNTh7p3_seq)
process.Flag_hfPosFilterNTh8 = cms.Path(process.hfPosFilterNTh8_seq)
process.Flag_hfPosFilterNTh10 = cms.Path(process.hfPosFilterNTh10_seq)
process.Flag_hfNegFilterNTh7 = cms.Path(process.hfNegFilterNTh7_seq)
process.Flag_hfNegFilterNTh7p6 = cms.Path(process.hfNegFilterNTh7p6_seq)
process.Flag_hfNegFilterNTh8 = cms.Path(process.hfNegFilterNTh8_seq)
process.Flag_hfNegFilterNTh10 = cms.Path(process.hfNegFilterNTh10_seq)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_clusterCompatibilityFilter , process.Flag_primaryVertexFilter , process.Flag_hfPosFilterNTh7 , process.Flag_hfPosFilterNTh7p3 , process.Flag_hfPosFilterNTh8 , process.Flag_hfPosFilterNTh10 , process.Flag_hfNegFilterNTh7 , process.Flag_hfNegFilterNTh7p6 , process.Flag_hfNegFilterNTh8 , process.Flag_hfNegFilterNTh10 ]

#! Adding the process schedule !!!!!!!!!!!!!!!!!!
for P in eventFilterPaths:
    process.schedule.insert(0, P)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!