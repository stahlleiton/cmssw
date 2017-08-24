import FWCore.ParameterSet.Config as cms

# pp iterative tracking modified for hiOffline reco (the vertex is the one reconstructed in HI)
################################### 0st step:pixel-triplet seeding, high-pT;
from RecoHI.HiTracking.HITrackingRegionProducer_cfi import *
hiRegitMuInitialStepTrackingRegions = HiTrackingRegionFactoryFromSTAMuonsEDProducer.clone(
    MuonTrackingRegionBuilder = dict(
        Pt_min          = 3.0,
        DeltaR          = 1, # default = 0.2
        DeltaZ          = 1, # this give you the length
        Rescale_Dz      = 4., # max(DeltaZ_Region,Rescale_Dz*vtx->zError())
    )
)

###################################  
from RecoTracker.IterativeTracking.InitialStep_cff import *

# SEEDING LAYERS
hiRegitMuInitialStepSeedLayers =  RecoTracker.IterativeTracking.InitialStep_cff.initialStepSeedLayers.clone()

# seeding
hiRegitMuInitialStepHitDoublets = RecoTracker.IterativeTracking.InitialStep_cff.initialStepHitDoublets.clone(
    seedingLayers = "hiRegitMuInitialStepSeedLayers",
    trackingRegions = "hiRegitMuInitialStepTrackingRegions",
    clusterCheck = "hiRegitMuClusterCheck"
)

hiRegitMuInitialStepHitTriplets = RecoTracker.IterativeTracking.InitialStep_cff.initialStepHitTriplets.clone(
    doublets = "hiRegitMuInitialStepHitDoublets"
)

hiRegitMuInitialStepHitQuadruplets = RecoTracker.IterativeTracking.InitialStep_cff.initialStepHitQuadruplets.clone(
    doublets = "hiRegitMuInitialStepHitDoublets"
)
from Configuration.Eras.Modifier_trackingPhase1QuadProp_cff import trackingPhase1QuadProp
trackingPhase1QuadProp.toModify(hiRegitMuInitialStepHitQuadruplets, triplets = "hiRegitMuInitialStepHitTriplets")

hiRegitMuInitialStepSeeds = RecoTracker.IterativeTracking.InitialStep_cff.initialStepSeeds.clone(
    seedingHitSets = "hiRegitMuInitialStepHitTriplets"
)
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
trackingPhase1.toModify(hiRegitMuInitialStepSeeds, seedingHitSets = "hiRegitMuInitialStepHitQuadruplets")


# building: feed the new-named seeds
hiRegitMuInitialStepTrajectoryFilterBase = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryFilterBase.clone()
hiRegitMuInitialStepTrajectoryFilterBase.minPt = 2.5 # after each new hit, apply pT cut for traj w/ at least minHitsMinPt = cms.int32(3),

hiRegitMuInitialStepTrajectoryFilter = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryFilter.clone()
hiRegitMuInitialStepTrajectoryFilter.filters = cms.VPSet(
      cms.PSet( refToPSet_ = cms.string('hiRegitMuInitialStepTrajectoryFilterBase')),
      #cms.PSet( refToPSet_ = cms.string('initialStepTrajectoryFilterShape'))
)

hiRegitMuInitialStepTrajectoryBuilder = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryBuilder.clone(
    trajectoryFilter = cms.PSet(
       refToPSet_ = cms.string('hiRegitMuInitialStepTrajectoryFilter')
       ),
)

hiRegitMuInitialStepTrajectoryFilterInOut = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrajectoryFilterInOut.clone()
trackingPhase1QuadProp.toModify(hiRegitMuInitialStepTrajectoryBuilder, inOutTrajectoryFilter = dict(refToPSet_ = "hiRegitMuInitialStepTrajectoryFilterInOut"))

# track candidates
hiRegitMuInitialStepTrackCandidates =  RecoTracker.IterativeTracking.InitialStep_cff.initialStepTrackCandidates.clone(
    src = cms.InputTag('hiRegitMuInitialStepSeeds'),
    TrajectoryBuilderPSet = cms.PSet(
       refToPSet_ = cms.string('hiRegitMuInitialStepTrajectoryBuilder')
       ),
    maxNSeeds = cms.uint32(1000000)
    )

# fitting: feed new-names
hiRegitMuInitialStepTracks = RecoTracker.IterativeTracking.InitialStep_cff.initialStepTracks.clone(
    src = 'hiRegitMuInitialStepTrackCandidates',
    AlgorithmName = cms.string('hiRegitMuInitialStep')
)


import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
import RecoHI.HiTracking.hiMultiTrackSelector_cfi
hiRegitMuInitialStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
    src                 ='hiRegitMuInitialStepTracks',
    vertices            = cms.InputTag("hiSelectedPixelVertex"),
    useAnyMVA = cms.bool(True),
    GBRForestLabel = cms.string('HIMVASelectorIter4'),
    GBRForestVars = cms.vstring(['chi2perdofperlayer', 'dxyperdxyerror', 'dzperdzerror', 'nhits', 'nlayers', 'eta']),
    trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
           name = 'hiRegitMuInitialStepLoose',
           min_nhits = cms.uint32(8)
            ), #end of pset
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
            name = 'hiRegitMuInitialStepTight',
            preFilterName = 'hiRegitMuInitialStepLoose',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(True),
            minMVA = cms.double(-0.38)
            ),
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
            name = 'hiRegitMuInitialStep',
            preFilterName = 'hiRegitMuInitialStepTight',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(True),
            minMVA = cms.double(-0.77)
            ),
        ) #end of vpset
    )
trackingPhase1.toModify(hiRegitMuInitialStepSelector, useAnyMVA = cms.bool(False))
trackingPhase1.toModify(hiRegitMuInitialStepSelector, trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
           name = 'hiRegitMuInitialStepLoose',
           min_nhits = cms.uint32(8)
            ), #end of pset
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
            name = 'hiRegitMuInitialStepTight',
            preFilterName = 'hiRegitMuInitialStepLoose',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(False),
            minMVA = cms.double(-0.38)
            ),
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
            name = 'hiRegitMuInitialStep',
            preFilterName = 'hiRegitMuInitialStepTight',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(False),
            minMVA = cms.double(-0.77)
            ),
        )
)

hiRegitMuonInitialStep = cms.Sequence(hiRegitMuInitialStepSeedLayers*
                                      hiRegitMuInitialStepTrackingRegions*
                                      hiRegitMuInitialStepHitDoublets*
                                      hiRegitMuInitialStepHitTriplets*
                                      hiRegitMuInitialStepSeeds*
                                      hiRegitMuInitialStepTrackCandidates*
                                      hiRegitMuInitialStepTracks*
                                      hiRegitMuInitialStepSelector)

_hiRegitMuonInitialStep_Phase1QuadProp = hiRegitMuonInitialStep.copy()
_hiRegitMuonInitialStep_Phase1QuadProp.replace(hiRegitMuInitialStepHitTriplets, hiRegitMuInitialStepHitTriplets+hiRegitMuInitialStepHitQuadruplets)
trackingPhase1QuadProp.toReplaceWith(hiRegitMuonInitialStep, _hiRegitMuonInitialStep_Phase1QuadProp)
_hiRegitMuonInitialStep_Phase1 = _hiRegitMuonInitialStep_Phase1QuadProp.copyAndExclude([hiRegitMuInitialStepHitTriplets])
trackingPhase1.toReplaceWith(hiRegitMuonInitialStep, _hiRegitMuonInitialStep_Phase1)
