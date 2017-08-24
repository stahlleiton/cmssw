import FWCore.ParameterSet.Config as cms

################################################################################### 
# pp iterative tracking modified for hiOffline reco (the vertex is the one reconstructed in HI)
################################### 3rd step: low pT and detached tracks from pixel quadruplets

from RecoHI.HiTracking.HITrackingRegionProducer_cfi import *
hiRegitMuDetachedQuadStepTrackingRegions = HiTrackingRegionFactoryFromSTAMuonsEDProducer.clone(
    MuonTrackingRegionBuilder = dict(
        # FIXME: NEED TO OPTIMIZE THIS PARAMETERS. THESE WERE TAKEN FROM PIXEL TRIPLETS
        Pt_min        = 0.9,
        DeltaR        = 2.0, # default = 0.2
        DeltaZ        = 2.0, # this give you the length
        Rescale_Dz    = 4., # max(DeltaZ_Region,Rescale_Dz*vtx->zError())
    )
)

###################################
from RecoTracker.IterativeTracking.DetachedQuadStep_cff import *

# NEW CLUSTERS (remove previously used clusters)
from RecoLocalTracker.SubCollectionProducers.trackClusterRemover_cfi import trackClusterRemover as _trackClusterRemover
hiRegitMuDetachedQuadStepClusters = _trackClusterRemover.clone(
    maxChi2                                  = 9.0,
    pixelClusters                            = "siPixelClusters",
    stripClusters                            = "siStripClusters",
    trajectories          		     = cms.InputTag("hiRegitMuPixelLessStepTracks"),
    overrideTrkQuals      		     = cms.InputTag('hiRegitMuPixelLessStepSelector','hiRegitMuPixelLessStep'),
    TrackQuality                             = 'tight',
    trackClassifier       		     = cms.InputTag(''),
    minNumberOfLayersWithMeasBeforeFiltering = 0
)


# SEEDING LAYERS
hiRegitMuDetachedQuadStepSeedLayers =  RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepSeedLayers.clone()
hiRegitMuDetachedQuadStepSeedLayers.BPix.skipClusters = cms.InputTag('hiRegitMuDetachedQuadStepClusters')
hiRegitMuDetachedQuadStepSeedLayers.FPix.skipClusters = cms.InputTag('hiRegitMuDetachedQuadStepClusters')

# seeding
hiRegitMuDetachedQuadStepHitDoublets = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepHitDoublets.clone(
    seedingLayers = "hiRegitMuDetachedQuadStepSeedLayers",
    trackingRegions = "hiRegitMuDetachedQuadStepTrackingRegions",
    clusterCheck = "hiRegitMuClusterCheck",
)

hiRegitMuDetachedQuadStepHitQuadruplets = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepHitQuadruplets.clone(
    doublets = "hiRegitMuDetachedQuadStepHitDoublets",
)

from Configuration.Eras.Modifier_trackingPhase1QuadProp_cff import trackingPhase1QuadProp
hiRegitMuDetachedQuadStepHitTriplets = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepHitTriplets.clone(
    doublets = "hiRegitMuDetachedQuadStepHitDoublets",
)
trackingPhase1QuadProp.toModify(hiRegitMuDetachedQuadStepHitQuadruplets, triplets = "hiRegitMuDetachedQuadStepHitTriplets")

hiRegitMuDetachedQuadStepSeeds = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepSeeds.clone(
    seedingHitSets = "hiRegitMuDetachedQuadStepHitQuadruplets"
)


# QUALITY CUTS DURING TRACK BUILDING
hiRegitMuDetachedQuadStepTrajectoryFilterBase = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepTrajectoryFilterBase.clone()
hiRegitMuDetachedQuadStepTrajectoryFilterBase.minPt = 0.8 # after each new hit, apply pT cut for traj w/ at least minHitsMinPt = cms.int32(3),

hiRegitMuDetachedQuadStepTrajectoryFilter = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepTrajectoryFilter.clone()
hiRegitMuDetachedQuadStepTrajectoryFilter.filters = cms.VPSet(
      cms.PSet( refToPSet_ = cms.string('hiRegitMuDetachedQuadStepTrajectoryFilterBase') ),
      #cms.PSet( refToPSet_ = cms.string('detachedQuadStepTrajectoryFilterShape') )
      )

hiRegitMuDetachedQuadStepTrajectoryBuilder = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepTrajectoryBuilder.clone(
    trajectoryFilter = dict(refToPSet_ = 'hiRegitMuDetachedQuadStepTrajectoryFilter'),
    #clustersToSkip = cms.InputTag('hiRegitMuDetachedQuadStepClusters')
)

hiRegitMuDetachedQuadStepTrackCandidates = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepTrackCandidates.clone(
    src = 'hiRegitMuDetachedQuadStepSeeds',
    clustersToSkip = cms.InputTag("hiRegitMuDetachedQuadStepClusters")
    TrajectoryBuilderPSet = dict(refToPSet_ = 'hiRegitMuDetachedQuadStepTrajectoryBuilder'),
    )

# fitting: feed new-names
hiRegitMuDetachedQuadStepTracks = RecoTracker.IterativeTracking.DetachedQuadStep_cff.detachedQuadStepTracks.clone(
    AlgorithmName = cms.string('hiRegitMuDetachedQuadStep'),
    src = 'hiRegitMuDetachedQuadStepTrackCandidates'
)


import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
import RecoHI.HiTracking.hiMultiTrackSelector_cfi
hiRegitMuDetachedQuadStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
    src       = 'hiRegitMuDetachedQuadStepTracks',
    vertices  = cms.InputTag("hiSelectedPixelVertex"),
    useAnyMVA = cms.bool(True),
    GBRForestLabel = cms.string('HIMVASelectorIter7'),
    GBRForestVars = cms.vstring(['chi2perdofperlayer', 'nhits', 'nlayers', 'eta']),
    trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
           name = 'hiRegitMuDetachedQuadStepLoose',
           min_nhits = cms.uint32(8)
            ),
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
            name = 'hiRegitMuDetachedQuadStepTight',
            preFilterName = 'hiRegitMuDetachedQuadStepLoose',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(True),
            minMVA = cms.double(-0.2)
            ),
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
            name = 'hiRegitMuDetachedQuadStep',
            preFilterName = 'hiRegitMuDetachedQuadStepTight',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(True),
            minMVA = cms.double(-0.09)
            )
        ) #end of vpset
    )
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
trackingPhase1.toModify(hiRegitMuDetachedQuadStepSelector, useAnyMVA = cms.bool(False))
trackingPhase1.toModify(hiRegitMuDetachedQuadStepSelector, trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
           name = 'hiRegitMuDetachedQuadStepLoose',
           min_nhits = cms.uint32(8)
            ),
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
            name = 'hiRegitMuDetachedQuadStepTight',
            preFilterName = 'hiRegitMuDetachedQuadStepLoose',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(False),
            minMVA = cms.double(-0.2)
            ),
        RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
            name = 'hiRegitMuDetachedQuadStep',
            preFilterName = 'hiRegitMuDetachedQuadStepTight',
            min_nhits = cms.uint32(8),
            useMVA = cms.bool(False),
            minMVA = cms.double(-0.09)
            )
        )
)

hiRegitMuonDetachedQuadStep = cms.Sequence(hiRegitMuDetachedQuadStepClusters*
                                           hiRegitMuDetachedQuadStepSeedLayers*
                                           hiRegitMuDetachedQuadStepTrackingRegions*
                                           hiRegitMuDetachedQuadStepHitDoublets*
                                           hiRegitMuDetachedQuadStepHitQuadruplets*
                                           hiRegitMuDetachedQuadStepSeeds*
                                           hiRegitMuDetachedQuadStepTrackCandidates*
                                           hiRegitMuDetachedQuadStepTracks*
                                           hiRegitMuDetachedQuadStepSelector
                                           )

_hiRegitMuonDetachedQuadStep_Phase1Prop = hiRegitMuonDetachedQuadStep.copy()
_hiRegitMuonDetachedQuadStep_Phase1Prop.replace(hiRegitMuDetachedQuadStepHitDoublets, hiRegitMuDetachedQuadStepHitDoublets+hiRegitMuDetachedQuadStepHitTriplets)
trackingPhase1QuadProp.toReplaceWith(hiRegitMuonDetachedQuadStep, _hiRegitMuonDetachedQuadStep_Phase1Prop)
