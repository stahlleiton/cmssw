import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedGenerator.trackerClusterCheck_cfi import trackerClusterCheck as _trackerClusterCheck
hiRegitMuClusterCheck = _trackerClusterCheck.clone(
    doClusterCheck = False # do not check for max number of clusters pixel or strips
)

from RecoHI.HiMuonAlgos.HiRegitMuonInitialStep_cff import *
from RecoHI.HiMuonAlgos.HiRegitMuonPixelPairStep_cff import *
from RecoHI.HiMuonAlgos.HiRegitMuonDetachedQuadStep_cff import *
from RecoHI.HiMuonAlgos.HiRegitMuonDetachedTripletStep_cff import *
from RecoHI.HiMuonAlgos.HiRegitMuonMixedTripletStep_cff import *
from RecoHI.HiMuonAlgos.HiRegitMuonPixelLessStep_cff import *
from RecoHI.HiMuonAlgos.HiRegitMuonSeededStep_cff import *

from RecoTracker.FinalTrackSelectors.trackAlgoPriorityOrder_cfi import trackAlgoPriorityOrder
import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
hiGeneralAndRegitMuTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = (cms.InputTag('hiRegitMuInitialStepTracks'),
                      cms.InputTag('hiRegitMuPixelPairStepTracks'),
                      cms.InputTag('hiRegitMuMixedTripletStepTracks'),
                      cms.InputTag('hiRegitMuPixelLessStepTracks'),
                      cms.InputTag('hiRegitMuDetachedTripletStepTracks'),
                      cms.InputTag('hiRegitMuonSeededTracksOutIn'),
                      cms.InputTag('hiRegitMuonSeededTracksInOut')
                      ),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("hiRegitMuInitialStepSelector","hiRegitMuInitialStepLoose"),
                                       cms.InputTag("hiRegitMuPixelPairStepSelector","hiRegitMuPixelPairStep"),
                                       cms.InputTag("hiRegitMuMixedTripletStepSelector","hiRegitMuMixedTripletStep"),
                                       cms.InputTag("hiRegitMuPixelLessStepSelector","hiRegitMuPixelLessStep"),
                                       cms.InputTag("hiRegitMuDetachedTripletStepSelector","hiRegitMuDetachedTripletStep"),
                                       cms.InputTag("hiRegitMuonSeededTracksOutInSelector","hiRegitMuonSeededTracksOutInHighPurity"),
                                       cms.InputTag("hiRegitMuonSeededTracksInOutSelector","hiRegitMuonSeededTracksInOutHighPurity")
                                       ),
    hasSelector=cms.vint32(1,1,1,1,1,1,1),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5,6), pQual=cms.bool(True))),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
from Configuration.Eras.Modifier_trackingPhase1QuadProp_cff import trackingPhase1QuadProp
_forPhase1 = dict(
    TrackProducers = (cms.InputTag('hiRegitMuInitialStepTracks'),
                      cms.InputTag('hiRegitMuPixelPairStepTracks'),
                      cms.InputTag('hiRegitMuMixedTripletStepTracks'),
                      cms.InputTag('hiRegitMuPixelLessStepTracks'),
                      cms.InputTag('hiRegitMuDetachedQuadStepTracks'),
                      cms.InputTag('hiRegitMuDetachedTripletStepTracks'),
                      cms.InputTag('hiRegitMuonSeededTracksOutIn'),
                      cms.InputTag('hiRegitMuonSeededTracksInOut')
                      ),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("hiRegitMuInitialStepSelector","hiRegitMuInitialStepLoose"),
                                       cms.InputTag("hiRegitMuPixelPairStepSelector","hiRegitMuPixelPairStep"),
                                       cms.InputTag("hiRegitMuMixedTripletStepSelector","hiRegitMuMixedTripletStep"),
                                       cms.InputTag("hiRegitMuPixelLessStepSelector","hiRegitMuPixelLessStep"),
                                       cms.InputTag("hiRegitMuDetachedQuadStepSelector","hiRegitMuDetachedQuadStep"),
                                       cms.InputTag("hiRegitMuDetachedTripletStepSelector","hiRegitMuDetachedTripletStep"),
                                       cms.InputTag("hiRegitMuonSeededTracksOutInSelector","hiRegitMuonSeededTracksOutInHighPurity"),
                                       cms.InputTag("hiRegitMuonSeededTracksInOutSelector","hiRegitMuonSeededTracksInOutHighPurity")
                                       ),
    hasSelector=cms.vint32(1,1,1,1,1,1,1,1),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5,6,7), pQual=cms.bool(True)))
)
trackingPhase1.toModify(hiGeneralAndRegitMuTracks, **_forPhase1)
trackingPhase1QuadProp.toModify(hiGeneralAndRegitMuTracks, **_forPhase1)

hiRegitMuTracking = cms.Sequence(hiRegitMuClusterCheck
                                 *hiRegitMuonInitialStep
                                 *hiRegitMuonPixelPairStep
                                 *hiRegitMuonMixedTripletStep
                                 *hiRegitMuonPixelLessStep
                                 *hiRegitMuonDetachedTripletStep
                                 *hiRegitMuonSeededStep
                                 )

_hiRegitMuTracking_Phase1 = hiRegitMuTracking.copy()
_hiRegitMuTracking_Phase1.replace(hiRegitMuonDetachedTripletStep, hiRegitMuonDetachedTripletStep+hiRegitMuonDetachedQuadStep)
trackingPhase1.toReplaceWith(hiRegitMuTracking, _hiRegitMuTracking_Phase1)
trackingPhase1QuadProp.toReplaceWith(hiRegitMuTracking, _hiRegitMuTracking_Phase1)

# Standalone muons
from RecoMuon.Configuration.RecoMuonPPonly_cff import *


hiRegitMuTrackingAndSta = cms.Sequence(standalonemuontracking
      *hiRegitMuTracking)

