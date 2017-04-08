import FWCore.ParameterSet.Config as cms

from RecoHI.HiTracking.hiMixedTripletStep_cff import *
from RecoTracker.IterativeTracking.PixelLessStep_cff import *
from RecoTracker.IterativeTracking.TobTecStep_cff import *
from RecoTracker.FinalTrackSelectors.MergeTrackCollections_cff import *
from RecoTracker.ConversionSeedGenerators.ConversionStep_cff import *
from RecoEgamma.EgammaPhotonProducers.conversionTracks_cff import *
from RecoEgamma.EgammaPhotonProducers.conversionTrackSequence_cff import *
from RecoEgamma.EgammaPhotonProducers.allConversions_cfi import *

pixelLessStepClusters.oldClusterRemovalInfo = cms.InputTag("hiMixedTripletClusters")
pixelLessStepClusters.trackClassifier       = cms.InputTag("")
pixelLessStepClusters.overrideTrkQuals      = cms.InputTag("hiMixedTripletStepSelector","hiMixedTripletStep")
pixelLessStepClusters.trajectories          = cms.InputTag("hiMixedTripletGlobalPrimTracks")

photonConvTrajSeedFromSingleLeg.TrackRefitter = cms.InputTag("hiGeneralTracks")
generalConversionTrackProducer.TrackProducer  = cms.string('hiGeneralTracks')

photonConvTrajSeedFromSingleLeg.primaryVerticesTag = cms.InputTag("hiSelectedVertex")
convStepSelector.vertices                          = cms.InputTag("hiFirstStepGoodPrimaryVertices")
pixelLessStepSelector.vertices                     = cms.InputTag("hiFirstStepGoodPrimaryVertices")
tobTecStepSelector.vertices                        = cms.InputTag("hiFirstStepGoodPrimaryVertices")
conv2StepSelector.vertices                         = cms.InputTag("hiFirstStepGoodPrimaryVertices")

allConversions.primaryVertexProducer = cms.InputTag("hiSelectedVertex")


hiConversionSequence = cms.Sequence(hiMixedTripletStep
                                    * PixelLessStep
                                    * TobTecStep 
                                    * ConvStep 
                                    * conversionStepTracks 
                                    * conversionTrackSequence 
                                    * allConversions
                                    )
