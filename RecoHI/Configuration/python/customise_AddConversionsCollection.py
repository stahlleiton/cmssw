import FWCore.ParameterSet.Config as cms

def customiseAddConversionsCollection(process):

    process.load("RecoHI.HiTracking.hiMixedTripletStep_cff")
    process.load("RecoTracker.IterativeTracking.PixelLessStep_cff")
    process.load("RecoTracker.IterativeTracking.TobTecStep_cff")
    process.load("RecoTracker.FinalTrackSelectors.MergeTrackCollections_cff")
    process.load("RecoTracker.ConversionSeedGenerators.ConversionStep_cff")
    process.load("RecoEgamma.EgammaPhotonProducers.conversionTracks_cff")
    process.load("RecoEgamma.EgammaPhotonProducers.conversionTrackSequence_cff")
    process.load("RecoEgamma.EgammaPhotonProducers.allConversions_cfi")
    
    process.pixelLessStepClusters.oldClusterRemovalInfo = cms.InputTag("hiMixedTripletClusters")
    process.pixelLessStepClusters.trackClassifier       = cms.InputTag("")
    process.pixelLessStepClusters.overrideTrkQuals      = cms.InputTag("hiMixedTripletStepSelector","hiMixedTripletStep")
    process.pixelLessStepClusters.trajectories          = cms.InputTag("hiMixedTripletGlobalPrimTracks")

    process.photonConvTrajSeedFromSingleLeg.TrackRefitter = cms.InputTag("hiGeneralTracks")
    process.generalConversionTrackProducer.TrackProducer  = cms.string('hiGeneralTracks')

    process.photonConvTrajSeedFromSingleLeg.primaryVerticesTag = cms.InputTag("hiSelectedVertex")
    process.convStepSelector.vertices                          = cms.InputTag("hiFirstStepGoodPrimaryVertices")
    process.pixelLessStepSelector.vertices                     = cms.InputTag("hiFirstStepGoodPrimaryVertices")
    process.tobTecStepSelector.vertices                        = cms.InputTag("hiFirstStepGoodPrimaryVertices")
    process.conv2StepSelector.vertices                         = cms.InputTag("hiFirstStepGoodPrimaryVertices")
    
    process.allConversions.primaryVertexProducer = cms.InputTag("hiSelectedVertex")

    process.hiConversionSequence = cms.Sequence(process.hiMixedTripletStep
                                                * process.PixelLessStep
                                                * process.TobTecStep 
                                                * process.ConvStep 
                                                * process.conversionStepTracks 
                                                * process.conversionTrackSequence 
                                                * process.allConversions
                                                )

    ###Add conversions in the sequence
    process.reconstruction_step += process.hiConversionSequence

    return process
