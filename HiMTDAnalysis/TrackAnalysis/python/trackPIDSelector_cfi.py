import FWCore.ParameterSet.Config as cms

trackPIDSelector = cms.EDProducer(
    'TrackPIDSelector',
    particle           = cms.string(""),
    maxInvBetaSignificance = cms.double(1),
    maxDeDxSignificanceInMTD = cms.double(0),
    maxDeDxSignificanceOutMTD = cms.double(2.5),
    trackSelection     = cms.string(""),
    recoTracksTag      = cms.InputTag("generalTracks"),
    primaryVertexTag   = cms.InputTag("offlinePrimaryVertices4D"),
    beamSpotTag        = cms.InputTag("offlineBeamSpot"),
    trackTMTDTag       = cms.InputTag("trackExtenderWithMTD:generalTracktmtd"),
    trackSigmaTMTDTag  = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd"),
    trackPathLengthTag = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength"),
    dEdxTag            = cms.InputTag("dedxPixelHarmonic2"),
)

pion = trackPIDSelector.clone(particle =  cms.string("Pion"))

kaon = trackPIDSelector.clone(particle =  cms.string("Kaon"))

proton = trackPIDSelector.clone(particle =  cms.string("Proton"))
