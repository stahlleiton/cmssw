import FWCore.ParameterSet.Config as cms

timeAna = cms.EDAnalyzer(
    'TimeAnalyzer',
    recoTracksTag      = cms.InputTag("generalTracks"),
    primaryVertexTag   = cms.InputTag("offlinePrimaryVertices"),
    beamSpotTag        = cms.InputTag("offlineBeamSpot"),
    trackBetaTag       = cms.InputTag("trackExtenderWithMTD:generalTrackBeta"),
    trackT0Tag         = cms.InputTag("trackExtenderWithMTD:generalTrackt0"),
    trackSigmaT0Tag    = cms.InputTag("trackExtenderWithMTD:generalTracksigmat0"),
    trackTMTDTag       = cms.InputTag("trackExtenderWithMTD:generalTracktmtd"),
    trackSigmaTMTDTag  = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd"),
    trackMomTag        = cms.InputTag("trackExtenderWithMTD:generalTrackp"),
    trackPathLengthTag = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength")
)

timeAnaSeq = cms.Sequence( timeAna )


