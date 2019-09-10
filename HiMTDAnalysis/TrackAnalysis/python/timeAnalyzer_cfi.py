import FWCore.ParameterSet.Config as cms

timeAna = cms.EDAnalyzer(
    'TimeAnalyzer',
    recoTracksTag      = cms.InputTag("generalTracks"),
    primaryVertexTag   = cms.InputTag("offlinePrimaryVertices4D"),
    beamSpotTag        = cms.InputTag("offlineBeamSpot"),
    trackBetaTag       = cms.InputTag("trackExtenderWithMTD:generalTrackBeta"),
    trackT0Tag         = cms.InputTag("trackExtenderWithMTD:generalTrackt0"),
    trackSigmaT0Tag    = cms.InputTag("trackExtenderWithMTD:generalTracksigmat0"),
    trackTMTDTag       = cms.InputTag("trackExtenderWithMTD:generalTracktmtd"),
    trackSigmaTMTDTag  = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd"),
    trackMomTag        = cms.InputTag("trackExtenderWithMTD:generalTrackp"),
    trackPathLengthTag = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength"),
    tofPIDT0Tag        = cms.InputTag("tofPID:t0"),
    tofPIDSigmaT0Tag   = cms.InputTag("tofPID:sigmat0"),
    dEdxTag            = cms.InputTag("dedxPixelHarmonic2"),
    genParticlesTag    = cms.InputTag("genParticles"),
    genT0Tag           = cms.InputTag("genParticles:t0"),
    centralitySrcTag   = cms.InputTag("hiCentrality"),
    centralityBinTag   = cms.InputTag("centralityBin:HFtowers"),
    doDiTracks         = cms.bool(False)
)

timeAnaSeq = cms.Sequence( timeAna )


