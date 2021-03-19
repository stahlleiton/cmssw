import FWCore.ParameterSet.Config as cms

timeAna = cms.EDAnalyzer(
    'TimeAnalyzer',
    recoTracksTag      = cms.InputTag("generalTracks"),
    primaryVertexTag   = cms.InputTag("offlinePrimaryVertices4D::RECO"),
    beamSpotTag        = cms.InputTag("offlineBeamSpot"),
    trackBetaTag       = cms.InputTag("trackExtenderWithMTD:generalTrackBeta:RECO"),
    trackT0Tag         = cms.InputTag("trackExtenderWithMTD:generalTrackt0:RECO"),
    trackSigmaT0Tag    = cms.InputTag("trackExtenderWithMTD:generalTracksigmat0:RECO"),
    trackTMTDTag       = cms.InputTag("trackExtenderWithMTD:generalTracktmtd:RECO"),
    trackSigmaTMTDTag  = cms.InputTag("trackExtenderWithMTD:generalTracksigmatmtd:RECO"),
    trackMomTag        = cms.InputTag("trackExtenderWithMTD:generalTrackp:RECO"),
    trackPathLengthTag = cms.InputTag("trackExtenderWithMTD:generalTrackPathLength:RECO"),
    tofPIDT0Tag        = cms.InputTag("tofPID:t0:RECO"),
    tofPIDSigmaT0Tag   = cms.InputTag("tofPID:sigmat0:RECO"),
    dEdxTags           = cms.VInputTag("dedxPixelHarmonic2"),
    genParticlesTag    = cms.InputTag("genParticles"),
    genT0Tag           = cms.InputTag("genParticles:t0"),
    centralitySrcTag   = cms.InputTag("hiCentrality"),
    centralityBinTag   = cms.InputTag("centralityBin:HFtowers"),
    associatorMap      = cms.InputTag("trackingParticleRecoTrackAsssociation"),
)

timeAnaSeq = cms.Sequence( timeAna )


