import FWCore.ParameterSet.Config as cms

convAna = cms.EDAnalyzer(
    'HiConversionAnalyzer',
    recoMuonsTag       = cms.InputTag("muons"),
    recoConversionsTag = cms.InputTag("allConversions"),
    genParticlesTag    = cms.InputTag("genParticles"),
    pfCandidatesTag    = cms.InputTag("particleFlow"),
    primaryVertexTag   = cms.InputTag("offlinePrimaryVertices"),
    beamSpotTag        = cms.InputTag("offlineBeamSpot"),
    doAll              = cms.bool(False)
    )

convAnaSeq = cms.Sequence( convAna )
