import FWCore.ParameterSet.Config as cms

muonAna = cms.EDAnalyzer(
    'HiMuonAnalyzer',
    patMuonsTag        = cms.InputTag("patMuonsWithTriggers"),
    recoMuonsTag       = cms.InputTag("muons"),
    genParticlesTag    = cms.InputTag("genParticles"),
    pfCandidatesTag    = cms.InputTag("particleFlow"),
    primaryVertexTag   = cms.InputTag("offlinePrimaryVertices"),
    beamSpotTag        = cms.InputTag("offlineBeamSpot"),
    doAll              = cms.bool(False),
    triggerResultsTag  = cms.InputTag("TriggerResults","","HLT"),
    triggerPathNames   = cms.vstring(
        "HLT_PAL2Mu12_v1",
        "HLT_PAL2Mu15_v1",
        "HLT_PAL3Mu3_v1",
        "HLT_PAL3Mu5_v3",
        "HLT_PAL3Mu7_v1",
        "HLT_PAL3Mu12_v1",
        "HLT_PAL3Mu15_v1"
        ),
    triggerFilterNames = cms.vstring(
        "hltL2fL1sSingleMu7BptxANDL1f0L2Filtered12",
        "hltL2fL1sSingleMu7BptxANDL1f0L2Filtered15",
        "hltL3fL1sSingleMu3BptxANDL1f0L2f0L3Filtered3",
        "hltL3fL1sSingleMu5BptxANDL1f0L2f0L3Filtered5",
        "hltL3fL1sSingleMu5BptxANDL1f0L2f0L3Filtered7",
        "hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered12",
        "hltL3fL1sSingleMu7BptxANDL1f0L2f0L3Filtered15"
        )                                
    )

muonAnaSeq = cms.Sequence( muonAna )
