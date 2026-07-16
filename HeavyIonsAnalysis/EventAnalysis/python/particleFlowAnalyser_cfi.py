import FWCore.ParameterSet.Config as cms

particleFlowAnalyser = cms.EDAnalyzer(
    "ParticleFlowAnalyser",
    pfCandidateSrc = cms.InputTag('packedPFCandidates'),
    addInfo = cms.bool(False),
    ptMin = cms.double(5.),
    absEtaMax = cms.double(5.),
    )
