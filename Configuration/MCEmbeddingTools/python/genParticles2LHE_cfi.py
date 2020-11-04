import FWCore.ParameterSet.Config as cms

genParticles2LHE = cms.EDProducer("GenParticles2LHEConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator"),
    fileName = cms.untracked.string(""),
    cmEnergy = cms.untracked.double(13000),
    beamID = cms.untracked.int32(2212)
)
