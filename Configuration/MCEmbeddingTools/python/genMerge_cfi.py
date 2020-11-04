import FWCore.ParameterSet.Config as cms

genMerge = cms.EDProducer("HepMCMerger",
    signal = cms.string("signal"),
    background = cms.string("background"),
)
