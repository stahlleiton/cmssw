import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.ZDC2018Gain_cfg import gainFill7435

QWzdcreco = cms.EDProducer('QWZDC2018RecHit',
        srcDetId = cms.untracked.InputTag("zdcdigi", "DetId"),
        srcfC = cms.untracked.InputTag("zdcdigi", "regularfC"),
        ZDCCalib = gainFill7435
        )

