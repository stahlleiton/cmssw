import FWCore.ParameterSet.Config as cms

from L1Trigger.L1TNtuples.l1UpgradeTree_cfi import *
from L1Trigger.L1TNtuples.l1EventTree_cfi import *

l1legacyMuonEmuTree = l1UpgradeTree.clone()
l1legacyMuonEmuTree.egToken   = cms.untracked.InputTag("caloStage1FinalDigis")
l1legacyMuonEmuTree.jetToken  = cms.untracked.InputTag("caloStage1FinalDigis")
l1legacyMuonEmuTree.sumToken  = cms.untracked.InputTag("caloStage1FinalDigis")
l1legacyMuonEmuTree.tauToken  = cms.untracked.InputTag("caloStage1FinalDigis","rlxTaus")
l1legacyMuonEmuTree.muonToken = cms.untracked.InputTag("simMuonLegacyInStage2FormatDigis","imdMuonsLegacy")

L1NtupleEMULegacy = cms.Sequence(
  l1EventTree
  +l1legacyMuonEmuTree
)

