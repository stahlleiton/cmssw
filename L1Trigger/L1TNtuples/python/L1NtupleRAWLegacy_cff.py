import FWCore.ParameterSet.Config as cms

from L1Trigger.L1TNtuples.l1EventTree_cfi import *
from L1Trigger.L1TNtuples.l1ExtraTree_cfi import *
from L1Trigger.L1TNtuples.l1UpgradeTree_cfi import *

l1legacyMuonInStage2FormatTree = l1UpgradeTree.clone()
l1legacyMuonInStage2FormatTree.egToken   = cms.untracked.InputTag("caloStage1FinalDigis")
l1legacyMuonInStage2FormatTree.jetToken  = cms.untracked.InputTag("caloStage1FinalDigis")
l1legacyMuonInStage2FormatTree.sumToken  = cms.untracked.InputTag("caloStage1FinalDigis")
l1legacyMuonInStage2FormatTree.tauToken  = cms.untracked.InputTag("caloStage1FinalDigis","rlxTaus")
l1legacyMuonInStage2FormatTree.muonToken = cms.untracked.InputTag("rawMuonLegacyInStage2FormatDigis","imdMuonsLegacy")

L1NtupleRAWLegacy = cms.Sequence(
  l1EventTree
  +l1ExtraTree
  +l1legacyMuonInStage2FormatTree
)
