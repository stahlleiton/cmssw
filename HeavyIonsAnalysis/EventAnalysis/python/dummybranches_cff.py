import FWCore.ParameterSet.Config as cms

dummy_branches_for_PbPb_2023_HLT = cms.vstring([
    'HLT_HIMinimumBiasHF1AND_v3',
    'HLT_HIMinimumBiasHF1ANDZDC1nOR_v1',
    'HLT_HIMinimumBiasHF1ANDZDC2nOR_v3',
    'HLT_HIEle20Gsf_v10',
    'HLT_HIL2SingleMu7_v3',
])

dummy_branches_for_PbPb_2024_HLT = cms.vstring([
    'HLT_HIMinimumBiasHF1AND_v7',
    'HLT_HIMinimumBiasHF1ANDZDC1nOR_v4',
    'HLT_HIGEDPhoton10_v14',
    'HLT_HIEle20Gsf_v14',
    'HLT_HIL2SingleMu7_v7',
])

dummy_branches_for_ppRef_2024_HLT = cms.vstring([
    'HLT_PPRefGEDPhoton30_v6',
    'HLT_PPRefEle20Gsf_v6',
    'HLT_PPRefL2SingleMu7_v6',
])

dummy_branches_for_pO_2025_HLT = cms.vstring([
    'HLT_MinimumBiasHF_OR_BptxAND_v1',
    'HLT_OxyL1SingleEG10_v1',
    'HLT_OxyL1SingleEG15_v1',
    'HLT_OxyL1SingleMuOpen_v1',
    'HLT_OxyL1SingleMu0_v1',
])
dummy_branches_for_OO_2025_HLT = dummy_branches_for_pO_2025_HLT.copy()
dummy_branches_for_NeNe_2025_HLT = dummy_branches_for_pO_2025_HLT.copy()

dummy_branches_for_PbPb_2025_HLT = cms.vstring([
    'HLT_HIMinimumBiasHF1AND_v8',
    'HLT_HIMinimumBiasHF1ANDZDC1nOR_v6',
    'HLT_HIGEDPhoton10_v16',
    'HLT_HIL2SingleMu7_v8',
])

dummy_branches_for_PbPb_2026_HLT = cms.vstring([
    'HLT_HIMinimumBiasHF1AND_v8',
    'HLT_HIMinimumBiasHF1ANDZDC1nOR_v6',
    'HLT_HIGEDPhoton10_v17',
    'HLT_HIL2SingleMu7_v9',
])
