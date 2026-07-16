import FWCore.ParameterSet.Config as cms

hltobject = cms.EDAnalyzer("TriggerObjectAnalyzer",
   processName = cms.string("HLT"),
   treeName = cms.string("JetTriggers"),
   triggerResults = cms.InputTag("TriggerResults","","HLT"),
   #triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
   triggerObjects = cms.InputTag("slimmedPatTrigger","")
)

trigger_list_data_2023_skimmed = cms.vstring(
   'HLT_HIEle20Gsf_v',
   'HLT_HIL2SingleMu7_v',
)

trigger_list_data_2024_skimmed = cms.vstring(
    'HLT_HIGEDPhoton10_v',
    'HLT_HIEle20Gsf_v',
    'HLT_HIL2SingleMu7_v',
)

trigger_list_data_2024_ppRef_skimmed = cms.vstring(
    'HLT_PPRefGEDPhoton30_v',
    'HLT_PPRefEle20Gsf_v',
    'HLT_PPRefL2SingleMu7_v',
)

trigger_list_data_2025_pO_skimmed = cms.vstring(
    'HLT_OxyL1SingleMuOpen_v',
    'HLT_OxyL1SingleMu0_v',
    'HLT_OxyL1SingleEG10_v',
    'HLT_OxyL1SingleEG15_v',
)
trigger_list_data_2025_OO_skimmed = cms.vstring(trigger_list_data_2025_pO_skimmed)
trigger_list_data_2025_NeNe_skimmed = cms.vstring(trigger_list_data_2025_pO_skimmed)

trigger_list_data_2025_skimmed = cms.vstring(
    'HLT_HIGEDPhoton10_v',
    'HLT_HIL2SingleMu7_v',
)

trigger_list_data_2026_skimmed = cms.vstring(
    'HLT_HIGEDPhoton10_v',
    'HLT_HIL2SingleMu7_v',
)
