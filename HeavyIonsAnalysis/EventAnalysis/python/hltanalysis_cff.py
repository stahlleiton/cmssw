import FWCore.ParameterSet.Config as cms

from HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi import *

hltanalysis = hltbitanalysis.clone(
    hltresults = cms.InputTag("TriggerResults::HLT"),
    UseTFileService = cms.untracked.bool(True),
    RunParameters = cms.PSet(
        isData = cms.untracked.bool(True)),

    OfflinePrimaryVertices0 = cms.InputTag(""),

    dummyBranches = cms.untracked.vstring()
    )




skimanalysis = cms.EDAnalyzer("FilterAnalyzer",
                              hltresults = cms.InputTag("TriggerResults","","HiForest"),
                              superFilters = cms.vstring("")
)




