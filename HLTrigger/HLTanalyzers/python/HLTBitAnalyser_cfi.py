import FWCore.ParameterSet.Config as cms

hltbitanalysis = cms.EDAnalyzer("HLTBitAnalyzer",
    ### General Settings
    UseL1Stage2                     = cms.untracked.bool(True),
    stageL1Trigger                  = cms.uint32(2),
    getPrescales                    = cms.untracked.bool(False),
    getL1InfoFromEventSetup         = cms.untracked.bool(False),

    ### L1 Legacy and Stage 1 objects
    l1GctHFBitCounts                = cms.InputTag("hltGctDigis"),
    l1GctHFRingSums                 = cms.InputTag("hltGctDigis"),
    l1GtReadoutRecord               = cms.InputTag("hltGtDigis::HLT"),
    l1extramc                       = cms.string('hltL1extraParticles'),
    l1extramu                       = cms.string('hltL1extraParticles'),

    ### L1 Stage 2 objects
    l1tAlgBlkInputTag               = cms.InputTag("hltGtStage2Digis"),  # Needed, fix bug of GlobalAlgBlk uninitialized token
    l1tExtBlkInputTag               = cms.InputTag("hltGtStage2Digis"),
    gObjectMapRecord                = cms.InputTag("hltGtStage2ObjectMap"),
    gmtStage2Digis                  = cms.string("hltGtStage2Digis"),
    caloStage2Digis                 = cms.string("hltGtStage2Digis"),

    ### HLT
    hltresults                      = cms.InputTag("TriggerResults::HLT"),
    HLTProcessName                  = cms.string("HLT"),

    ### GEN objects
    mctruth                         = cms.InputTag("genParticles::HLT"),
    genEventInfo                    = cms.InputTag("generator::SIM"),

    ### SIM objects
    simhits                         = cms.InputTag("g4SimHits"),

    ## reco vertices
    OfflinePrimaryVertices0         = cms.InputTag('offlinePrimaryVertices'),

    ### Run parameters
    RunParameters = cms.PSet(
        HistogramFile = cms.untracked.string('hltbitanalysis.root'),
        isData         = cms.untracked.bool(False),
        Monte          = cms.bool(True),
        GenTracks      = cms.bool(False),
    )

)
