import FWCore.ParameterSet.Config as cms

DisplacedL2DimuonPSet = cms.PSet(
    hltPathsToCheck = cms.vstring(
        "HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v", # Run2 
        "HLT_L2DoubleMu23_NoVertex_v", # Run2 
        "HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_v" # Run2 
        ),
    #recMuonLabel  = cms.InputTag("muons"),
    recMuonTrkLabel  = cms.InputTag("refittedStandAloneMuons"),
    # -- Analysis specific cuts
    minCandidates = cms.uint32(2),
    # -- Analysis specific binnings
    parametersDxy      = cms.vdouble(50, -50, 50),
    parametersTurnOn = cms.vdouble( 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                                    60, 70, 80, 100
                                   ),
    dropPt3 = cms.bool(True),
    )
