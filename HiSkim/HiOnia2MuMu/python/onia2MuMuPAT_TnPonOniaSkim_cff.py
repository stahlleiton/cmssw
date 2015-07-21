import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.helpers import *

def tnpOnOniaSkim(process, GlobalTag, MC=False, HLT='HLT', Filter=True):
    # Setup the process

    # PbPb
    IN_ACCEPTANCE = '( (abs(eta)<1.0 && pt>=3.4) || (1.0<=abs(eta)<1.5 && pt>=5.8-2.4*abs(eta)) || (1.5<=abs(eta)<2.4 && pt>=3.3667-7.0/9.0*abs(eta)) )'
    # several selection cuts
    TRACK_CUTS    = "isTrackerMuon && track.numberOfValidHits > 10 && track.normalizedChi2 < 4 && track.hitPattern.pixelLayersWithMeasurement > 0 && abs(dB) < 3 && abs(track.dz) < 15  && muonID('TrackerMuonArbitrated')"
    GLB_CUTS      = "isGlobalMuon && globalTrack.normalizedChi2 < 20"
    QUALITY_CUTS  =  "(" + GLB_CUTS + ' && ' + TRACK_CUTS + ")"
    
    # Make Tag and Probe pairs for efficiency measurements
    process.tagMuonsSglTrg = cms.EDFilter("PATMuonSelector",
                                          src = cms.InputTag("patMuonsWithTrigger"),
                                          cut = cms.string(QUALITY_CUTS + ' && ' + IN_ACCEPTANCE + " && (!triggerObjectMatchesByPath('HLT_HIL2Mu3_NHitQ_v*').empty() || !triggerObjectMatchesByPath('HLT_HIL2Mu7_v*').empty() || !triggerObjectMatchesByPath('HLT_HIL2Mu15_v*').empty())")
                                          )

    # produce patMuons that use the STA momentum information
    process.patMuonsWithTriggerSta = cms.EDProducer("RedefineMuonP4FromTrackPAT",
                                                    src   = cms.InputTag("patMuonsWithTrigger"),
                                                    track = cms.string("outer")
                                                    )

    # tracking efficiency: probe (STA), and tag-probe pair
    process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
                                         src = cms.InputTag("patMuonsWithTriggerSta"),
                                         cut = cms.string("outerTrack.isNonnull")
                                         )

    process.tpPairsSta = cms.EDProducer("CandViewShallowCloneCombiner",
                                        cut = cms.string('1.0 < mass < 5.0'),
                                        decay = cms.string('tagMuonsSglTrg@+ probeMuonsSta@-')
                                        )

    # trigger efficiency: probe (muons that pass all the analyis cuts), and tag-probe pair
    process.probeMuons = cms.EDFilter("PATMuonSelector",
                                      src = cms.InputTag("patMuonsWithTrigger"),
                                      cut = cms.string(QUALITY_CUTS + ' && ' + IN_ACCEPTANCE)
                                      )

    process.tpPairsTrig = cms.EDProducer("CandViewShallowCloneCombiner",
                                     cut = cms.string('2.6 < mass < 4.0'),
                                     decay = cms.string('tagMuonsSglTrg@+ probeMuons@-')
                                     )

    # muID efficiency: probe (global muons), and tag-probe pair
    process.probeMuonsTrk = cms.EDFilter("PATMuonSelector",
                                         src = cms.InputTag("patMuonsWithTrigger"),
                                         cut = cms.string("isGlobalMuon && " + IN_ACCEPTANCE)
                                         )

    process.tpPairsTracks = cms.EDProducer("CandViewShallowCloneCombiner",
                                           cut = cms.string('2. < mass < 4.0'),
                                           decay = cms.string('tagMuonsSglTrg@+ probeMuonsTrk@-')
                                           )

    # check if there is at least one Tag and Probe pair
    process.tpPairsStaFilter = cms.EDFilter("CandViewCountFilter",
                                            src = cms.InputTag('tpPairsSta'),
                                            minNumber = cms.uint32(1),
                                            )

    process.tpPairsTrigFilter = cms.EDFilter("CandViewCountFilter",
                                             src = cms.InputTag('tpPairsTrig'),
                                             minNumber = cms.uint32(1),
                                             )

    process.tpPairsTracksFilter = cms.EDFilter("CandViewCountFilter",
                                               src = cms.InputTag('tpPairsTracks'),
                                               minNumber = cms.uint32(1),
                                               )
    
    ##### the Tag and Probe sequences: make the probe, the tag-probe pairs, and make sure have at least 1 per/event
    process.tnpSta = cms.Sequence(
        process.probeMuonsSta *
        process.tpPairsSta *
        process.tpPairsStaFilter
        )

    process.tnpTrig = cms.Sequence(
        process.probeMuons *
        process.tpPairsTrig *
        process.tpPairsTrigFilter
        )

    process.tnpMuID = cms.Sequence(
        process.probeMuonsTrk *
        process.tpPairsTracks *
        process.tpPairsTracksFilter
        )

    #####  create first the tag, then proceed with probe and pairs
    # inner track reco efficiency: 
    process.TagAndProbeSta = cms.Path(
        process.tagMuonsSglTrg *
        process.tnpSta
    )

    # muon reco and ID efficiency
    process.TagAndProbeMuID = cms.Path(
        process.tagMuonsSglTrg *
        process.tnpMuID
    )

    # muon trigger efficiency
    process.TagAndProbeTrig = cms.Path(
        process.tagMuonsSglTrg *
        process.tnpTrig
    )
    
    if MC:
        process.tagMuonsSglTrgMCMatch = process.muonMatch.clone(src = "tagMuonsSglTrg") # tag
        process.probeMuonsStaMCMatch  = process.muonMatch.clone(src = "probeMuonsSta") # inner tracking eff
        process.probeMuonsTrkMCMatch  = process.muonMatch.clone(src = "probeMuonsTrk") # Muon reco and ID eff
        process.probeMuonsMCMatch     = process.muonMatch.clone(src = "probeMuons") # muon trigger eff
        process.TagAndProbeSta.replace(process.tpPairsSta, process.tagMuonsSglTrgMCMatch * process.probeMuonsStaMCMatch * process.tpPairsSta)
        process.TagAndProbeMuID.replace(process.tpPairsTracks, process.tagMuonsSglTrgMCMatch * process.probeMuonsTrkMCMatch * process.tpPairsTracks)
        process.TagAndProbeTrig.replace(process.tpPairsTrig, process.tagMuonsSglTrgMCMatch * process.probeMuonsMCMatch * process.tpPairsTrig)
    
    # output
    process.load('Configuration.EventContent.EventContent_cff')
    process.load("Configuration.EventContent.EventContentHeavyIons_cff")

    process.outTnP = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('tnp.root'),
        outputCommands = cms.untracked.vstring('keep *'),
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('TagAndProbeSta','TagAndProbeMuID','TagAndProbeTrig') ) if Filter else cms.untracked.PSet()
    )

    process.e = cms.EndPath(process.outTnP)

