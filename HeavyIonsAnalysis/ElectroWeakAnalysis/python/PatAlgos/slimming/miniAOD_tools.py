import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.slimming.miniAOD_tools import *

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import overrideJEC_pPb8TeV


def includeMETFilters(process):
    # Add New Bad Muon MET Filters
    process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
    process.load("RecoMET.METFilters.badGlobalMuonTaggersAOD_cff")
    process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")

    # individual filters
    process.HBHENoiseFilterRun1                       = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun1') )
    process.HBHENoiseFilterRun2Loose                  = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Loose') )
    process.HBHENoiseFilterRun2Tight                  = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight') )
    process.BadChargedCandidateBooleanFilter          = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('BadChargedCandidateFilter') )
    process.BadPFMuonBooleanFilter                    = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('BadPFMuonFilter') )
    process.BadChargedCandidateSummer16BooleanFilter  = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('BadChargedCandidateSummer16Filter') )
    process.BadPFMuonSummer16BooleanFilter            = process.HBHENoiseFilter.clone( inputLabel = cms.InputTag('BadPFMuonSummer16Filter') )
    process.badTrackerMuonProducer  = process.selectedPatMuons.clone( cut = cms.string("isTrackerMuon() && pt()>20. && !( numberOfMatchedStations()>0 && ( innerTrack().numberOfValidHits()>=10 || ( innerTrack().numberOfValidHits()>=7 && innerTrack().numberOfLostHits()==0 ) ) )") )
    process.badTrackerMuonFilter = cms.EDFilter("PATCandViewCountFilter", minNumber = cms.uint32(0), maxNumber = cms.uint32(0), src = cms.InputTag("badTrackerMuonProducer") )

    process.Flag_collisionEventSelectionPA = cms.Path( process.collisionEventSelectionPA )
    process.Flag_collisionEventSelectionPA_rejectPU = cms.Path( process.collisionEventSelectionPA_rejectPU )
    process.Flag_HBHENoiseFilterRun1 = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilterRun1)
    process.Flag_HBHENoiseFilterRun2Loose = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilterRun2Loose)
    process.Flag_HBHENoiseFilterRun2Tight = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilterRun2Tight)
    process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter * process.BadChargedCandidateBooleanFilter)
    process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter * process.BadPFMuonBooleanFilter)
    process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter * process.BadChargedCandidateSummer16BooleanFilter)
    process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter * process.BadPFMuonSummer16BooleanFilter)
    process.Flag_noBadMuons = cms.Path(process.noBadGlobalMuons)
    process.Flag_badMuons = cms.Path(~process.badGlobalMuonTagger)
    process.Flag_badTrackerMuons = cms.Path(process.badTrackerMuonProducer * process.badTrackerMuonFilter)
    process.Flag_duplicateMuons = cms.Path(~process.cloneGlobalMuonTagger)

    metFilterPaths = [ process.Flag_collisionEventSelectionPA , process.Flag_collisionEventSelectionPA_rejectPU , process.Flag_HBHENoiseFilter , 
                       process.Flag_HBHENoiseFilterRun1 , process.Flag_HBHENoiseFilterRun2Loose , process.Flag_HBHENoiseFilterRun2Tight,
                       process.Flag_HBHENoiseIsoFilter, process.Flag_CSCTightHaloFilter, process.Flag_CSCTightHaloTrkMuUnvetoFilter , process.Flag_CSCTightHalo2015Filter, 
                       process.Flag_globalTightHalo2016Filter , process.Flag_globalSuperTightHalo2016Filter , process.Flag_HcalStripHaloFilter , process.Flag_hcalLaserEventFilter, 
                       process.Flag_EcalDeadCellTriggerPrimitiveFilter , process.Flag_EcalDeadCellBoundaryEnergyFilter , process.Flag_goodVertices , process.Flag_trackingFailureFilter,
                       process.Flag_eeBadScFilter , process.Flag_chargedHadronTrackResolutionFilter, process.Flag_muonBadTrackFilter , process.Flag_BadChargedCandidateFilter , 
                       process.Flag_BadPFMuonFilter , process.Flag_BadChargedCandidateSummer16Filter, process.Flag_BadPFMuonSummer16Filter , process.Flag_noBadMuons, process.Flag_badMuons, 
                       process.Flag_duplicateMuons , process.Flag_ecalLaserCorrFilter , process.Flag_badTrackerMuons ]

    for P in metFilterPaths:
        process.schedule.insert(0, P)


def applyTriggerMatching(process):
    #
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTriggerStandAlone , switchOnTriggerMatchEmbedding
    switchOnTriggerStandAlone( process, outputModule = '' )
    #process.patTrigger.packTriggerPathNames = cms.bool(True)
    # Trigger Matching for HLT
    process.triggerMatchDRDPt = cms.EDProducer("PATTriggerMatcherDRDPtLessByR", src = cms.InputTag( "" ), matched = cms.InputTag( 'selectedPatTrigger' ),  matchedCuts = cms.string(""), maxDPtRel = cms.double( 0.5 ), maxDeltaR = cms.double( 0.5 ), resolveAmbiguities = cms.bool( False ), resolveByMatchQuality = cms.bool( True ))
    process.triggerMatchDR    = cms.EDProducer("PATTriggerMatcherDRLessByR", src = cms.InputTag( "" ), matched = cms.InputTag( 'selectedPatTrigger' ),  matchedCuts = cms.string(""), maxDeltaR = cms.double( 0.5 ), resolveAmbiguities = cms.bool( False ), resolveByMatchQuality = cms.bool( True ))
    # Muon HLT Trigger Matching
    process.muonMatchHLTL2   = process.triggerMatchDRDPt.clone(src = cms.InputTag( "patMuons" ), matchedCuts = cms.string('coll("hltL2MuonCandidates")'),   maxDeltaR = 0.3, maxDPtRel = 10.0)  
    process.muonMatchHLTL3   = process.triggerMatchDRDPt.clone(src = cms.InputTag( "patMuons" ), matchedCuts = cms.string('coll("hltL3MuonCandidates")'),   maxDeltaR = 0.1, maxDPtRel = 10.0)  
    process.muonMatchHLTHIL3 = process.triggerMatchDRDPt.clone(src = cms.InputTag( "patMuons" ), matchedCuts = cms.string('coll("hltHIL3MuonCandidates")'), maxDeltaR = 0.1, maxDPtRel = 10.0)
    process.patMuonsWithTriggers = cms.EDProducer("PATTriggerMatchMuonEmbedder", src = cms.InputTag( "patMuons" ),  matches = cms.VInputTag('muonMatchHLTL2', 'muonMatchHLTL3', 'muonMatchHLTHIL3') )
    # Electron HLT Trigger Matching
    #process.electronMatchHLT = process.triggerMatchDRDPt.clone(src = cms.InputTag( "patElectrons" ), matchedCuts = cms.string('coll("hltEgammaCandidates")'),   maxDeltaR = 0.3, maxDPtRel = 10.0)  #FIXME
    #process.patElectronsWithTriggers = cms.EDProducer("PATTriggerMatchElectronEmbedder", src = cms.InputTag( "patElectrons" ),  matches = cms.VInputTag('electronMatchHLT') )
    # Photon HLT Trigger Matching
    #process.photonMatchHLT   = process.triggerMatchDRDPt.clone(src = cms.InputTag( "patPhotons" ), matchedCuts = cms.string('coll("hltEgammaCandidates")'),                maxDeltaR = 0.3, maxDPtRel = 10.0)  #FIXME
    #process.photonMatchHIHLT = process.triggerMatchDRDPt.clone(src = cms.InputTag( "patPhotons" ), matchedCuts = cms.string('coll("hltRecoHIEcalWithCleaningCandidate")'), maxDeltaR = 0.3, maxDPtRel = 10.0)  #FIXME
    #process.patPhotonsWithTriggers = cms.EDProducer("PATTriggerMatchPhotonEmbedder", src = cms.InputTag( "patPhotons" ),  matches = cms.VInputTag('photonMatchHLT' , 'photonMatchHIHLT') )
    # Jet HLT Trigger Matching
    #process.jetMatchHLTPFJets             = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPAAK4PFJetsCorrected")'),                            maxDeltaR = 0.4)  
    #process.jetMatchHLTPFBJets            = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPACombinedSecondaryVertexBJetTagsPF")'),             maxDeltaR = 0.4)  
    #process.jetMatchHLTPFBjetsCommTrack   = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltBTagPFCSVp016SingleWithMatching40CommonTracking")'), maxDeltaR = 0.4) 
    #process.jetMatchHLTCaloJets           = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPAAK4CaloJetsCorrectedIDPassed")'),                  maxDeltaR = 0.4) 
    #process.jetMatchHLTPFCaloJets30Eta2p1 = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPAAK4PFJetsCorrectedMatchedToCaloJets30Eta2p1")'),   maxDeltaR = 0.4) 
    #process.jetMatchHLTPFCaloJets15Eta5p1 = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPAAK4PFJetsCorrectedMatchedToCaloJets15Eta5p1")'),   maxDeltaR = 0.4)  
    #process.jetMatchHLTPFCaloJets30Eta5p1 = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPAAK4PFJetsCorrectedMatchedToCaloJets30Eta5p1")'),   maxDeltaR = 0.4)  
    #process.jetMatchHLTPFCaloJets50Eta5p1 = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPAAK4PFJetsCorrectedMatchedToCaloJets50Eta5p1")'),   maxDeltaR = 0.4)  
    #process.jetMatchHLTCaloBJets          = process.triggerMatchDR.clone(src = cms.InputTag( "patJets" ), matchedCuts = cms.string('coll("hltPABLifetimeL3FilterCSVCaloJet40Eta2p1")'),           maxDeltaR = 0.4)
    #process.patJetsWithTriggers = cms.EDProducer("PATTriggerMatchJetEmbedder", src = cms.InputTag( "patJets" ),  matches = cms.VInputTag('jetMatchHLTPFJets', 'jetMatchHLTCaloJets', 'jetMatchHLTPFCaloJets30Eta2p1', 'jetMatchHLTPFCaloJets15Eta5p1', 'jetMatchHLTPFCaloJets30Eta5p1', 'jetMatchHLTPFCaloJets50Eta5p1', 'jetMatchHLTCaloBJets', 'jetMatchHLTPFBJets', 'jetMatchHLTPFBjetsCommTrack') )
    # Keep HLT Trigger Matching Objects
    process.MicroEventContent.outputCommands.append( 'keep *_*WithTriggers_*_*' )



def miniAOD_ForHiEWQ_customizeCommon(process, isData):
    # Add JEC for pPb 8 TeV
    overrideJEC_pPb8TeV(process)
    # Add MET Filters 
    includeMETFilters(process)
    # Apply trigger matching
    applyTriggerMatching(process)

    process.patMuons.isoDeposits = cms.PSet()
    process.patElectrons.isoDeposits = cms.PSet()
    process.patTaus.isoDeposits = cms.PSet()
    process.patPhotons.isoDeposits = cms.PSet()
    #
    process.patMETs.srcJetResPhi = cms.string('AK4PF_phi')
    process.patMETs.srcJetResPt = cms.string('AK4PF_pt')
    process.patMETs.srcJetSF = cms.string('AK4PF_offline')
    #
    process.patMuons.embedTrack         = True  # used for IDs
    process.patMuons.embedCombinedMuon  = True  # used for IDs
    process.patMuons.embedMuonBestTrack = True  # used for IDs
    process.patMuons.embedStandAloneMuon = True # maybe?
    process.patMuons.embedPickyMuon = False   # no, use best track
    process.patMuons.embedTpfmsMuon = False   # no, use best track
    process.patMuons.embedDytMuon   = False   # no, use best track
    #
    # disable embedding of electron and photon associated objects already stored by the ReducedEGProducer
    process.patElectrons.embedGsfElectronCore = False  ## process.patElectrons.embed in AOD externally stored gsf electron core
    process.patElectrons.embedSuperCluster    = False  ## process.patElectrons.embed in AOD externally stored supercluster
    process.patElectrons.embedPflowSuperCluster         = False  ## process.patElectrons.embed in AOD externally stored supercluster
    process.patElectrons.embedSeedCluster               = False  ## process.patElectrons.embed in AOD externally stored the electron's seedcluster
    process.patElectrons.embedBasicClusters             = False  ## process.patElectrons.embed in AOD externally stored the electron's basic clusters
    process.patElectrons.embedPreshowerClusters         = False  ## process.patElectrons.embed in AOD externally stored the electron's preshower clusters
    process.patElectrons.embedPflowBasicClusters        = False  ## process.patElectrons.embed in AOD externally stored the electron's pflow basic clusters
    process.patElectrons.embedPflowPreshowerClusters    = False  ## process.patElectrons.embed in AOD externally stored the electron's pflow preshower clusters
    process.patElectrons.embedRecHits         = False  ## process.patElectrons.embed in AOD externally stored the RecHits - can be called from the PATElectronProducer
    process.patElectrons.electronSource = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.patElectrons.electronIDSources = cms.PSet(
            # configure many IDs as InputTag <someName> = <someTag> you
            # can comment out those you don't want to save some disk space
            eidRobustLoose      = cms.InputTag("reducedEgamma","eidRobustLoose"),
            eidRobustTight      = cms.InputTag("reducedEgamma","eidRobustTight"),
            eidLoose            = cms.InputTag("reducedEgamma","eidLoose"),
            eidTight            = cms.InputTag("reducedEgamma","eidTight"),
            eidRobustHighEnergy = cms.InputTag("reducedEgamma","eidRobustHighEnergy"),
            )
    process.patElectrons.addPFClusterIso = cms.bool(True)
    process.patElectrons.ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "eleEcalPFClusIso")
    process.patElectrons.hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "eleHcalPFClusIso")

    process.elPFIsoDepositChargedPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositChargedAllPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositNeutralPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositGammaPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositPUPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    #
    process.patPhotons.embedSuperCluster = False ## whether to process.patPhotons.embed in AOD externally stored supercluster
    process.patPhotons.embedSeedCluster               = False  ## process.patPhotons.embed in AOD externally stored the photon's seedcluster
    process.patPhotons.embedBasicClusters             = False  ## process.patPhotons.embed in AOD externally stored the photon's basic clusters
    process.patPhotons.embedPreshowerClusters         = False  ## process.patPhotons.embed in AOD externally stored the photon's preshower clusters
    process.patPhotons.embedRecHits         = False  ## process.patPhotons.embed in AOD externally stored the RecHits - can be called from the PATPhotonProducer
    process.patPhotons.addPFClusterIso = cms.bool(True)
    process.patPhotons.ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "phoEcalPFClusIso")
    process.patPhotons.hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "phoHcalPFClusIso")
    process.patPhotons.photonSource = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.patPhotons.electronSource = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.patPhotons.photonIDSources = cms.PSet(
                PhotonCutBasedIDLoose = cms.InputTag('reducedEgamma',
                                                      'PhotonCutBasedIDLoose'),
                PhotonCutBasedIDTight = cms.InputTag('reducedEgamma',
                                                      'PhotonCutBasedIDTight')
              )

    process.phPFIsoDepositChargedPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositChargedAllPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositNeutralPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositGammaPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositPUPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    #

    if isData:
        process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_DataPPb")
    else:
        process.load("HeavyIonsAnalysis.JetAnalysis.FullJetSequence_JECPPb")

        
    # make ak4PF jets the baseline jets
    process.patJets = process.ak4PFpatJetsWithBtagging.clone()

    #
    process.selectedPatJets.cut = cms.string("pt > 10")
    process.selectedPatMuons.cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))")
    process.selectedPatElectrons.cut = cms.string("")
    process.selectedPatTaus.cut = cms.string("pt > 18. && tauID('decayModeFindingNewDMs')> 0.5")
    process.selectedPatPhotons.cut = cms.string("")

    #
    # apply type I + other PFMEt corrections to pat::MET object
    # and estimate systematic uncertainties on MET

    # pfMET from ak4PF jets
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncForMiniAODProduction
    runMetCorAndUncForMiniAODProduction(process, metType="PF",
                                        jetCollUnskimmed="ak4PFpatJetsWithBtagging",
                                        jetFlavor="AK4PF_offline",
                                        postfix=""
                                        )

    #caloMET computation
    from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
    addMETCollection(process,
                     labelName = "patCaloMet",
                     metSource = "caloMetM"
                     )

    
    #noHF pfMET =========
    process.noHFCands = cms.EDFilter("GenericPFCandidateSelector",
                                     src=cms.InputTag("particleFlow"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )
    runMetCorAndUncForMiniAODProduction(process,
                                        pfCandColl=cms.InputTag("noHFCands"),
                                        recoMetFromPFCs=True, #needed for HF removal
                                        jetCollUnskimmed="ak4PFpatJetsWithBtagging",
                                        jetFlavor="AK4PF_offline",
                                        jetSelection="pt>15 && abs(eta)<3.",
                                        postfix="NoHF"
                                        )
    process.corrPfMetType1NoHF.src = cms.InputTag("ak4PFJets")

    process.load('PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi')
    process.slimmedMETsNoHF = process.slimmedMETs.clone()
    process.slimmedMETsNoHF.src = cms.InputTag("patMETsNoHF")
    process.slimmedMETsNoHF.rawVariation =  cms.InputTag("patPFMetNoHF")
    process.slimmedMETsNoHF.t1Uncertainties = cms.InputTag("patPFMetT1%sNoHF") 
    process.slimmedMETsNoHF.t01Variation = cms.InputTag("patPFMetT0pcT1NoHF")
    process.slimmedMETsNoHF.t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%sNoHF")
    process.slimmedMETsNoHF.tXYUncForRaw = cms.InputTag("patPFMetTxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT1 = cms.InputTag("patPFMetT1TxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT01 = cms.InputTag("patPFMetT0pcT1TxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT1Smear = cms.InputTag("patPFMetT1SmearTxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT01Smear = cms.InputTag("patPFMetT0pcT1SmearTxyNoHF")
    del process.slimmedMETsNoHF.caloMET
    # ================== NoHF pfMET
    #keep this after all addJetCollections otherwise it will attempt computing them also for stuf with no taginfos
    #Some useful BTAG vars
    if not hasattr( process, 'pfImpactParameterTagInfos' ):
        process.load('RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi')
    if not hasattr( process, 'pfSecondaryVertexTagInfos' ):
        process.load('RecoBTag.SecondaryVertex.pfSecondaryVertexTagInfos_cfi')
    process.patJets.userData.userFunctions = cms.vstring(
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).p4.M):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).numberOfSourceCandidatePtrs):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").flightDistance(0).value):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").flightDistance(0).significance):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).p4.x):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).p4.y):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).p4.z):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).vertex.x):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).vertex.y):(0)',
    '?(tagInfoCandSecondaryVertex("pfSecondaryVertex").nVertices()>0)?(tagInfoCandSecondaryVertex("pfSecondaryVertex").secondaryVertex(0).vertex.z):(0)',
    )
    process.patJets.userData.userFunctionLabels = cms.vstring('vtxMass','vtxNtracks','vtx3DVal','vtx3DSig','vtxPx','vtxPy','vtxPz','vtxPosX','vtxPosY','vtxPosZ')
    process.patJets.tagInfoSources = cms.VInputTag(cms.InputTag("pfSecondaryVertexTagInfos"))
    process.patJets.addTagInfos = cms.bool(True)
    
    #EGM object modifications
    from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma_modifications
    process.slimmedElectrons.modifierConfig.modifications = egamma_modifications
    process.slimmedPhotons.modifierConfig.modifications   = egamma_modifications

    #VID Electron IDs
    electron_ids = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
    switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
    process.egmGsfElectronIDs.physicsObjectSrc = \
        cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.electronMVAValueMapProducer.src = \
        cms.InputTag('reducedEgamma','reducedGedGsfElectrons')
    process.electronRegressionValueMapProducer.src = \
        cms.InputTag('reducedEgamma','reducedGedGsfElectrons')
    for idmod in electron_ids:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection,None,False)

    #VID Photon IDs
    photon_ids = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                  'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
    switchOnVIDPhotonIdProducer(process,DataFormat.AOD)
    process.egmPhotonIsolation.srcToIsolate = \
        cms.InputTag("reducedEgamma","reducedGedPhotons")  
    for iPSet in process.egmPhotonIsolation.isolationConeDefinitions:
        iPSet.particleBasedIsolation = cms.InputTag("reducedEgamma","reducedPhotonPfCandMap") 
 
    process.egmPhotonIDs.physicsObjectSrc = \
        cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.photonIDValueMapProducer.src = \
        cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.photonRegressionValueMapProducer.src = \
        cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.photonIDValueMapProducer.particleBasedIsolation = \
        cms.InputTag("reducedEgamma","reducedPhotonPfCandMap")
    process.photonMVAValueMapProducer.src = \
        cms.InputTag('reducedEgamma','reducedGedPhotons')
    for idmod in photon_ids:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection,None,False)


    # REMOVE PUPPI FROM HELL !!!!!!!!!!!!
    process.packedPFCandidates.PuppiSrc = ''
    process.packedPFCandidates.PuppiNoLepSrc = ''
    # REMOVE AK8 JETS !!!!!!!!!!!!
    del process.slimmedJetsAK8
    #---------------------------------------------------------------------------
    #Ignore  Boosted Subjets taus
    #from RecoTauTag.Configuration.boostedHPSPFTaus_cfi import addBoostedTaus
    #addBoostedTaus(process)
    #---------------------------------------------------------------------------
    #Ignore noHF pfMET =========
    ## Ignore Legacy tight b-tag track selection

    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag , massSearchReplaceParam
    listSeq = [cms.Sequence( process.patCaloMet + process.ak4PFL1OffsetCorrector ) , process.pfBTagging , process.patCandidates , process.ak4PFL1FastL2L3ResidualCorrectorChain, process.patMETCorrections, process.patMETCorrectionsNoHF, process.fullPatMetSequence , process.fullPatMetSequenceNoHF ]
    for seq in listSeq:
        massSearchReplaceAnyInputTag(seq,cms.InputTag("ak4PFJetsCHS"),cms.InputTag("ak4PFJets"), skipLabelTest=True)
        massSearchReplaceAnyInputTag(seq,cms.InputTag("AK4PFchs"),cms.InputTag("AK4PF_offline"), skipLabelTest=True)
        massSearchReplaceAnyInputTag(seq,cms.InputTag("AK4PF"),cms.InputTag("AK4PF_offline"), skipLabelTest=True)
        massSearchReplaceParam(seq,"algorithm", cms.string('AK4PF'), cms.string('AK4PF_offline'))
        massSearchReplaceParam(seq,"algorithm", cms.string('AK4PFchs'), cms.string('AK4PF_offline'))
        massSearchReplaceParam(seq,"payload", cms.string('AK4PFchs'), cms.string('AK4PF_offline'))
        massSearchReplaceParam(seq,"payload", cms.string('AK4PF'), cms.string('AK4PF_offline'))
        massSearchReplaceAnyInputTag(seq,cms.InputTag("L2L3Residual"),cms.InputTag("L3Absolute"), skipLabelTest=True)
        massSearchReplaceAnyInputTag(seq,cms.InputTag("L1FastJet"),cms.InputTag(''), skipLabelTest=True)
        massSearchReplaceAnyInputTag(seq,cms.InputTag("ak4PFL1FastL2L3ResidualCorrector"),cms.InputTag('ak4PFL1FastL2L3Corrector'), skipLabelTest=True)
        massSearchReplaceParam(seq,"addResidualJES", cms.bool(True), cms.bool(False))
        massSearchReplaceAnyInputTag(seq,cms.InputTag("AK4PFchs_pt"),cms.InputTag("AK4PF_pt"), skipLabelTest=True)
        massSearchReplaceAnyInputTag(seq,cms.InputTag("AK4PFchs_phi"),cms.InputTag("AK4PF_phi"), skipLabelTest=True)
        massSearchReplaceParam(seq,"srcJetResPt",cms.string('AK4PFchs_pt'),cms.string('AK4PF_pt'))
        massSearchReplaceParam(seq,"algopt",cms.string('AK4PFchs_pt'),cms.string('AK4PF_pt'))
        massSearchReplaceParam(seq,"srcJetResPhi",cms.string('AK4PFchs_phi'),cms.string('AK4PF_phi'))
        massSearchReplaceParam(seq,"jetCorrPayloadName", cms.string('AK4PF_offline'), cms.string("AK4PF"))
        massSearchReplaceParam(seq,"srcJetSF", cms.string('AK4PFchs'), cms.string('AK4PF'))
        massSearchReplaceParam(seq,"algo", cms.string('AK4PF'), cms.string('AK4PF_offline'))
        if isData: massSearchReplaceParam(seq,"algo", cms.string('AK4PFchs'), cms.string('AK4PF_offline'))
        else:      massSearchReplaceParam(seq,"algo", cms.string('AK4PFchs'), cms.string('AK4PF'))

    del process.towerMaker

def miniAOD_ForHiEWQ_customizeMC(process):
    process.load('PhysicsTools.PatAlgos.slimming.slimmedAddPileupInfo_cfi')
    process.muonMatch.matched = "prunedGenParticles"
    process.electronMatch.matched = "prunedGenParticles"
    process.electronMatch.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.photonMatch.matched = "prunedGenParticles"
    process.photonMatch.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.tauMatch.matched = "prunedGenParticles"
    process.tauGenJets.GenParticles = "prunedGenParticles"
    process.patJetPartons.particles = "prunedGenParticles"
    process.patJetPartonMatch.matched = "prunedGenParticles"
    process.patJetPartonMatch.mcStatus = [ 3, 23 ]
    process.patJetGenJetMatch.matched = "slimmedGenJets"
    process.patElectrons.embedGenMatch = False
    process.patPhotons.embedGenMatch = False
    process.patTaus.embedGenMatch = False
    process.patJets.embedGenPartonMatch = False
    #also jet flavour must be switched
    process.patJetFlavourAssociation.rParam = 0.4
    # For Onia Tree
    process.patMuons.embedGenMatch = True
    process.muonMatch.resolveByMatchQuality = True # NEW WITH RESPECT TO BEFORE

    

def miniAOD_ForHiEWQ_customizeAllData(process):
    miniAOD_ForHiEWQ_customizeCommon(process, isData=True)
    miniAOD_customizeData(process)
    return process

def miniAOD_ForHiEWQ_customizeAllMC(process):
    miniAOD_ForHiEWQ_customizeCommon(process, isData=False)
    miniAOD_ForHiEWQ_customizeMC(process)
    return process
