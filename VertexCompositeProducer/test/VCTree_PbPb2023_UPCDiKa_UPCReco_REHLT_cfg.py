import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023_UPC)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.numberOfThreads=cms.untracked.uint32(1)

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/965305ee-df8f-493c-b380-f477a5c10547.root"),
    secondaryFileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/19368c6e-3485-4b23-8950-2e5d542b186c.root"),
    lumisToProcess = cms.untracked.VLuminosityBlockRange('375455:162')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v7')


## ##############################################################################################################################
## Variables Production #########################################################################################################

#* Set ZDC information
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    authenticationMethod = cms.untracked.uint32(1)
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("VertexCompositeAnalysis/VertexCompositeProducer/data/emap_2023_newZDC_v3.txt")))
)

#* cent_seq: Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.cent_seq = cms.Sequence(process.centralityBin)

#* Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

# DiKa selection
kaonSelection = cms.string("")#(pt > 0.0 && abs(eta) < 3.0) && quality(\"highPurity\")")
kaonFinalSelection = cms.string("")#abs(userFloat(\"dzSig\"))<3.0 && abs(userFloat(\"dxySig\"))<3.0")
diKaSelection = cms.string("charge==0")
process.diKa = generalParticles.clone(
    pdgId = cms.uint32(333),
    primaryVertices = "primaryVertexRecoveryForUPC",
    preSelection = diKaSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.uint32(321), charge = cms.int32(+1), selection = kaonSelection, finalSelection = kaonFinalSelection),
        cms.PSet(pdgId = cms.uint32(321), charge = cms.int32(-1), selection = kaonSelection, finalSelection = kaonFinalSelection),
    ]),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
    # dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2', 'energyLossProducer:energyLossAllHits', 'energyLossProducer:energyLossPixHits', 'energyLossProducer:energyLossStrHits')
)
process.oneDiKa = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("diKa"), minNumber = cms.uint32(1))

# Add diKa event selection
process.twoTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("generalTracks"), minNumber = cms.uint32(2))
process.hpTracks = cms.EDFilter("TrackSelector", src = cms.InputTag("generalTracks"), cut = cms.string("quality(\"highPurity\")"))
process.hpCands = cms.EDProducer("ChargedCandidateProducer", src = cms.InputTag("hpTracks"), particleType = cms.string('pi+'))
process.maxTwoHPCands = cms.EDFilter("PATCandViewCountFilter", src = cms.InputTag("hpCands"), minNumber = cms.uint32(0), maxNumber = cms.uint32(2))
process.goodTracks = cms.EDFilter("TrackSelector",
            src = cms.InputTag("generalTracks"),
            cut = kaonSelection,
            )
process.twoGoodTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("goodTracks"), minNumber = cms.uint32(2))
process.goodKaons = cms.EDProducer("ChargedCandidateProducer",
            src = cms.InputTag("goodTracks"),
            particleType = cms.string('pi+')
            )
process.goodDiKaons = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = diKaSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodKaons@+ goodKaons@-')
            )
process.oneGoodDiKa = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDiKaons"), minNumber = cms.uint32(1))
process.diKaEvtSel = cms.Sequence(process.twoTracks * process.hpTracks * process.hpCands * process.maxTwoHPCands * process.goodTracks * process.twoGoodTracks * process.goodKaons * process.goodDiKaons * process.oneGoodDiKa)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # UPC ZB triggers
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.primaryVertexRecoveryForUPC_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hiClusterCompatibility)
process.primaryVertexFilterRecoveryForUPC = process.primaryVertexFilter.clone(src = "primaryVertexRecoveryForUPC")

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.colEvtSel *
    process.diKaEvtSel *
    process.primaryVertexRecoveryForUPC
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.diKa_rereco_step = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh10_seq * process.hfNegFilterNTh10_seq * process.diKa * process.oneDiKa * process.cent_seq)

## Adding the VertexComposite tree ################################################################################################

event_filter = cms.untracked.vstring(
        "Flag_colEvtSel",
        "Flag_clusterCompatibilityFilter",
        "Flag_primaryVertexFilter",
        "Flag_primaryVertexFilterRecoveryForUPC",
        "Flag_hfPosFilterNTh7",
        "Flag_hfPosFilterNTh7p3",
        "Flag_hfPosFilterNTh8",
        "Flag_hfPosFilterNTh10",
        "Flag_hfNegFilterNTh7",
        "Flag_hfNegFilterNTh7p6",
        "Flag_hfNegFilterNTh8",
        "Flag_hfNegFilterNTh10",
    )

trig_info = cms.untracked.VPSet([
    # UPC ZB triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*'), filter = cms.string('hltSinglePixelTrackLowPtForUPC'), minN = cms.int32(1)),
  ])

# Add trigger objects
from HLTrigger.Configuration.HLT_HIon_cff import fragment
for module in ["hltTriggerType", "HLTSiStripClusterChargeCutNone", "hltPixelTracksCleanerBySharedHits", "hltESPPixelCPEGeneric", "hltESPTTRHBuilderPixelOnly", "hltOnlineBeamSpotESProducer", "hltESPStripCPEfromTrackAngle", "hltESPPixelCPEFastHIon", "hltESPMeasurementTracker", "hltOnlineBeamSpot", "hltPixelActivityFilterMaxClusters4e2", "hltSiStripRawToDigi", "hltSiStripZeroSuppression", "hltSiPixelClustersGPUPPOnAA", "hltOnlineBeamSpotToGPU", "hltSiPixelDigiErrorsSoAPPOnAA", "hltSiPixelDigisLegacyPPOnAA", "hltSiPixelDigisSoAPPOnAA", "hltSiPixelDigisFromSoAPPOnAA",                         "hltSiPixelDigisPPOnAA", "hltSiPixelClustersLegacyPPOnAA", "hltSiPixelClustersGPUPPOnAA", "hltSiPixelClustersFromSoAPPOnAA", "hltSiPixelClustersFromSoAPPOnAA", "hltSiPixelClustersPPOnAA",
        "hltSiPixelClustersCachePPOnAA", "hltSiPixelRecHitsFromLegacyPPOnAA", "hltSiPixelRecHitsGPUPPOnAA", "hltSiPixelRecHitsFromGPUPPOnAA", "hltSiPixelRecHitsPPOnAA", "hltSiPixelRecHitsSoAFromGPUPPOnAA", "hltSiPixelRecHitsSoAPPOnAA", "hltSiStripExcludedFEDListProducer",                         "hltHITrackingSiStripRawToClustersFacilityZeroSuppression", "hltMeasurementTrackerEventPPOnAA", "hltPixelLayerTripletsForUPCPPOnAA", "hltPixelTracksForUPCFilterLowPtPPOnAA", "hltPixelTracksForUPCFitterPPOnAA",
        "hltPixelTracksTrackingRegionsLowPtForUPCPPOnAA", "hltPixelClusterCheckForUPCPPOnAA", "hltPixelTracksHitDoubletsLowPtForUPCPPOnAA", "hltPixelTracksHitTripletsLowPtForUPCPPOnAA", "hltPixelTracksLowPtForUPCPPOnAA", "hltPixelCandsLowPtForUPCPPOnAA", "hltSinglePixelTrackLowPtForUPC",         "hltTriggerSummaryAOD", "hltBoolEnd"]:
    setattr(process, module, getattr(fragment, module).clone())

process.hltSiStripRawToDigi.ProductLabel = "rawDataRepacker"
process.hltSiPixelDigisLegacyPPOnAA.InputLabel = "rawDataRepacker"
process.hltSiStripExcludedFEDListProducer.ProductLabel = "rawDataRepacker"
process.hltSinglePixelTrackLowPtForUPC = cms.EDFilter("HLTPixelTrackFilter2", saveTags = cms.bool(True), pixelTracks = cms.InputTag( "hltPixelCandsLowPtForUPCPPOnAA" ), minPixelTracks = cms.uint32( 1 ) )

process.HLTBeginSequence = cms.Sequence( process.hltTriggerType + process.hltOnlineBeamSpot )
process.HLTDoSiStripZeroSuppression = cms.Sequence( process.hltSiStripRawToDigi + process.hltSiStripZeroSuppression )
process.HLTDoLocalPixelPPOnAATask = cms.ConditionalTask( process.hltOnlineBeamSpotToGPU , process.hltSiPixelDigiErrorsSoAPPOnAA , process.hltSiPixelDigisLegacyPPOnAA , process.hltSiPixelDigisSoAPPOnAA , process.hltSiPixelDigisFromSoAPPOnAA , process.hltSiPixelDigisPPOnAA , process.               hltSiPixelClustersLegacyPPOnAA , process.hltSiPixelClustersGPUPPOnAA , process.hltSiPixelClustersFromSoAPPOnAA , process.hltSiPixelClustersPPOnAA , process.hltSiPixelClustersCachePPOnAA , process.hltSiPixelRecHitsFromLegacyPPOnAA , process.hltSiPixelRecHitsGPUPPOnAA , process.                    hltSiPixelRecHitsFromGPUPPOnAA , process.hltSiPixelRecHitsPPOnAA , process.hltSiPixelRecHitsSoAFromGPUPPOnAA , process.hltSiPixelRecHitsSoAPPOnAA )
process.HLTDoLocalPixelPPOnAASequence = cms.Sequence( process.HLTDoLocalPixelPPOnAATask )
process.HLTDoLocalStripSequencePPOnAA = cms.Sequence( process.hltSiStripExcludedFEDListProducer + process.hltHITrackingSiStripRawToClustersFacilityZeroSuppression + process.hltMeasurementTrackerEventPPOnAA )
process.HLTRecopixelvertexingSequencePPOnAAForUPCLowPt = cms.Sequence( process.hltPixelLayerTripletsForUPCPPOnAA + process.hltPixelTracksForUPCFilterLowPtPPOnAA + process.hltPixelTracksForUPCFitterPPOnAA + process.hltPixelTracksTrackingRegionsLowPtForUPCPPOnAA + process.                          hltPixelClusterCheckForUPCPPOnAA + process.hltPixelTracksHitDoubletsLowPtForUPCPPOnAA + process.hltPixelTracksHitTripletsLowPtForUPCPPOnAA + process.hltPixelTracksLowPtForUPCPPOnAA )
process.HLTEndSequence = cms.Sequence( process.hltBoolEnd )

for path in ["HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v", "HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v"]:
    setattr(process, f"hltFilter_{path}", process.hltFilter.clone(HLTPaths = [f"{path}*"]))

process.HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v = cms.Path( process.eventFilter_HM + process.hltFilter_HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v + process.HLTBeginSequence + process.HLTDoSiStripZeroSuppression + process.HLTDoLocalPixelPPOnAASequence +    process.HLTDoLocalStripSequencePPOnAA + process.HLTRecopixelvertexingSequencePPOnAAForUPCLowPt + process.hltPixelCandsLowPtForUPCPPOnAA + process.hltSinglePixelTrackLowPtForUPC + process.HLTEndSequence )
process.HLTriggerFinalPath = cms.Path( process.eventFilter_HM + process.hltTriggerSummaryAOD )


from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.diKaAna = particleAna.clone(
  recoParticles = cms.InputTag("diKa"),
  primaryVertices = cms.InputTag("primaryVertexRecoveryForUPC"),
  selectEvents = cms.string("diKa_rereco_step"),
  eventFilterNames = event_filter,
  addTrgObj = cms.untracked.bool(True),
  triggerInfo = trig_info,
  triggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD::ANASKIM"),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('diKa_ana.root'))
# process.p = cms.EndPath(process.diKaAna * process.generalTracksAna * process.hiConformalPixelTracksAna)
process.p = cms.EndPath(process.diKaAna)

#! Define the process schedule !!!!!!!!!!!!!!!!!!
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.diKa_rereco_step,
    process.HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v,
    process.HLTriggerFinalPath,
    process.p
)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Add the event selection filters ###############################################################################################
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.hiClusterCompatibility)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_primaryVertexFilterRecoveryForUPC = cms.Path(process.eventFilter_HM * process.primaryVertexFilterRecoveryForUPC)
process.Flag_hfPosFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7_seq)
process.Flag_hfPosFilterNTh7p3 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7p3_seq)
process.Flag_hfPosFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh8_seq)
process.Flag_hfPosFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh10_seq)
process.Flag_hfNegFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7_seq)
process.Flag_hfNegFilterNTh7p6 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7p6_seq)
process.Flag_hfNegFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh8_seq)
process.Flag_hfNegFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh10_seq)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_clusterCompatibilityFilter , process.Flag_primaryVertexFilter , process.Flag_primaryVertexFilterRecoveryForUPC , process.Flag_hfPosFilterNTh7 , process.Flag_hfPosFilterNTh7p3 , process.Flag_hfPosFilterNTh8 , process.Flag_hfPosFilterNTh10 , process.Flag_hfNegFilterNTh7 , process.Flag_hfNegFilterNTh7p6 , process.Flag_hfNegFilterNTh8 , process.Flag_hfNegFilterNTh10 ]

#! Adding the process schedule !!!!!!!!!!!!!!!!!!
for P in eventFilterPaths:
    process.schedule.insert(0, P)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
