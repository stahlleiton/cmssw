import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run3_pp_on_PbPb)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('L1Trigger.L1TGlobal.PrescalesVetosFract_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('file:804d2352-a415-4455-b5c9-fff299c5506c.root'),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('125X_dataRun3_relval_v6')

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.centralityBin)

# Add the VertexComposite producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles
pionSelection = cms.string("(pt > 0.0 && abs(eta) < 3.0) && quality(\"highPurity\")")
diPiSelection = cms.string("charge==0")
process.diPi = generalParticles.clone(
    pdgId = cms.int32(113),
    #postSelection = cms.string("userFloat('vertexProb') >= 1E-6"),
    finalSelection = diPiSelection,
    # daughter information
    daughterInfo = cms.VPSet([
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(+1), selection = pionSelection),
        cms.PSet(pdgId = cms.int32(211), charge = cms.int32(-1), selection = pionSelection),
    ]),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
)
process.onlyOneDiPi = cms.EDFilter("PATCandViewCountFilter", src = cms.InputTag("diPi"), minNumber = cms.uint32(1), maxNumber = cms.uint32(1))

# Add event selection
process.twoTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("generalTracks"), minNumber = cms.uint32(2))
process.goodTracks = cms.EDFilter("TrackSelector",
            src = cms.InputTag("generalTracks"),
            cut = pionSelection,
            )
process.twoGoodTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("goodTracks"), minNumber = cms.uint32(2))
process.goodPions = cms.EDProducer("ChargedCandidateProducer",
            src = cms.InputTag("goodTracks"),
            particleType = cms.string('pi+')
            )
process.twoGoodPions = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodPions"), minNumber = cms.uint32(2))
process.goodDiPions = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = diPiSelection,
            checkCharge = cms.bool(True),
            decay = cms.string('goodPions@- goodPions@+')
            )
process.oneGoodDiPi = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDiPions"), minNumber = cms.uint32(1))
process.diPiEvtSel = cms.Sequence(process.twoTracks * process.goodTracks * process.twoGoodTracks * process.goodPions * process.twoGoodPions * process.goodDiPions * process.oneGoodDiPi)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*',
    'HLT_HIZeroBias_v*'
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.clusterCompatibilityFilter_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.diPiEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Luminosity producer
process.lumiInfo = cms.EDProducer('LumiProducerFromBrilcalc',
                                  lumiFile = cms.string("./lumiData.csv"),
                                  throwIfNotFound = cms.bool(False),
                                  doBunchByBunch = cms.bool(False))
process.lumi_seq = cms.Sequence(process.lumiInfo)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq * process.lumi_seq)
process.diPi_rereco_step = cms.Path(process.eventFilter_HM * process.diPi * process.onlyOneDiPi)

# Add the VertexComposite tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna
process.diPiAna = particleAna.clone(
  recoParticles = cms.InputTag("diPi"),
  selectEvents = cms.string("diPi_rereco_step"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_colEvtSel',
      'Flag_hfCoincFilter2Th4',
      'Flag_primaryVertexFilter',
      'Flag_clusterCompatibilityFilter',
      'Flag_hfPosFilterNTh3',
      'Flag_hfNegFilterNTh3',
      'Flag_hfPosFilterNTh4',
      'Flag_hfNegFilterNTh4',
      'Flag_hfPosFilterNTh5',
      'Flag_hfNegFilterNTh5',
      'Flag_hfPosFilterNTh6',
      'Flag_hfNegFilterNTh6',
      'Flag_hfPosFilterNTh7',
      'Flag_hfNegFilterNTh7',
      'Flag_hfPosFilterNTh8',
      'Flag_hfNegFilterNTh8',
      'Flag_hfPosFilterNTh7p3',
      'Flag_hfNegFilterNTh7p6'
  ),
  triggerInfo = cms.untracked.VPSet([
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIZeroBias_v*')),
  ]),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('diPi_ana.root'))
process.p = cms.EndPath(process.diPiAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pcentandep_step,
    process.diPi_rereco_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.clusterCompatibilityFilter)
process.Flag_hfPosFilterNTh3 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh3_seq)
process.Flag_hfPosFilterNTh4 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh4_seq)
process.Flag_hfPosFilterNTh5 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh5_seq)
process.Flag_hfPosFilterNTh6 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh6_seq)
process.Flag_hfPosFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7_seq)
process.Flag_hfPosFilterTh8 = cms.Path(process.eventFilter_HM * process.hfPosFilterTh8_seq)
process.Flag_hfPosFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh8_seq)
process.Flag_hfPosFilterNTh7p3 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7p3_seq)
process.Flag_hfPosFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh10_seq)
process.Flag_hfNegFilterNTh3 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh3_seq)
process.Flag_hfNegFilterNTh4 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh4_seq)
process.Flag_hfNegFilterNTh5 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh5_seq)
process.Flag_hfNegFilterNTh6 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh6_seq)
process.Flag_hfNegFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7_seq)
process.Flag_hfNegFilterTh8 = cms.Path(process.eventFilter_HM * process.hfNegFilterTh8_seq)
process.Flag_hfNegFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh8_seq)
process.Flag_hfNegFilterNTh7p6 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7p6_seq)
process.Flag_hfNegFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh10_seq)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter, process.Flag_clusterCompatibilityFilter, process.Flag_hfPosFilterNTh3, process.Flag_hfNegFilterNTh3,process.Flag_hfPosFilterNTh4, process.Flag_hfNegFilterNTh4, process.                Flag_hfPosFilterNTh5, process.Flag_hfNegFilterNTh5, process.Flag_hfPosFilterNTh6, process.Flag_hfNegFilterNTh6, process.Flag_hfPosFilterNTh7, process.Flag_hfNegFilterNTh7, process.Flag_hfPosFilterTh8, process.Flag_hfPosFilterNTh8, process.Flag_hfNegFilterTh8, process.Flag_hfNegFilterNTh8,        process.Flag_hfPosFilterNTh7p3, process.Flag_hfNegFilterNTh7p6, process.Flag_hfPosFilterNTh10, process.Flag_hfNegFilterNTh10 ]

process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.clusterCompatibilityFilter *
    process.hfPosFilterNTh10_seq *
    process.hfNegFilterNTh10_seq *
    #process.diPiEvtSel *
    process.primaryVertexFilter
)

for P in eventFilterPaths:
    process.schedule.insert(0, P)
