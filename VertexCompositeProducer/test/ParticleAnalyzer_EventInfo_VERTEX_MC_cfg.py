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
process.MessageLogger.cerr.threshold = 'ERROR'

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/phys_heavyions/jiazhao/STARlight/2023Run3/Reco/STARlight_CohPhiToKK_Reco_132X_240125_044529/STARlight/CohPhiToKK/240125_034539/0000/step3_STARlight_Reco_1.root'
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_mcRun3_2023_realistic_HI_v9')

# Set ZDC information
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    authenticationMethod = cms.untracked.uint32(1)
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("emap_2023_newZDC_v3.txt")))
)

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db"),
        label = cms.untracked.string("HFtowers")
        )
    ]
)
process.cent_seq = cms.Sequence(process.centralityBin)

# Add primary vertex
from RecoVertex.PrimaryVertexProducer.TkClusParameters_cff import DA_vectParameters
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVertices
process.offlinePrimaryVerticesPP = offlinePrimaryVertices.clone(
        verbose = cms.untracked.bool(False),
        TrackLabel = cms.InputTag("generalTracks"),
        beamSpotLabel = cms.InputTag("offlineBeamSpot"),

        TkFilterParameters = cms.PSet(
            algorithm=cms.string('filter'),
            maxNormalizedChi2 = cms.double(10.0),
            minPixelLayersWithHits=cms.int32(2),
            minSiliconLayersWithHits = cms.int32(5),
            maxD0Significance = cms.double(4.0), 
            maxD0Error = cms.double(1.0), 
            maxDzError = cms.double(1.0), 
            minPt = cms.double(0.0),
            maxEta = cms.double(2.4),
            trackQuality = cms.string("any")
        ),

        TkClusParameters = cms.PSet(
                algorithm   = cms.string("DA_vect"),
                TkDAClusParameters = cms.PSet(
                coolingFactor = cms.double(0.6),  # moderate annealing speed
                zrange = cms.double(4.),          # consider only clusters within 4 sigma*sqrt(T) of a track
                delta_highT = cms.double(1.e-2),  # convergence requirement at high T
                delta_lowT = cms.double(1.e-3),   # convergence requirement at low T
                convergence_mode = cms.int32(0),  # 0 = two steps, 1 = dynamic with sqrt(T)
                Tmin = cms.double(2.0),           # end of vertex splitting
                Tpurge = cms.double(2.0),         # cleaning 
                Tstop = cms.double(0.5),          # end of annealing 
                vertexSize = cms.double(0.006),   # added in quadrature to track-z resolutions        
                d0CutOff = cms.double(3.),        # downweight high IP tracks 
                dzCutOff = cms.double(3.),        # outlier rejection after freeze-out (T<Tmin)       
                zmerge = cms.double(1e-2),        # merge intermediat clusters separated by less than zmerge
                uniquetrkweight = cms.double(0.8),# require at least two tracks with this weight at T=Tpurge
                uniquetrkminp = cms.double(0.0),  # minimal a priori track weight for counting unique tracks
                runInBlocks = cms.bool(False),    # activate the DA running in blocks of z sorted tracks
                block_size = cms.uint32(10000),   # block size in tracks
                overlap_frac = cms.double(0.0)    # overlap between consecutive blocks (blocks_size*overlap_frac)
            )
        ),


        vertexCollections = cms.VPSet(
            [cms.PSet(label=cms.string(""),
               algorithm=cms.string("AdaptiveVertexFitter"),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(0.0),
               useBeamConstraint = cms.bool(False),
               maxDistanceToBeam = cms.double(1.0)
               ),
            cms.PSet(label=cms.string("WithBS"),
               algorithm = cms.string('AdaptiveVertexFitter'),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(2.0),
               useBeamConstraint = cms.bool(True),
               maxDistanceToBeam = cms.double(1.0),
               )
            ]
        ),

        isRecoveryIteration = cms.bool(False),
        recoveryVtxCollection = cms.InputTag("")
)
process.offlinePrimaryVerticesHI = offlinePrimaryVertices.clone(
        TkFilterParameters = dict(
            algorithm="filterWithThreshold",
            maxD0Significance = 2.0,
            maxD0Error = 10.0,
            maxDzError = 10.0,
            minPixelLayersWithHits = 3,
            minSiliconLayersWithHits = 5,
            maxEta = 2.4,
            maxNormalizedChi2 = 10.0,
            minPt = 0.7,
            trackQuality = "highPurity",
            numTracksThreshold = cms.int32(10),
            maxNumTracksThreshold = cms.int32(1000),
            minPtTight = cms.double(1.0)
        ),
        TkClusParameters = cms.PSet(
            algorithm = cms.string("gap"),
            TkGapClusParameters = cms.PSet(
                zSeparation = cms.double(1.0)
            )
        ),
        vertexCollections = {
            0: dict(chi2cutoff = 4.0, minNdof = -1.1),
            1: dict(chi2cutoff = 4.0, minNdof = -2.0),
        }
)
process.load('VertexCompositeAnalysis.VertexCompositeProducer.primaryVertexRecoveryForUPC_cfi')
process.nTracks = cms.EDFilter("TrackCountFilter", src = cms.InputTag("generalTracks"), minNumber = cms.uint32(2))
process.primaryVertexPath = cms.Path(process.nTracks + process.offlinePrimaryVerticesPP + process.offlinePrimaryVerticesHI + process.primaryVertexRecoveryForUPC)

# Add the Particle tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc

process.eventAna = particleAna_mc.clone(
  recoParticles = cms.InputTag(""),
  primaryVertices = cms.InputTag("primaryVertexRecoveryForUPC"),
  selectEvents = cms.string("primaryVertexPath"),
  eventFilterNames = cms.untracked.vstring(
      'Flag_primaryVertexFilter',
      'Flag_primaryVertexFilterPP',
      'Flag_primaryVertexFilterHI',
      'Flag_primaryVertexFilterRecoveryForUPC',
      'Flag_hiClusterCompatibility',
  ),
  triggerInfo = cms.untracked.VPSet([
    # UPC muon triggers
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_BptxAND_MaxPixelCluster1000_v*')), #0
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2AND_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2AND_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2OR_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2OR_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2OR_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2OR_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_OR_SingleMuCosmic_EMTF_BptxAND_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_OR_SingleMuCosmic_EMTF_NotMBHF2AND_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_OR_SingleMuCosmic_EMTF_NotMBHF2AND_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_OR_SingleMuCosmic_EMTF_NotMBHF2OR_MaxPixelCluster1000_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_OR_SingleMuCosmic_EMTF_NotMBHF2OR_v*')),
    # UPC egamma triggers
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleEG2_BptxAND_SinglePixelTrack_MaxPixelTrack_v*')), #15
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleEG2_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleEG2_NotMBHF2AND_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleEG5_BptxAND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleEG5_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_DoubleEG5_NotMBHF2AND_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG2_NotMBHF2AND_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_BptxAND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2AND_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2OR_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG5_BptxAND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG5_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG5_NotMBHF2AND_v*')),
    # UPC zero bias triggers
    cms.PSet(path = cms.string('HLT_HIZeroBias_v*')), #30
    cms.PSet(path = cms.string('HLT_HIZeroBias_HighRate_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # UPC ZDC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrack_MaxPixelTrack_v*')), #35
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*')),
    cms.PSet(path = cms.string('HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*')),
    # Empty BX triggers
    cms.PSet(path = cms.string('HLT_HIL1NotBptxOR_v*')), #38
    cms.PSet(path = cms.string('HLT_HIL1UnpairedBunchBptxMinus_v*')), #39
    cms.PSet(path = cms.string('HLT_HIL1UnpairedBunchBptxPlus_v*')), #40
  ]),
  dataset = cms.untracked.string("Dataset_HIForward*"),
)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('event_ana.root'))
process.p = cms.EndPath(process.eventAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.primaryVertexPath,
    process.p
)

# Add the event selection filters
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.primaryVertexFilterPP = process.primaryVertexFilter.clone(src = cms.InputTag("offlinePrimaryVerticesPP"))
process.primaryVertexFilterHI = process.primaryVertexFilter.clone(src = cms.InputTag("offlinePrimaryVerticesHI"))
process.primaryVertexFilterRecoveryForUPC = process.primaryVertexFilter.clone(src = cms.InputTag("primaryVertexRecoveryForUPC"))
process.Flag_primaryVertexFilter = cms.Path(process.nTracks * process.primaryVertexFilter)
process.Flag_primaryVertexFilterPP = cms.Path(process.nTracks * process.offlinePrimaryVerticesPP * process.primaryVertexFilterPP)
process.Flag_primaryVertexFilterHI = cms.Path(process.nTracks * process.offlinePrimaryVerticesHI * process.primaryVertexFilterHI)
process.Flag_primaryVertexFilterRecoveryForUPC = cms.Path(process.nTracks * process.primaryVertexRecoveryForUPC * process.primaryVertexFilterRecoveryForUPC)
process.Flag_hiClusterCompatibility = cms.Path(process.hiClusterCompatibility)

eventFilterPaths = [ process.Flag_primaryVertexFilter , process.Flag_primaryVertexFilterPP , process.Flag_primaryVertexFilterHI , process.Flag_primaryVertexFilterRecoveryForUPC , process.Flag_hiClusterCompatibility ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)
