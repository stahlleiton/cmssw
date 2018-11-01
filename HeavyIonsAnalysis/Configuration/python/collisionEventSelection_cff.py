import FWCore.ParameterSet.Config as cms

# Coincidence of HF towers above threshold
from HeavyIonsAnalysis.Configuration.hfCoincFilter_cff import *

# Selection of at least a two-track fitted vertex
primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True), # otherwise it won't filter the events
)

beamScrapingFilter = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

# Cluster-shape filter re-run offline
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from HLTrigger.special.hltPixelClusterShapeFilter_cfi import *
hltPixelClusterShapeFilter.inputTag = "siPixelRecHits"

# Cluster-shape filter re-run offline from ClusterCompatibility object
from HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi import *

collisionEventSelection = cms.Sequence(
    hfCoincFilter3 *
    primaryVertexFilter *
    siPixelRecHits *
    hltPixelClusterShapeFilter)

collisionEventSelectionAOD = cms.Sequence(
    hfCoincFilter3 *
    primaryVertexFilter *
    clusterCompatibilityFilter)
