import FWCore.ParameterSet.Config as cms

# Tracker local reco
from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *

# Refitter
from RecoTracker.TrackProducer.TrackRefitters_cff import *

# Track refitter
from RecoTracker.TrackProducer.TrackRefitter_cfi import *
refitterForEnergyLoss     = TrackRefitter.clone()
refitterForEnergyLoss.src = 'generalTracks'

refitterForEnergyLoss.MeasurementTrackerEvent = cms.InputTag('MeasurementTrackerEvent')

# Energy Loss
energyLossProducer   = cms.EDProducer("EnergyLossProducer",
  trackProducer      = cms.InputTag('refitterForEnergyLoss'),
  trackProducer2     = cms.InputTag('generalTracks'),
  trajectoryProducer = cms.InputTag('refitterForEnergyLoss'),
  clusterShapeCache  = cms.InputTag('siPixelClusterShapeCache'),
  dedxProducer       = cms.InputTag('dedxPixelHarmonic2'),
  tag = cms.string('PbPb2023'),
  outFile = cms.string('result.dat.gz')
)

# Paths
produceEnergyLoss = cms.Path(refitterForEnergyLoss
                           * energyLossProducer)
