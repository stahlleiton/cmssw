import FWCore.ParameterSet.Config as cms

unpackedMuons = cms.EDProducer('MuonUnpacker',
  muons = cms.InputTag('slimmedMuons'),
  tracks = cms.InputTag('unpackedTracksAndVertices'),
  primaryVertices = cms.InputTag('unpackedTracksAndVertices'),
  beamSpot = cms.InputTag('offlineBeamSpot'),
  muonSelectors = cms.vstring(
    'AllTrackerMuons',
    'TMOneStationTight'
  ),
  addPropToMuonSt = cms.bool(False),
  useTrack = cms.string('tracker'),
  useState = cms.string('atVertex'),
  useSimpleGeometry = cms.bool(True),
  useStation2 = cms.bool(True),
  fallbackToME1 = cms.bool(True),
  useMB2InOverlap = cms.bool(True),
  cosmicPropagationHypothesis = cms.bool(False),
  propagatorAlong = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorAlong'),
  propagatorAny = cms.ESInputTag('', 'SteppingHelixPropagatorAny'),
  propagatorOpposite = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorOpposite'),
  mightGet = cms.optional.untracked.vstring
)
