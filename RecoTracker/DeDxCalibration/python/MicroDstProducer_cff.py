import FWCore.ParameterSet.Config as cms

from RecoTracker.DeDx.dedxHitCalibrator_cfi import dedxHitCalibrator as _dedxHitCalibrator
from SimGeneral.MixingModule.SiStripSimParameters_cfi import SiStripSimBlock as _SiStripSimBlock
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import siPixelClusters as _siPixelClusters

# Micro Dst
microDstProducer = cms.EDAnalyzer("MicroDstProducer",
  MeVPerElectron = cms.double(1000*_SiStripSimBlock.GevPerElectron.value()),
  VCaltoElectronGain = _siPixelClusters.VCaltoElectronGain,
  VCaltoElectronGain_L1 = _siPixelClusters.VCaltoElectronGain_L1,
  VCaltoElectronOffset = _siPixelClusters.VCaltoElectronOffset,
  VCaltoElectronOffset_L1 = _siPixelClusters.VCaltoElectronOffset_L1,
  pixelSaturationThr = cms.int32(254),
  trackProducer = cms.InputTag('generalTracks'),
  dedxHitInfo   = cms.InputTag('dedxEstimator'),
  dedxMomentum = cms.InputTag('dedxEstimator:momentumAtHit')
)

# Paths
produceMicroDst   = cms.Path(microDstProducer)
