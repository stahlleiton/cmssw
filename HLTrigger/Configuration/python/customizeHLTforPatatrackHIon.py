import copy
import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA
from HLTrigger.Configuration.common import *
from HLTrigger.Configuration.customizeHLTforPatatrack import *


# customisation for running the "Patatrack" pixel local reconstruction
def customisePixelLocalReconstructionHIon(process):

    if not 'HLTDoLocalPixelSequence' in process.__dict__:
        return process


    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    process.HLTDoLocalPixelSequencePPOnAA = cms.Sequence()


    # Event Setup

    process.load("CalibTracker.SiPixelESProducers.siPixelGainCalibrationForHLTGPU_cfi")         # this should be used only on GPUs, will crash otherwise
    process.load("RecoLocalTracker.SiPixelClusterizer.siPixelFedCablingMapGPUWrapper_cfi")      # this should be used only on GPUs, will crash otherwise
    process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEFastESProducer_cfi")
    process.PixelCPEFastESProducer.DoLorentz = True


    # Modules and EDAliases

    # referenced in HLTDoLocalPixelTask

    # transfer the beamspot to the gpu
    from RecoVertex.BeamSpotProducer.offlineBeamSpotCUDA_cfi import offlineBeamSpotCUDA as _offlineBeamSpotCUDA
    process.hltOnlineBeamSpotCUDA = _offlineBeamSpotCUDA.clone(
        src = "hltOnlineBeamSpot"
    )

    # reconstruct the pixel digis and clusters on the gpu
    from RecoLocalTracker.SiPixelClusterizer.siPixelRawToClusterCUDA_cfi import siPixelRawToClusterCUDA as _siPixelRawToClusterCUDA
    process.hltSiPixelClustersCUDA = _siPixelRawToClusterCUDA.clone()

    # copy the pixel digis errors to the host
    from EventFilter.SiPixelRawToDigi.siPixelDigiErrorsSoAFromCUDA_cfi import siPixelDigiErrorsSoAFromCUDA as _siPixelDigiErrorsSoAFromCUDA
    process.hltSiPixelDigiErrorsSoA = _siPixelDigiErrorsSoAFromCUDA.clone(
        src = "hltSiPixelClustersCUDA"
    )

    # convert the pixel digis errors to the legacy format
    from EventFilter.SiPixelRawToDigi.siPixelDigiErrorsFromSoA_cfi import siPixelDigiErrorsFromSoA as _siPixelDigiErrorsFromSoA
    process.hltSiPixelDigiErrors = _siPixelDigiErrorsFromSoA.clone(
        digiErrorSoASrc = "hltSiPixelDigiErrorsSoA",
        UsePhase1 = True
    )

    # copy the pixel digis (except errors) and clusters to the host
    from EventFilter.SiPixelRawToDigi.siPixelDigisSoAFromCUDA_cfi import siPixelDigisSoAFromCUDA as _siPixelDigisSoAFromCUDA
    process.hltSiPixelDigisSoA = _siPixelDigisSoAFromCUDA.clone(
        src = "hltSiPixelClustersCUDA"
    )

    # convert the pixel digis (except errors) and clusters to the legacy format
    from RecoLocalTracker.SiPixelClusterizer.siPixelDigisClustersFromSoA_cfi import siPixelDigisClustersFromSoA as _siPixelDigisClustersFromSoA
    process.hltSiPixelDigisClusters = _siPixelDigisClustersFromSoA.clone(
        src = "hltSiPixelDigisSoA"
    )

    # SwitchProducer wrapping the legacy pixel digis producer or an alias combining the pixel digis information converted from SoA
    process.hltSiPixelDigis = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltSiPixelDigis,
        # alias used to access products from multiple conversion modules
        cuda = cms.EDAlias(
            hltSiPixelDigisClusters = cms.VPSet(
                cms.PSet(type = cms.string("PixelDigiedmDetSetVector"))
            ),
            hltSiPixelDigiErrors = cms.VPSet(
                cms.PSet(type = cms.string("DetIdedmEDCollection")),
                cms.PSet(type = cms.string("SiPixelRawDataErroredmDetSetVector")),
                cms.PSet(type = cms.string("PixelFEDChanneledmNewDetSetVector"))
            )
        )
    )

    # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
    process.hltSiPixelClustersPPOnAA = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltSiPixelClustersPPOnAA,
        # alias used to access products from multiple conversion modules
        cuda = cms.EDAlias(
            hltSiPixelDigisClusters = cms.VPSet(
                cms.PSet(type = cms.string("SiPixelClusteredmNewDetSetVector"))
            )
        )
    )

    # reconstruct the pixel rechits on the gpu
    from RecoLocalTracker.SiPixelRecHits.siPixelRecHitCUDA_cfi import siPixelRecHitCUDA as _siPixelRecHitCUDA
    process.hltSiPixelRecHitsCUDA = _siPixelRecHitCUDA.clone(
        src = "hltSiPixelClustersCUDA",
        beamSpot = "hltOnlineBeamSpotCUDA"
    )

    # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA
    from RecoLocalTracker.SiPixelRecHits.siPixelRecHitFromSOA_cfi import siPixelRecHitFromSOA as _siPixelRecHitFromSOA
    process.hltSiPixelRecHitsPPOnAA = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltSiPixelRecHitsPPOnAA,
        # converter to legacy format
        cuda = _siPixelRecHitFromSOA.clone(
            pixelRecHitSrc = "hltSiPixelRecHitsCUDA",
            src = "hltSiPixelClustersPPOnAA"
        )
    )
    process.hltSiPixelRecHitsPPOnAAForLowPt = process.hltSiPixelRecHitsPPOnAA.clone()


    # Tasks and Sequences

    process.HLTDoLocalPixelTask = cms.Task(
          process.hltOnlineBeamSpotCUDA,                    # transfer the beamspot to the gpu
          process.hltSiPixelClustersCUDA,                   # reconstruct the pixel digis and clusters on the gpu
          process.hltSiPixelRecHitsCUDA,                    # reconstruct the pixel rechits on the gpu
          process.hltSiPixelDigisSoA,                       # copy the pixel digis (except errors) and clusters to the host
          process.hltSiPixelDigisClusters,                  # convert the pixel digis (except errors) and clusters to the legacy format
          process.hltSiPixelDigiErrorsSoA,                  # copy the pixel digis errors to the host
          process.hltSiPixelDigiErrors,                     # convert the pixel digis errors to the legacy format
          process.hltSiPixelDigis,                          # SwitchProducer wrapping the legacy pixel digis producer or an alias combining the pixel digis information converted from SoA
    )

    process.HLTDoLocalPixelTaskPPOnAA = cms.Task(
          process.HLTDoLocalPixelTask,
          process.hltSiPixelClustersPPOnAA,                 # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
          process.hltSiPixelClustersCachePPOnAA,            # legacy module, used by the legacy pixel quadruplet producer
          process.hltSiPixelRecHitsPPOnAA)                  # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA

    process.HLTDoLocalPixelSequencePPOnAA = cms.Sequence(process.HLTDoLocalPixelTaskPPOnAA)

    process.HLTDoLocalPixelTaskPPOnAAForLowPt = cms.Task(
          process.HLTDoLocalPixelTask,
          process.hltSiPixelClustersPPOnAAForLowPt,         # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
          process.hltSiPixelClustersCachePPOnAAForLowPt,    # legacy module, used by the legacy pixel quadruplet producer
          process.hltSiPixelRecHitsPPOnAAForLowPt)          # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA

    process.HLTDoLocalPixelSequencePPOnAAForLowPt = cms.Sequence(process.HLTDoLocalPixelTaskPPOnAAForLowPt)


    # done
    return process


# customisation for running the "Patatrack" pixel track reconstruction
def customisePixelTrackReconstructionHIon(process):

    if not 'HLTRecoPixelTracksSequence' in process.__dict__:
        return process


    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    process.HLTRecoPixelTracksPPOnAASequence = cms.Sequence()
    process.HLTRecopixelvertexingSequencePPOnAA = cms.Sequence()
    process.HLTDoLocalPixelSequencePPOnAAForLowPt = cms.Sequence()


    # Modules and EDAliases

    # referenced in process.HLTRecoPixelTracksTask

    # cpu only: convert the pixel rechits from legacy to SoA format
    from RecoLocalTracker.SiPixelRecHits.siPixelRecHitHostSoA_cfi import siPixelRecHitHostSoA as _siPixelRecHitHostSoA
    process.hltSiPixelRecHitSoAPPOnAA = _siPixelRecHitHostSoA.clone(
        src = "hltSiPixelClustersPPOnAA",
        beamSpot = "hltOnlineBeamSpot",
        convertToLegacy = True
    )

    # build pixel ntuplets and pixel tracks in SoA format on gpu
    from RecoPixelVertexing.PixelTriplets.caHitNtupletCUDA_cfi import caHitNtupletCUDA as _caHitNtupletCUDA
    process.hltPixelTracksCUDA = _caHitNtupletCUDA.clone(
        idealConditions = False,
        pixelRecHitSrc = "hltSiPixelRecHitsCUDA",
        onGPU = True
    )

    # SwitchProducer providing the pixel tracks in SoA format on cpu
    process.hltPixelTracksSoAPPOnAA = SwitchProducerCUDA(
        # build pixel ntuplets and pixel tracks in SoA format on cpu
        cpu = _caHitNtupletCUDA.clone(
            idealConditions = False,
            pixelRecHitSrc = "hltSiPixelRecHitSoAPPOnAA",
            onGPU = False
        ),
        # transfer the pixel tracks in SoA format to the host
        cuda = cms.EDProducer("PixelTrackSoAFromCUDA",
            src = cms.InputTag("hltPixelTracksCUDA")
        )
    )

    # convert the pixel tracks from SoA to legacy format
    from RecoPixelVertexing.PixelTrackFitting.pixelTrackProducerFromSoA_cfi import pixelTrackProducerFromSoA as _pixelTrackProducerFromSoA
    process.hltPixelTracksPPOnAA = _pixelTrackProducerFromSoA.clone(
        beamSpot = "hltOnlineBeamSpot",
        pixelRecHitLegacySrc = "hltSiPixelRecHitsPPOnAA",
        trackSrc = "hltPixelTracksSoAPPOnAA"
    )


    # referenced in process.HLTRecopixelvertexingTask

    # build pixel vertices in SoA format on gpu
    from RecoPixelVertexing.PixelVertexFinding.pixelVertexCUDA_cfi import pixelVertexCUDA as _pixelVertexCUDA
    process.hltPixelVerticesCUDA = _pixelVertexCUDA.clone(
        pixelTrackSrc = "hltPixelTracksCUDA",
        onGPU = True
    )

    # build or transfer pixel vertices in SoA format on cpu
    process.hltPixelVerticesSoAPPOnAA = SwitchProducerCUDA(
        # build pixel vertices in SoA format on cpu
        cpu = _pixelVertexCUDA.clone(
            pixelTrackSrc = "hltPixelTracksSoAPPOnAA",
            onGPU = False
        ),
        # transfer the pixel vertices in SoA format to cpu
        cuda = cms.EDProducer("PixelVertexSoAFromCUDA",
            src = cms.InputTag("hltPixelVerticesCUDA")
        )
    )

    # convert the pixel vertices from SoA to legacy format
    from RecoPixelVertexing.PixelVertexFinding.pixelVertexFromSoA_cfi import pixelVertexFromSoA as _pixelVertexFromSoA
    process.hltPixelVerticesPPOnAA = _pixelVertexFromSoA.clone(
        src = "hltPixelVerticesSoAPPOnAA",
        TrackCollection = "hltPixelTracksPPOnAA",
        beamSpot = "hltOnlineBeamSpot"
    )


    # Tasks and Sequences

    process.HLTRecoPixelTracksTask = cms.Task(
          process.hltPixelTracksTrackingRegionsPPOnAA,      # from the original sequence
          process.hltSiPixelRecHitSoAPPOnAA,                # pixel rechits on cpu, converted to SoA
          process.hltPixelTracksCUDA,                       # pixel ntuplets on gpu, in SoA format
          process.hltPixelTracksSoAPPOnAA,                  # pixel ntuplets on cpu, in SoA format
          process.hltPixelTracksPPOnAA)                     # pixel tracks on cpu, in legacy format


    process.HLTRecoPixelTracksSequencePPOnAA = cms.Sequence(process.HLTRecoPixelTracksTask)
    process.HLTRecoPixelTracksPPOnAASequence = cms.Sequence(process.HLTRecoPixelTracksSequencePPOnAA)

    process.HLTRecopixelvertexingTask = cms.Task(
          process.HLTRecoPixelTracksTask,
          process.hltPixelVerticesCUDA,                     # pixel vertices on gpu, in SoA format
          process.hltPixelVerticesSoAPPOnAA,                # pixel vertices on cpu, in SoA format
          process.hltPixelVerticesPPOnAA,                   # pixel vertices on cpu, in legacy format
          process.hltTrimmedPixelVerticesPPOnAA)            # from the original sequence

    process.HLTRecopixelvertexingSequencePPOnAA = cms.Sequence(
          process.hltPixelTracksFitter +                    # not used here, kept for compatibility with legacy sequences
          process.hltPixelTracksFilter,                     # not used here, kept for compatibility with legacy sequences
          process.HLTRecopixelvertexingTask)
    process.HLTPixelVertexingPPOnAASequence = cms.Sequence(process.HLTRecopixelvertexingSequencePPOnAA)
    process.HLTPixelVertexingSequencePPOnAA = cms.Sequence(process.HLTRecopixelvertexingSequencePPOnAA)


    # done
    return process


# customisation for offloading the HCAL local reconstruction via CUDA if a supported gpu is present
def customiseHcalLocalReconstructionHIon(process):
    process = customiseHcalLocalReconstruction(process)
    process.HLTDoLocalHcalWithTowerSequence = cms.Sequence(process.hltTowerMakerForAll, process.HLTDoLocalHcalTask)


    # done
    return process


# customisation for running the Patatrack reconstruction, with automatic offload via CUDA when a supported gpu is available
def customizeHLTforPatatrackHIon(process):
    process = customiseCommon(process)
    process = customisePixelLocalReconstructionHIon(process)
    process = customisePixelTrackReconstructionHIon(process)
    process = customiseEcalLocalReconstruction(process)
    process = customiseHcalLocalReconstructionHIon(process)
    return process
