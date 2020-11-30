import copy
import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA
from HLTrigger.Configuration.common import *
from HLTrigger.Configuration.customizeHLTforPatatrack import *
from HLTrigger.Configuration.customizeHLTforPatatrack import consumeCPULegacyProducts as consumeCPULegacyProductsForPatatrack


# customisation for running the "Patatrack" pixel local reconstruction
def customisePixelLocalReconstructionHIon(process):
    process = customisePixelLocalReconstruction(process)

    if not 'HLTDoLocalPixelSequence' in process.__dict__:
        return process


    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    process.HLTDoLocalPixelSequencePPOnAA = cms.Sequence()
    process.HLTDoLocalPixelSequencePPOnAAForLowPt = cms.Sequence()


    # Modules and EDAliases

    # referenced in HLTDoLocalPixelTask

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
    process.hltSiPixelClustersPPOnAAForLowPt = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltSiPixelClustersPPOnAAForLowPt,
        # alias used to access products from multiple conversion modules
        cuda = cms.EDAlias(
            hltSiPixelDigisClusters = cms.VPSet(
                cms.PSet(type = cms.string("SiPixelClusteredmNewDetSetVector"))
            )
        )
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
    process.hltSiPixelRecHitsPPOnAAForLowPt = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltSiPixelRecHitsPPOnAAForLowPt,
        # converter to legacy format
        cuda = _siPixelRecHitFromSOA.clone(
            pixelRecHitSrc = "hltSiPixelRecHitsCUDA",
            src = "hltSiPixelClustersPPOnAAForLowPt"
        )
    )


    # Tasks and Sequences

    process.HLTDoLocalPixelTaskPP = process.HLTDoLocalPixelTask.clone()

    process.HLTDoLocalPixelTaskPPOnAA = cms.Task(
          process.HLTDoLocalPixelTaskPP,
          process.hltSiPixelClustersPPOnAA,                 # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
          process.hltSiPixelClustersCachePPOnAA,            # legacy module, used by the legacy pixel quadruplet producer
          process.hltSiPixelRecHitsPPOnAA)                  # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA

    process.HLTDoLocalPixelSequencePPOnAA = cms.Sequence(process.HLTDoLocalPixelTaskPPOnAA)

    process.HLTDoLocalPixelTaskPPOnAAForLowPt = cms.Task(
          process.HLTDoLocalPixelTaskPP,
          process.hltSiPixelClustersPPOnAAForLowPt,         # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
          process.hltSiPixelClustersCachePPOnAAForLowPt,    # legacy module, used by the legacy pixel quadruplet producer
          process.hltSiPixelRecHitsPPOnAAForLowPt)          # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA

    process.HLTDoLocalPixelSequencePPOnAAForLowPt = cms.Sequence(process.HLTDoLocalPixelTaskPPOnAAForLowPt)

    process.HLTDoLocalPixelTask = cms.Task(
        process.HLTDoLocalPixelTaskPPOnAA,
        process.HLTDoLocalPixelTaskPPOnAAForLowPt
    )


    # done
    return process


# customisation for running the "Patatrack" pixel track reconstruction
def customisePixelTrackReconstructionHIon(process):
    process = customisePixelTrackReconstruction(process)

    if not 'HLTRecoPixelTracksSequence' in process.__dict__:
        return process


    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    process.HLTRecoPixelTracksPPOnAASequence = cms.Sequence()
    process.HLTRecoPixelTracksSequencePPOnAA = cms.Sequence()
    process.HLTRecopixelvertexingSequencePPOnAA = cms.Sequence()
    process.HLTPixelVertexingPPOnAASequence = cms.Sequence()
    process.HLTPixelVertexingSequencePPOnAA = cms.Sequence()


    # Modules and EDAliases

    # referenced in process.HLTRecoPixelTracksTask

    # cpu only: convert the pixel rechits from legacy to SoA format
    process.hltSiPixelRecHitSoA.src = "hltSiPixelClustersPPOnAA"

    # convert the pixel tracks from SoA to legacy format
    process.hltPixelTracksPPOnAA.pixelRecHitLegacySrc = "hltSiPixelRecHitsPPOnAA"


    # referenced in process.HLTRecopixelvertexingTask

    # convert the pixel vertices from SoA to legacy format
    process.hltPixelVerticesPPOnAA.TrackCollection = "hltPixelTracksPPOnAA"


    # Tasks and Sequences

    process.HLTRecoPixelTracksTask = cms.Task(
          process.hltPixelTracksTrackingRegionsPPOnAA,      # from the original sequence
          process.hltSiPixelRecHitSoA,                      # pixel rechits on cpu, converted to SoA
          process.hltPixelTracksCUDA,                       # pixel ntuplets on gpu, in SoA format
          process.hltPixelTracksSoA,                        # pixel ntuplets on cpu, in SoA format
          process.hltPixelTracksPPOnAA)                     # pixel tracks on cpu, in legacy format


    process.HLTRecoPixelTracksSequencePPOnAA = cms.Sequence(process.HLTRecoPixelTracksTask)
    process.HLTRecoPixelTracksPPOnAASequence = cms.Sequence(process.HLTRecoPixelTracksTask)

    process.HLTRecopixelvertexingTask = cms.Task(
          process.HLTRecoPixelTracksTask,
          process.hltPixelVerticesCUDA,                     # pixel vertices on gpu, in SoA format
          process.hltPixelVerticesSoA,                      # pixel vertices on cpu, in SoA format
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


def consumeCPULegacyProducts(process):
    process = consumeCPULegacyProductsForPatatrack(process)
    process.hltPixelConsumer.eventProducts = cms.untracked.vstring( 'hltPixelTracksPPOnAA', 'hltPixelVerticesPPOnAA')

    # done
    return process
