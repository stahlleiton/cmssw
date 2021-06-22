import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.customizeHLTforPatatrack import *
from HLTrigger.Configuration.customizeHLTforPatatrack import consumeCPULegacyProducts as consumeCPULegacyProductsForPatatrack


# customisation for running the "Patatrack" pixel local reconstruction
def customisePixelLocalReconstructionHIon(process):
    process.HLTDoLocalPixelSequence = cms.Sequence()
    from HLTrigger.Configuration.HLT_FULL_cff import fragment
    process.hltSiPixelDigis = fragment.hltSiPixelDigis.clone()
    process.hltSiPixelClusters = fragment.hltSiPixelClusters.clone()
    process.hltSiPixelClustersCache = fragment.hltSiPixelClustersCache.clone()
    process.hltSiPixelRecHits = fragment.hltSiPixelRecHits.clone()
    process = customisePixelLocalReconstruction(process)

    if not 'HLTDoLocalPixelSequence' in process.__dict__:
        return process


    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    process.HLTDoLocalPixelSequencePPOnAA = cms.Sequence()
    process.HLTDoLocalPixelSequencePPOnAAForLowPt = cms.Sequence()
    process.HLTHIDoLocalPixelSequence = cms.Sequence()


    # Modules and EDAliases

    # referenced in HLTDoLocalPixelTask

    # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
    process.hltSiPixelClustersPPOnAA = process.hltSiPixelClusters.clone(cpu = process.hltSiPixelClustersPPOnAA)
    process.hltSiPixelClustersPPOnAAForLowPt = process.hltSiPixelClusters.clone(cpu = process.hltSiPixelClustersPPOnAAForLowPt)
    process.hltHISiPixelDigis = SwitchProducerCUDA(
        # legacy producer
        cpu = process.hltHISiPixelDigis,
        # alias used to access products from multiple conversion modules
        cuda = cms.EDAlias(
            hltSiPixelDigisClusters = cms.VPSet(
                cms.PSet(type = cms.string("PixelDigiedmDetSetVector"))
            )
        )
    )
    process.hltHISiPixelClusters = process.hltSiPixelClusters.clone(cpu = process.hltHISiPixelClusters)

    # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA
    process.hltSiPixelRecHitsPPOnAA = process.hltSiPixelRecHits.clone(cpu = process.hltSiPixelRecHitsPPOnAA)
    process.hltSiPixelRecHitsPPOnAA.cuda.src = "hltSiPixelClustersPPOnAA"
    process.hltSiPixelRecHitsPPOnAAForLowPt = process.hltSiPixelRecHits.clone(cpu = process.hltSiPixelRecHitsPPOnAAForLowPt)
    process.hltSiPixelRecHitsPPOnAAForLowPt.cuda.src = "hltSiPixelClustersPPOnAAForLowPt"
    process.hltHISiPixelRecHits = process.hltSiPixelRecHits.clone(cpu = process.hltHISiPixelRecHits)
    process.hltHISiPixelRecHits.cuda.src = "hltHISiPixelClusters"


    # Tasks and Sequences

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

    process.HLTHIDoLocalPixelSequence = cms.Sequence(
          process.hltHISiPixelDigis *                       # SwitchProducer wrapping the legacy pixel digis producer or an alias combining the pixel digis information converted from SoA
          process.hltHISiPixelClusters *                    # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
          process.hltHISiPixelClustersCache *               # legacy module, used by the legacy pixel quadruplet producer
          process.hltHISiPixelRecHits,			    # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA
          process.HLTDoLocalPixelTask)


    # done
    return process


# customisation for running the "Patatrack" pixel track reconstruction
def customisePixelTrackReconstructionHIon(process):
    process.HLTRecoPixelTracksSequence = cms.Sequence()
    from HLTrigger.Configuration.HLT_FULL_cff import fragment
    process.hltPixelTracksTrackingRegions = fragment.hltPixelTracksTrackingRegions.clone()
    process.hltTrimmedPixelVertices = fragment.hltTrimmedPixelVertices.clone()
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
    process.hltPixelTracksPPOnAA = process.hltPixelTracks.clone(pixelRecHitLegacySrc = "hltSiPixelRecHitsPPOnAA")


    # referenced in process.HLTRecopixelvertexingTask

    # convert the pixel vertices from SoA to legacy format
    process.hltPixelVerticesPPOnAA = process.hltPixelVertices.clone(TrackCollection = "hltPixelTracksPPOnAA")


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

    if not 'HLTDoLocalHcalSequence' in process.__dict__:
        return process

    # add HLT tower maker
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


# customisation for running the Patatrack triplets reconstruction, with automatic offload via CUDA when a supported gpu is available
def customizeHLTforPatatrackTripletsHIon(process):
    process = customizeHLTforPatatrackHIon(process)
    process = enablePatatrackPixelTriplets(process)
    return process


def consumeCPULegacyProducts(process):
    process = consumeCPULegacyProductsForPatatrack(process)
    process.hltPixelConsumer.eventProducts = cms.untracked.vstring( 'hltPixelTracksPPOnAA', 'hltPixelVerticesPPOnAA')

    # done
    return process
