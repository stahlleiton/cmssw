import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.customizeHLTforPatatrack import *
from HLTrigger.Configuration.customizeHLTforPatatrack import consumeCPULegacyProducts as consumeCPULegacyProductsForPatatrack
from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
from HLTrigger.Configuration.HLT_FULL_cff import fragment


# customisation for running the "Patatrack" pixel local reconstruction
def customisePixelLocalReconstructionHIon(process):

    hiSeq = ['HLTDoLocalPixelSequencePPOnAA', 'HLTDoLocalPixelSequencePPOnAAForLowPt', 'HLTHIDoLocalPixelSequence']
    hasHISeq = any(seq in process.__dict__ for seq in hiSeq)

    if not (hasHISeq or 'HLTDoLocalPixelSequence' in process.__dict__):
        return process

    if not 'HLTDoLocalPixelSequence' in process.__dict__:
    	process.HLTDoLocalPixelSequence = cms.Sequence()
    	process.hltSiPixelClusters = fragment.hltSiPixelClusters.clone()
    	process.hltSiPixelClustersCache = fragment.hltSiPixelClustersCache.clone()
    	process.hltSiPixelRecHits = fragment.hltSiPixelRecHits.clone()

    process = customisePixelLocalReconstruction(process)

    # apply customisation for PbPb collisions
    # maximum fed for phase1 is 150 and max word for HIon is 8000
    process.hltSiPixelClustersCUDA.MaxFEDWords = cms.uint32(150 * 8000)
    pp_on_AA_2018.toModify(process.hltSiPixelClustersCUDA, isRun2 = True)

    if not hasHISeq:
        return process

    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    for seq in hiSeq:
        if seq in process.__dict__:
            setattr(process, seq, cms.Sequence())


    # Modules and EDAliases

    # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
    if 'HLTDoLocalPixelSequencePPOnAA' in process.__dict__:
    	process.hltSiPixelClustersPPOnAA = process.hltSiPixelClusters.clone(cpu = process.hltSiPixelClustersPPOnAA)
    if 'HLTDoLocalPixelSequencePPOnAAForLowPt' in process.__dict__:
    	process.hltSiPixelClustersPPOnAAForLowPt = process.hltSiPixelClusters.clone(cpu = process.hltSiPixelClustersPPOnAAForLowPt)
    if 'HLTHIDoLocalPixelSequence' in process.__dict__:
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
    if 'HLTDoLocalPixelSequencePPOnAA' in process.__dict__:
    	process.hltSiPixelRecHitsPPOnAA = process.hltSiPixelRecHits.clone(cpu = process.hltSiPixelRecHitsPPOnAA)
    	process.hltSiPixelRecHitsPPOnAA.cuda.src = "hltSiPixelClustersPPOnAA"
    if 'HLTDoLocalPixelSequencePPOnAAForLowPt' in process.__dict__:
    	process.hltSiPixelRecHitsPPOnAAForLowPt = process.hltSiPixelRecHits.clone(cpu = process.hltSiPixelRecHitsPPOnAAForLowPt)
    	process.hltSiPixelRecHitsPPOnAAForLowPt.cuda.src = "hltSiPixelClustersPPOnAAForLowPt"
    if 'HLTHIDoLocalPixelSequence' in process.__dict__:
    	process.hltHISiPixelRecHits = process.hltSiPixelRecHits.clone(cpu = process.hltHISiPixelRecHits)
    	process.hltHISiPixelRecHits.cuda.src = "hltHISiPixelClusters"


    # Tasks and Sequences

    if 'HLTDoLocalPixelSequencePPOnAA' in process.__dict__:
    	process.HLTDoLocalPixelTaskPPOnAA = cms.Task(
              process.HLTDoLocalPixelTask,
              process.hltSiPixelClustersPPOnAA,                 # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
              process.hltSiPixelClustersCachePPOnAA,            # legacy module, used by the legacy pixel quadruplet producer
              process.hltSiPixelRecHitsPPOnAA)                  # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA

    	process.HLTDoLocalPixelSequencePPOnAA = cms.Sequence(process.HLTDoLocalPixelTaskPPOnAA)

    if 'HLTDoLocalPixelSequencePPOnAAForLowPt' in process.__dict__:
    	process.HLTDoLocalPixelTaskPPOnAAForLowPt = cms.Task(
              process.HLTDoLocalPixelTask,
              process.hltSiPixelClustersPPOnAAForLowPt,         # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
              process.hltSiPixelClustersCachePPOnAAForLowPt,    # legacy module, used by the legacy pixel quadruplet producer
              process.hltSiPixelRecHitsPPOnAAForLowPt)          # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA

    	process.HLTDoLocalPixelSequencePPOnAAForLowPt = cms.Sequence(process.HLTDoLocalPixelTaskPPOnAAForLowPt)

    if 'HLTHIDoLocalPixelSequence' in process.__dict__:
    	process.HLTHIDoLocalPixelSequence = cms.Sequence(
              process.hltHISiPixelDigis*                        # SwitchProducer wrapping the legacy pixel digis producer or an alias combining the pixel digis information converted from SoA
              process.hltHISiPixelClusters*                     # SwitchProducer wrapping the legacy pixel cluster producer or an alias for the pixel clusters information converted from SoA
              process.hltHISiPixelClustersCache*                # legacy module, used by the legacy pixel quadruplet producer
              process.hltHISiPixelRecHits,     		        # SwitchProducer wrapping the legacy pixel rechit producer or the transfer of the pixel rechits to the host and the conversion from SoA
              process.HLTDoLocalPixelTask)


    # done
    return process


# customisation for running the "Patatrack" pixel track reconstruction
def customisePixelTrackReconstructionHIon(process):

    hiHLTPixelTracksSeq = ['HLTRecoPixelTracksPPOnAASequence', 'HLTRecoPixelTracksSequencePPOnAA']
    hasHIHLTPixelTracksSeq = any(seq in process.__dict__ for seq in hiHLTPixelTracksSeq)

    if not (hasHIHLTPixelTracksSeq or 'HLTRecoPixelTracksSequence' in process.__dict__):
        return process

    if not 'HLTRecoPixelTracksSequence' in process.__dict__:
    	process.HLTRecoPixelTracksSequence = cms.Sequence()
    	process.hltPixelTracksTrackingRegions = fragment.hltPixelTracksTrackingRegions.clone()

    hiHLTPixelVertexRecoSeq = ['HLTRecopixelvertexingSequencePPOnAA', 'HLTPixelVertexingPPOnAASequence', 'HLTPixelVertexingSequencePPOnAA']
    hasHIHLTPixelVertexRecoSeq = any(seq in process.__dict__ for seq in hiHLTPixelVertexRecoSeq)

    if hasHIHLTPixelVertexRecoSeq and not 'HLTRecopixelvertexingSequence' in process.__dict__:
        process.HLTRecopixelvertexingSequence = cms.Sequence()
        process.hltTrimmedPixelVertices = fragment.hltTrimmedPixelVertices.clone()

    process = customisePixelTrackReconstruction(process)

    # apply customisation for PbPb collisions
    pp_on_AA_2018.toModify(process.hltPixelTracksCUDA, idealConditions = False)
    pp_on_AA_2018.toModify(process.hltPixelTracksSoA.cpu, idealConditions = False)

    if not hasHIHLTPixelTracksSeq:
        return process

    # FIXME replace the Sequences with empty ones to avoid exanding them during the (re)definition of Modules and EDAliases

    for seq in hiHLTPixelTracksSeq:
        if seq in process.__dict__:
            setattr(process, seq, cms.Sequence())

    for seq in hiHLTPixelVertexRecoSeq:
    	if seq in process.__dict__:
            setattr(process, seq, cms.Sequence())


    # Modules and EDAliases

    # referenced in process.HLTRecoPixelTracksTask

    # cpu only: convert the pixel rechits from legacy to SoA format
    process.hltSiPixelRecHitSoA.src = "hltSiPixelClustersPPOnAA"

    # convert the pixel tracks from SoA to legacy format
    process.hltPixelTracksPPOnAA = process.hltPixelTracks.clone(pixelRecHitLegacySrc = "hltSiPixelRecHitsPPOnAA")


    # referenced in process.HLTRecopixelvertexingTask
    if hasHIHLTPixelVertexRecoSeq:

    	# convert the pixel vertices from SoA to legacy format
    	process.hltPixelVerticesPPOnAA = process.hltPixelVertices.clone(TrackCollection = "hltPixelTracksPPOnAA")


    # Tasks and Sequences

    process.HLTRecoPixelTracksTask = cms.Task(
          process.hltPixelTracksTrackingRegionsPPOnAA,      # from the original sequence
          process.hltSiPixelRecHitSoA,                      # pixel rechits on cpu, converted to SoA
          process.hltPixelTracksCUDA,                       # pixel ntuplets on gpu, in SoA format
          process.hltPixelTracksSoA,                        # pixel ntuplets on cpu, in SoA format
          process.hltPixelTracksPPOnAA)                     # pixel tracks on cpu, in legacy format

    if 'HLTRecoPixelTracksSequencePPOnAA' in process.__dict__:
    	process.HLTRecoPixelTracksSequencePPOnAA = cms.Sequence(process.HLTRecoPixelTracksTask)
    if 'HLTRecoPixelTracksPPOnAASequence' in process.__dict__:
    	process.HLTRecoPixelTracksPPOnAASequence = cms.Sequence(process.HLTRecoPixelTracksTask)

    if hasHIHLTPixelVertexRecoSeq:
        process.HLTRecopixelvertexingTask = cms.Task(
              process.HLTRecoPixelTracksTask,
              process.hltPixelVerticesCUDA,                     # pixel vertices on gpu, in SoA format
              process.hltPixelVerticesSoA,                      # pixel vertices on cpu, in SoA format
              process.hltPixelVerticesPPOnAA,                   # pixel vertices on cpu, in legacy format
              process.hltTrimmedPixelVerticesPPOnAA)		# from the original sequence

        process.HLTRecopixelvertexingSequence = cms.Sequence(
              process.hltPixelTracksFitter +                    # not used here, kept for compatibility with legacy sequences
              process.hltPixelTracksFilter,                     # not used here, kept for compatibility with legacy sequences
              process.HLTRecopixelvertexingTask)

        if 'HLTRecopixelvertexingSequencePPOnAA' in process.__dict__:
            process.HLTRecopixelvertexingSequencePPOnAA = cms.Sequence(process.HLTRecopixelvertexingSequence)
        if 'HLTPixelVertexingSequencePPOnAA' in process.__dict__:
            process.HLTPixelVertexingSequencePPOnAA = cms.Sequence(process.HLTRecopixelvertexingSequence)
        if 'HLTPixelVertexingPPOnAASequence' in process.__dict__:
            process.HLTPixelVertexingPPOnAASequence = cms.Sequence(process.HLTRecopixelvertexingSequence)


    # done
    return process


# customisation for offloading the ECAL local reconstruction via CUDA if a supported gpu is present
def customiseEcalLocalReconstructionHIon(process):
    process = customiseEcalLocalReconstruction(process)

    if not any(seq in process.__dict__ for seq in ['HLTDoFullUnpackingEgammaEcalMFSequence', 'HLTDoFullUnpackingEgammaEcalSequence', 'HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence']):
        return process

    # apply customisation for PbPb collisions
    # maximum number of bytes per fed for HIon is 32 * 1024
    process.hltEcalDigisGPU.maxFedSize = cms.uint32(32 * 1024)


    # done
    return process


# customisation for offloading the HCAL local reconstruction via CUDA if a supported gpu is present
def customiseHcalLocalReconstructionHIon(process):
    process = customiseHcalLocalReconstruction(process)

    if not 'HLTDoLocalHcalWithTowerSequence' in process.__dict__:
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
    process = customiseEcalLocalReconstructionHIon(process)
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
