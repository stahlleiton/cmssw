import FWCore.ParameterSet.Config as cms

def customiseAddConversionsCollection(process):

    process.load("RecoHI.HiTracking.hiMixedTripletStep_cff")
    process.load("RecoHI.HiTracking.hiPixelLessStep_cff")
    process.load("RecoHI.HiTracking.hiTobTecStep_cff")

## ----------------------------------------------------------------------------- ##
##
## HI Version of : RecoTracker/ConversionSeedGenerators/python/ConversionStep_cff.py
##

    process.hiConvClusters = cms.EDProducer("HITrackClusterRemover",
        clusterLessSolution = cms.bool(True),
        trajectories = cms.InputTag("hiTobTecStepTracks"),
        overrideTrkQuals = cms.InputTag('hiTobTecStepSelector','hiTobTecStep'),
        oldClusterRemovalInfo = cms.InputTag("hiTobTecStepClusters"),
        TrackQuality = cms.string('highPurity'),
        pixelClusters = cms.InputTag("siPixelClusters"),
        stripClusters = cms.InputTag("siStripClusters"),
        Common = cms.PSet(
            maxChi2 = cms.double(30.0),
            ),
        Strip = cms.PSet(
            #Yen-Jie's mod to preserve merged clusters
            maxSize = cms.uint32(2),
            maxChi2 = cms.double(30.0) #?????
            )
        )

    from RecoTracker.ConversionSeedGenerators.ConversionStep_cff import convLayerPairs
    process.hiConvLayerPairs = convLayerPairs.clone()
    process.hiConvLayerPairs.BPix.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.FPix.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TIB1.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TIB2.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TIB3.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TIB4.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TID1.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TID2.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TID3.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TEC.skipClusters  = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TOB1.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TOB2.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TOB3.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TOB4.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TOB5.skipClusters = cms.InputTag('hiConvClusters')
    process.hiConvLayerPairs.TOB6.skipClusters = cms.InputTag('hiConvClusters')


    from RecoTracker.ConversionSeedGenerators.PhotonConversionTrajectorySeedProducerFromSingleLeg_cfi import photonConvTrajSeedFromSingleLeg
    process.hiPhotonConvTrajSeedFromSingleLeg  = photonConvTrajSeedFromSingleLeg.clone(
        TrackRefitter = cms.InputTag('hiGeneralTracks'),
        primaryVerticesTag   = cms.InputTag("hiOfflinePrimaryVertices"),
        newSeedCandidates    = cms.string("hiConvSeedCandidates")
        )
    process.hiPhotonConvTrajSeedFromSingleLeg.OrderedHitsFactoryPSet.SeedingLayers = cms.InputTag('hiConvLayerPairs')
    process.hiPhotonConvTrajSeedFromSingleLeg.ClusterCheckPSet.cut = "strip < 400000 && pixel < 40000 && (strip < 60000 + 7.0*pixel) && (pixel < 8000 + 0.14*strip)"
    process.hiPhotonConvTrajSeedFromSingleLeg.SeedCreatorPSet.ComponentName = cms.string('SeedForHIPhotonConversion1Leg')

    
    # TRACKER DATA CONTROL
    
    # QUALITY CUTS DURING TRACK BUILDING
    import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff
    process.hiConvCkfTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
        maxLostHits = 1,
        minimumNumberOfHits = 3,
        minPt = 0.1
        )
    
    import RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi
    process.hiConvStepChi2Est = RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi.Chi2ChargeMeasurementEstimator.clone(
        ComponentName = cms.string('hiConvStepChi2Est'),
        nSigma = cms.double(3.0),
        MaxChi2 = cms.double(30.0),
        MaxDisplacement = cms.double(100),
        MaxSagitta = cms.double(-1.),
        clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight'))
        )
    

    # TRACK BUILDING
    import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
    process.hiConvCkfTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
        trajectoryFilter = cms.PSet(refToPSet_ = cms.string('hiConvCkfTrajectoryFilter')),
        minNrOfHitsForRebuild = 3,
        maxCand = 1,
        estimator = cms.string('hiConvStepChi2Est')
        )
    
    # MAKING OF TRACK CANDIDATES
    import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
    process.hiConvTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
        src = cms.InputTag('hiPhotonConvTrajSeedFromSingleLeg:hiConvSeedCandidates'),
        clustersToSkip = cms.InputTag('hiConvClusters'),
        TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('hiConvCkfTrajectoryBuilder'))
        )
    
    import TrackingTools.TrackFitters.RungeKuttaFitters_cff
    process.hiConvStepFitterSmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.KFFittingSmootherWithOutliersRejectionAndRK.clone(
        ComponentName = 'hiConvStepFitterSmoother',
        EstimateCut = 30,
        Smoother = cms.string('hiConvStepRKSmoother')
        )
    
    process.hiConvStepRKTrajectorySmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectorySmoother.clone(
        ComponentName = cms.string('hiConvStepRKSmoother'),
        errorRescaling = 10.0
        )
    
    
    # TRACK FITTING
    import RecoTracker.TrackProducer.TrackProducer_cfi
    process.hiConvStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
        src = 'hiConvTrackCandidates',
        AlgorithmName = cms.string('hiConversionStep'),
        Fitter = 'hiConvStepFitterSmoother',
        )


    # FINAL SELECTION
    import RecoHI.HiTracking.hiMultiTrackSelector_cfi
    process.hiConvStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
        src='hiConvStepTracks',
        useAnyMVA = cms.bool(True),
        GBRForestLabel = cms.string('HIMVASelectorIter13'),
        GBRForestVars = cms.vstring(['chi2perdofperlayer', 'nhits', 'nlayers', 'eta']),
        trackSelectors= cms.VPSet(
            RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiLooseMTS.clone(
                name = 'hiConvStepLoose',
                applyAdaptedPVCuts = cms.bool(False),
                useMVA = cms.bool(False),
                ), #end of pset
            RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
                name = 'hiConvStepTight',
                preFilterName = 'hiConvStepLoose',
                applyAdaptedPVCuts = cms.bool(False),
                useMVA = cms.bool(False),
                minMVA = cms.double(-0.2)
                ),
            RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
                name = 'hiConvStep',
                preFilterName = 'hiConvStepTight',
                applyAdaptedPVCuts = cms.bool(False),
                useMVA = cms.bool(False),
                minMVA = cms.double(-0.09)
                ),
            ) #end of vpset
        ) #end of clone


    import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
    process.hiConversionStepTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
        TrackProducers = cms.VInputTag(cms.InputTag('hiConvStepTracks')),
        hasSelector=cms.vint32(1),
        selectedTrackQuals = cms.VInputTag(cms.InputTag("hiConvStepSelector","hiConvStep")),
        copyExtras = True,
        makeReKeyedSeeds = cms.untracked.bool(False),
        setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(1), pQual=cms.bool(True) ) ),
        )

    # SEQUENCE
    process.hiConvStep = cms.Sequence(process.hiConvClusters 
                                      * process.hiConvLayerPairs
                                      * process.hiPhotonConvTrajSeedFromSingleLeg 
                                      * process.hiConvTrackCandidates
                                      * process.hiConvStepTracks
                                      * process.hiConvStepSelector
                                      * process.hiConversionStepTracks
                                      )


## ----------------------------------------------------------------------------- ##
##
## HI Version of : RecoEgamma/EgammaPhotonProducers/python/conversionTracks_cff.py
##

    # Conversion Track candidate producer 
    process.load("RecoEgamma.EgammaPhotonProducers.trajectoryBuilderForConversions_cfi")
    process.load("RecoEgamma.EgammaPhotonProducers.trajectoryFilterForConversions_cfi")
    from RecoEgamma.EgammaPhotonProducers.conversionTrackCandidates_cfi import conversionTrackCandidates
    process.hiConversionTrackCandidates = conversionTrackCandidates.clone(
        inOutTrackCandidateCollection = cms.string('inOutTracksFromConversions'),
        outInTrackCandidateCollection = cms.string('outInTracksFromConversions')
        )
    # Conversion Track producer  ( final fit )
    from RecoEgamma.EgammaPhotonProducers.ckfOutInTracksFromConversions_cfi import ckfOutInTracksFromConversions
    process.hiCkfOutInTracksFromConversions = ckfOutInTracksFromConversions.clone(
        src = cms.InputTag("hiConversionTrackCandidates","outInTracksFromConversions"),
        producer = cms.string('hiConversionTrackCandidates'),
        ComponentName = cms.string('hiCkfOutInTracksFromConversions')
        )

    from RecoEgamma.EgammaPhotonProducers.ckfInOutTracksFromConversions_cfi import ckfInOutTracksFromConversions
    process.hiCkfInOutTracksFromConversions = ckfInOutTracksFromConversions.clone(
        src = cms.InputTag("hiConversionTrackCandidates","inOutTracksFromConversions"),
        producer = cms.string('hiConversionTrackCandidates'),
        ComponentName = cms.string('hiCkfInOutTracksFromConversions')
        )

    process.hiCkfTracksFromConversions = cms.Sequence(process.hiConversionTrackCandidates
                                                      * process.hiCkfOutInTracksFromConversions
                                                      * process.hiCkfInOutTracksFromConversions
                                                      )

## ----------------------------------------------------------------------------- ##
##
## HI Version of : RecoEgamma/EgammaPhotonProducers/python/conversionTrackSequence_cff.py
##

    ## Conversion Prducer
    import RecoEgamma.EgammaPhotonProducers.conversionTrackProducer_cfi
    
    # producer from general tracks collection, set tracker only, merged arbitrated, merged arbitrated ecal/general flags
    process.hiGeneralConversionTrackProducer = RecoEgamma.EgammaPhotonProducers.conversionTrackProducer_cfi.conversionTrackProducer.clone(
        TrackProducer = cms.string('hiGeneralTracks'),
        setTrackerOnly = cms.bool(True),
        setArbitratedMergedEcalGeneral = cms.bool(True),
        )

    # producer from conversionStep tracks collection, set tracker only, merged arbitrated, merged arbitrated ecal/general flags
    process.hiConversionStepConversionTrackProducer = RecoEgamma.EgammaPhotonProducers.conversionTrackProducer_cfi.conversionTrackProducer.clone(
        TrackProducer = cms.string('hiConversionStepTracks'),
        setTrackerOnly = cms.bool(True),
        setArbitratedMergedEcalGeneral = cms.bool(True),
        )

    # producer from inout ecal seeded tracks, set arbitratedecalseeded, mergedarbitratedecalgeneral and mergedarbitrated flags
    process.hiInOutConversionTrackProducer = RecoEgamma.EgammaPhotonProducers.conversionTrackProducer_cfi.conversionTrackProducer.clone(
        TrackProducer = cms.string('hiCkfInOutTracksFromConversions'),
        setArbitratedEcalSeeded = cms.bool(True),
        setArbitratedMergedEcalGeneral = cms.bool(True),    
        )

    # producer from outin ecal seeded tracks, set arbitratedecalseeded, mergedarbitratedecalgeneral and mergedarbitrated flags
    process.hiOutInConversionTrackProducer = RecoEgamma.EgammaPhotonProducers.conversionTrackProducer_cfi.conversionTrackProducer.clone(
        TrackProducer = cms.string('hiCkfOutInTracksFromConversions'),
        setArbitratedEcalSeeded = cms.bool(True),
        setArbitratedMergedEcalGeneral = cms.bool(True),    
        )

    # producer from gsf tracks, set only mergedarbitrated flag (default behaviour)
    process.gsfConversionTrackProducer = RecoEgamma.EgammaPhotonProducers.conversionTrackProducer_cfi.conversionTrackProducer.clone(
        TrackProducer = cms.string('electronGsfTracks'),
        filterOnConvTrackHyp = cms.bool(False),
        )

    # SEQUENCE
    process.hiConversionTrackProducers = cms.Sequence(process.hiGeneralConversionTrackProducer
                                                    * process.hiConversionStepConversionTrackProducer
                                                    * process.hiInOutConversionTrackProducer
                                                    * process.hiOutInConversionTrackProducer
                                                    * process.gsfConversionTrackProducer
                                                    )

    ## Conversion Merger
    import RecoEgamma.EgammaPhotonProducers.conversionTrackMerger_cfi

    # merge generalTracks and conversionStepTracks collections, with arbitration by nhits then chi^2/ndof for ecalseededarbitrated, mergedarbitratedecalgeneral and mergedarbitrated flags
    process.hiGeneralConversionStepConversionTrackMerger = RecoEgamma.EgammaPhotonProducers.conversionTrackMerger_cfi.conversionTrackMerger.clone(
        TrackProducer1 = cms.InputTag('hiGeneralConversionTrackProducer'),
        TrackProducer2 = cms.InputTag('hiGonversionStepConversionTrackProducer'),
        # 3: arbitrate output/flag (remove duplicates by shared hits), arbitration first by number of hits, second by chisq/ndof  
        arbitratedMergedPreferCollection = cms.int32(3),
        arbitratedMergedEcalGeneralPreferCollection = cms.int32(3),        
        )

    # merge two ecal-seeded collections, with arbitration by nhits then chi^2/ndof for ecalseededarbitrated, mergedarbitratedecalgeneral and mergedarbitrated flags
    process.hiInOutOutInConversionTrackMerger = RecoEgamma.EgammaPhotonProducers.conversionTrackMerger_cfi.conversionTrackMerger.clone(
        TrackProducer1 = cms.InputTag('hiInOutConversionTrackProducer'),
        TrackProducer2 = cms.InputTag('hiOutInConversionTrackProducer'),
        # 3: arbitrate output/flag (remove duplicates by shared hits), arbitration first by number of hits, second by chisq/ndof  
        arbitratedEcalSeededPreferCollection = cms.int32(3),
        arbitratedMergedPreferCollection = cms.int32(3),
        arbitratedMergedEcalGeneralPreferCollection = cms.int32(3),        
        )

    process.hiGeneralInOutOutInConversionTrackMerger = RecoEgamma.EgammaPhotonProducers.conversionTrackMerger_cfi.conversionTrackMerger.clone(
        TrackProducer1 = cms.InputTag('hiInOutOutInConversionTrackMerger'),
        TrackProducer2 = cms.InputTag('hiGeneralConversionStepConversionTrackMerger'),
        arbitratedMergedPreferCollection = cms.int32(3),
        arbitratedMergedEcalGeneralPreferCollection = cms.int32(2),        
        )

    process.hiGsfGeneralInOutOutInConversionTrackMerger = RecoEgamma.EgammaPhotonProducers.conversionTrackMerger_cfi.conversionTrackMerger.clone(
        TrackProducer1 = cms.InputTag('hiGeneralInOutOutInConversionTrackMerger'),
        TrackProducer2 = cms.InputTag('gsfConversionTrackProducer'),
        arbitratedMergedPreferCollection = cms.int32(2),
        )

    # final output collection contains combination of generaltracks, ecal seeded tracks and gsf tracks, with overlaps removed by shared hits
    # precedence is given first to gsf tracks, then to the combination of ecal seeded and general tracks
    # overlaps between the ecal seeded track collections and between ecal seeded and general tracks are arbitrated first by nhits then by chi^2/dof
    # (logic and much of the code is adapted from FinalTrackSelectors)
    
    process.hiConversionTrackMergers = cms.Sequence(process.hiInOutOutInConversionTrackMerger
                                                    * process.hiGeneralConversionStepConversionTrackMerger
                                                    * process.hiGeneralInOutOutInConversionTrackMerger
                                                    * process.hiGsfGeneralInOutOutInConversionTrackMerger
                                                    )


## ----------------------------------------------------------------------------- ##
##
## CONVERSION TRACK SEQUENCE
##

    process.hiConversionTrackSequence = cms.Sequence(process.hiCkfTracksFromConversions
                                                   * process.hiConversionTrackProducers
                                                   * process.hiConversionTrackMergers
                                                   )


## ----------------------------------------------------------------------------- ##
##
## HI Version of : RecoEgamma/EgammaPhotonProducers/python/allConversions_cfi.py
##

    from RecoEgamma.EgammaPhotonProducers.allConversions_cfi import allConversions
    process.hiAllConversions = allConversions.clone(
        src = cms.InputTag("hiGsfGeneralInOutOutInConversionTrackMerger"),
        primaryVertexProducer = cms.InputTag("hiOfflinePrimaryVertices")
        )





    process.hiConversionSequence = cms.Sequence(process.hiMixedTripletStep
                                                * process.hiPixelLessStep
                                                * process.hiTobTecStep 
                                                * process.hiConvStep 
                                                * process.hiConversionTrackSequence 
                                                * process.hiAllConversions
                                                )

    ###Add conversions in the sequence
    process.reconstruction_step += process.hiConversionSequence


    ## Include in event content
    for p in ["RAWAODSIMoutput" , "AODSIMoutput" , "RECOoutput" , "AODoutput"] :
        if hasattr(process, p) : getattr(process, p).outputCommands.extend(['keep *_hiAllConversions_*_*'])


    return process
