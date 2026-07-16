import FWCore.ParameterSet.Config as cms

def candidateBtaggingMiniAOD(process, isMC = True, jetPtMin = 15, jetR = 0.4, addNegTag = True):
    # DeepNtuple settings
    R = str(int(jetR*10))
    jetCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] + ([] if isMC else ['L2L3Residual'])
    jetCorrectionsAK = ('AK4PFchs', jetCorrLevels, 'None')

    bTagInfos = [
        'pfDeepCSVTagInfos',
        'pfDeepFlavourTagInfos',
        'pfImpactParameterTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'pfParticleTransformerAK4TagInfos',
        'pfUnifiedParticleTransformerAK4TagInfos'
    ]

    bTagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probudsg',
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:probc',
        'pfDeepFlavourJetTags:probg',
        'pfDeepFlavourJetTags:problepb',
        'pfDeepFlavourJetTags:probuds',
        'pfParticleTransformerAK4JetTags:probb',
        'pfParticleTransformerAK4JetTags:probbb',
        'pfParticleTransformerAK4JetTags:probc',
        'pfParticleTransformerAK4JetTags:probg',
        'pfParticleTransformerAK4JetTags:problepb',
        'pfParticleTransformerAK4JetTags:probuds',
    ]
    from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4JetTags_cfi import pfUnifiedParticleTransformerAK4JetTags as _UParTJetTags
    bTagDiscriminators += ['pfUnifiedParticleTransformerAK4JetTags:' + f for f in _UParTJetTags.flav_names]
    if addNegTag:
        bTagDiscriminators += ['pfNegativeOnlyJetProbabilityBJetTags']

    # Create gen-level information
    if isMC:
        process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
            src = cms.InputTag("packedGenParticles"),
            cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
        )
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        setattr(process,f'ak{R}GenJetsRecluster', ak4GenJets.clone(
            src = 'packedGenParticlesForJetsNoNu',
            rParam = jetR
        ))
        setattr(process,f'genAK{R}Task', cms.Task(process.packedGenParticlesForJetsNoNu, getattr(process,f'ak{R}GenJetsRecluster')))

    # Add negative secondary vertices
    process.svTask = cms.Task()
    if addNegTag:
        from RecoVertex.AdaptiveVertexFinder.inclusiveNegativeVertexing_cff import inclusiveCandidateNegativeVertexFinder, candidateNegativeVertexMerger, candidateNegativeVertexArbitrator, inclusiveCandidateNegativeSecondaryVertices
        process.inclusiveCandidateNegativeVertexFinder = inclusiveCandidateNegativeVertexFinder.clone(
            tracks = "packedPFCandidates",
            primaryVertices = "offlineSlimmedPrimaryVertices",
        )
        process.candidateNegativeVertexMerger = candidateNegativeVertexMerger.clone()
        process.candidateNegativeVertexArbitrator = candidateNegativeVertexArbitrator.clone(
            tracks = "packedPFCandidates",
            primaryVertices = "offlineSlimmedPrimaryVertices"
        )
        process.inclusiveCandidateNegativeSecondaryVertices = inclusiveCandidateNegativeSecondaryVertices.clone()
        for mod in ["inclusiveCandidateNegativeVertexFinder","candidateNegativeVertexMerger","candidateNegativeVertexArbitrator","inclusiveCandidateNegativeSecondaryVertices"]:
            process.svTask.add(getattr(process, mod))

    # Create CHS subtracted reco jets
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    addJetCollection(
        process,
        postfix            = "",
        labelName          = f"AK{R}PFJetsCHS",
        jetSource          = cms.InputTag(f"ak{R}PFJetsCHS"),
        algo               = "ak", #name of algo must be in this format
        rParam             = jetR,
        pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
        pfCandidates       = cms.InputTag("packedPFCandidates"),
        svSource           = cms.InputTag("slimmedSecondaryVertices"),
        muSource           = cms.InputTag("slimmedMuons"),
        elSource           = cms.InputTag("slimmedElectrons"),
        getJetMCFlavour    = isMC,
        genJetCollection   = cms.InputTag(f"ak{R}GenJetsRecluster" if isMC else ""),
        genParticles       = cms.InputTag("prunedGenParticles" if isMC else ""),
        jetCorrections     = jetCorrectionsAK,
    )
    getattr(process,f'patJetPartonMatchAK{R}PFJetsCHS').maxDeltaR = jetR

    if not isMC:
        for label in [f"patJetsAK{R}PFJetsCHS"]:
            getattr(process, label).addGenJetMatch = False
            getattr(process, label).addGenPartonMatch = False
            getattr(process, label).embedGenJetMatch = False
            getattr(process, label).embedGenPartonMatch = False
            getattr(process, label).genJetMatch = ""
            getattr(process, label).genPartonMatch = ""
    else:
        getattr(process,f'patJetPartonAssociationLegacyAK{R}PFJetsCHS').coneSizeToAssociate = min(jetR, 0.3)

    from CommonTools.ParticleFlow.pfCHS_cff import pfCHS
    process.pfCHS = pfCHS.clone()
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJetsCHS
    setattr(process,f'ak{R}PFJetsCHS', ak4PFJetsCHS.clone(
        src = 'pfCHS',
        rParam = jetR,
        jetPtMin = jetPtMin
    ))
    for mod in ['pfCHS',f"ak{R}PFJetsCHS"]:
        process.patAlgosToolsTask.add(getattr(process, mod))

    # Create b-tagging sequence ----------------
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        labelName = f"AK{R}PFJetsCHSDeepFlavour",
        jetSource = cms.InputTag(f'patJetsAK{R}PFJetsCHS'),
        jetCorrections = jetCorrectionsAK,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag("slimmedSecondaryVertices"),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
    )

    # Add UParT v2 model
    getattr(process,f"pfUnifiedParticleTransformerAK4JetTagsAK{R}PFJetsCHSDeepFlavour").model_path = 'RecoBTag/Combined/data/UParTAK4/PUPPI/V01/UParTAK4_v2.onnx'
    getattr(process,f"pfUnifiedParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour").fix_lt_sorting = True

    getattr(process,f'pfImpactParameterTagInfosAK{R}PFJetsCHSDeepFlavour').maxDeltaR = jetR
    taginfos = [f"pfDeepFlavourTagInfosAK{R}PFJetsCHSDeepFlavour", f"pfParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour", f"pfUnifiedParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour"]
    for taginfo in taginfos:
        getattr(process, taginfo).jet_radius = jetR

    if hasattr(process,f'updatedPatJetsTransientCorrectedAK{R}PFJetsCHSDeepFlavour'):
        getattr(process,f'updatedPatJetsTransientCorrectedAK{R}PFJetsCHSDeepFlavour').addTagInfos = True
        getattr(process,f'updatedPatJetsTransientCorrectedAK{R}PFJetsCHSDeepFlavour').addBTagInfo = True
    else:
        raise ValueError(f'I could not find updatedPatJetsTransientCorrectedAK{R}PFJetsCHSDeepFlavour to embed the tagInfos, please check the cfg')

    # Add extra b tagging algos
    from RecoBTag.ImpactParameter.pfJetProbabilityBJetTags_cfi import pfJetProbabilityBJetTags
    setattr(process,f'pfJetProbabilityBJetTagsAK{R}PFJetsCHSDeepFlavour', pfJetProbabilityBJetTags.clone(tagInfos = [f"pfImpactParameterTagInfosAK{R}PFJetsCHSDeepFlavour"]))
    process.patAlgosToolsTask.add(getattr(process,f'pfJetProbabilityBJetTagsAK{R}PFJetsCHSDeepFlavour'))

    #Add negative taggers
    if addNegTag:
        setattr(process,f'pfNegativeUnifiedParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour', getattr(process,f'pfUnifiedParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour').clone(
            flip = True,
            secondary_vertices = 'inclusiveCandidateNegativeSecondaryVertices',
        ))
        setattr(process,f'pfNegativeUnifiedParticleTransformerAK4JetTagsAK{R}PFJetsCHSDeepFlavour', getattr(process,f'pfUnifiedParticleTransformerAK4JetTagsAK{R}PFJetsCHSDeepFlavour').clone(
            src = f'pfNegativeUnifiedParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour',
        ))
        process.patAlgosToolsTask.add(getattr(process,f'pfNegativeUnifiedParticleTransformerAK4JetTagsAK{R}PFJetsCHSDeepFlavour'))
        process.patAlgosToolsTask.add(getattr(process,f'pfNegativeUnifiedParticleTransformerAK4TagInfosAK{R}PFJetsCHSDeepFlavour'))

    # Associate to forest sequence
    if isMC:
        process.forest.associate(getattr(process,f'genAK{R}Task'))
    process.forest.associate(process.svTask)
    process.forest.associate(process.patAlgosToolsTask)
