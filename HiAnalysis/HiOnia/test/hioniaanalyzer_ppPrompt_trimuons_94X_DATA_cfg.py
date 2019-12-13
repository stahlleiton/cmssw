import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


#----------------------------------------------------------------------------

# Setup Settings for ONIA TREE: PbPb 2018

HLTProcess     = "HLT" # Name of HLT process 
isMC           = False # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "TwoGlbAmongThree" # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker), TwoGlbAmongThree (which requires two isGlobal for a trimuon, and one isGlobal for a dimuon) are available
applyEventSel  = True # Only apply Event Selection if the required collections are present
OnlySoftMuons  = True # Keep only isSoftMuon's (with highPurity because this is pp config) from the beginning of HiSkim. In any case, if applyCuts=True, isSoftMuon is required at HiAnalysis level for muons of selected dimuons.
applyCuts      = True # At HiAnalysis level, apply kinematic acceptance cuts + identification cuts (isSoftMuon or isTightMuon, depending on TightGlobalMuon flag) for muons from selected di(tri)muons + hard-coded cuts on the di(tri)muon that you would want to add (but recommended to add everything in LateDimuonSelection, applied at the end of HiSkim)
SofterSgMuAcceptance = True # Whether to accept muons with a softer acceptance cuts than the usual (pt>3.5GeV at central eta, pt>1.8 at high |eta|). Applies when applyCuts=True
doTrimuons     = True # Make collections of trimuon candidates in addition to dimuons, and keep only events with >0 trimuons
doDimuonTrk    = False # Make collections of Jpsi+track candidates in addition to dimuons
atLeastOneCand = True # Keep only events that have one selected dimuon (or at least one trimuon if doTrimuons = true). BEWARE this can cause trouble in .root output if no event is selected by onia2MuMuPatGlbGlbFilter!
OneMatchedHLTMu = 4   # Keep only di(tri)muons of which the one(two) muon(s) are matched to the HLT Filter of this number. You can get the desired number in the output of oniaTree. Set to -1 for no matching. WARNING: it is the trigger bit+1 !
muonLessPV     = True  # Recalculate the PV without the two muons from the selected dimuon
#############################################################################
keepExtraColl  = False # General Tracks + Stand Alone Muons + Converted Photon collections
useSVfinder    = False # External SV finder to check if the muons are from a resolved SV
#saveHLTBit     = False # for trigger analysis
#saveHLTobj     = False # For trigger analysis
#----------------------------------------------------------------------------

# Print Onia Tree settings:
print( " " )
print( "[INFO] Settings used for ONIA TREE: " )
print( "[INFO] isMC                 = " + ("True" if isMC else "False") )
print( "[INFO] applyEventSel        = " + ("True" if applyEventSel else "False") )
print( "[INFO] keepExtraColl        = " + ("True" if keepExtraColl else "False") )
print( "[INFO] SofterSgMuAcceptance = " + ("True" if SofterSgMuAcceptance else "False") )
print( "[INFO] muonSelection        = " + muonSelection )
print( "[INFO] onlySoftMuons        = " + ("True" if OnlySoftMuons else "False") )
print( "[INFO] doTrimuons           = " + ("True" if doTrimuons else "False") )
print( " " )

# set up process
process = cms.Process("HIOnia")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "Oniatree.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles =[#'file:/home/llr/cms/falmagne/production/pp2017/BcTrimu/CMSSW_9_4_14/src/pp5TeV_TrimuonSkim.root',
    '/store/user/gfalmagn/AOD/DoubleMu_Run2017G_AOD_Run_306546_306826_trimuonSkim_05082019/DoubleMuon/crab_DoubleMu_Run2017G_AOD_Run_306546_306826_trimuonSkim_05082019/190805_182650/0000/pp5TeV_TrimuonSkim_288.root',
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/90000/B6F95CE3-5B2E-E811-9F6F-A0369FD20D28.root',
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/90000/B68C2814-912D-E811-9A7E-A0369FC5252C.root',
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/90000/B2518E21-2538-E811-B08A-0025905C9726.root',
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/90000/B03BCD95-892E-E811-B8C5-A0369FC513DC.root',
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/90000/AC904279-852E-E811-AD85-5C260AFFFC17.root',
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/90000/AC288DD0-552E-E811-A74F-3417EBE52915.root'
    #'/store/data/Run2017G/DoubleMuon/AOD/17Nov2017-v1/00000/E2A9F70B-8042-E811-9D04-FA163E74586C.root'
    ]
options.maxEvents = 100 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

triggerList    = {
    # Double Muon Trigger List
    'DoubleMuonTrigger' : cms.vstring(
        "HLT_HIL1DoubleMuOpen_v1",
        "HLT_HIL1DoubleMuOpen_OS_v1",
        "HLT_HIL1DoubleMuOpen_SS_v1",
        "HLT_HIL1DoubleMu0_v1",
        "HLT_HIL1DoubleMu0_HighQ_v1",
        "HLT_HIL1DoubleMu10_v1",
        "HLT_HIL2DoubleMu0_v1",
        "HLT_HIL2DoubleMu10_v1",
        "HLT_HIL3DoubleMu0_v1",
        "HLT_HIL3DoubleMu10_v1",
        ),
    # Double Muon Filter List
    'DoubleMuonFilter'  : cms.vstring(
        "hltL1fL1sDoubleMuOpenL1Filtered0",
        "hltL1fL1sDoubleMuOpenOSL1Filtered0",
        "hltL1fL1sDoubleMuOpenSSL1Filtered0",
        "hltL1fL1sDoubleMu0L1Filtered0",
        "hltL1fL1sDoubleMu0L1HighQFiltered0",
        "hltL1fL1sDoubleMu10L1Filtered0",
        "hltL2fL1sDoubleMu0L1f0L2Filtered0",
        "hltL2fL1sDoubleMu10L1f0L2Filtered10",
        "hltL3fL1sDoubleMu0L1f0L2f0L3Filtered0",
        "hltL3fL1sDoubleMu10L1f0L2f0L3Filtered10",
        ),
    # Single Muon Trigger List
    'SingleMuonTrigger' : cms.vstring(
        "HLT_HIL1Mu12_v1",
        "HLT_HIL1Mu16_v1",
        "HLT_HIL2Mu7_v1",
        "HLT_HIL2Mu12_v1",
        "HLT_HIL2Mu15_v1",
        "HLT_HIL2Mu20_v1",
        "HLT_HIL3Mu3_v1",
        "HLT_HIL3Mu5_v1",
        "HLT_HIL3Mu7_v1",
        "HLT_HIL3Mu12_v1",
        "HLT_HIL3Mu15_v1",
        "HLT_HIL3Mu20_v1",
        "HLT_HIL2Mu3_NHitQ10_v1",
        "HLT_HIL3Mu3_NHitQ10_v1",
        "HLT_HIL2Mu5_NHitQ10_v1",
        "HLT_HIL3Mu5_NHitQ10_v1",
        ),
    # Single Muon Filter List
    'SingleMuonFilter'  : cms.vstring(
        "hltL1fL1sSingleMu12L1Filtered0",
        "hltL1fL1sSingleMu16L1Filtered0",
        "hltL2fL1sSingleMu3OR5L1f0L2Filtered7",
        "hltL2fL1sSingleMu7L1f0L2Filtered12",
        "hltL2fL1sSingleMu7L1f0L2Filtered15",
        "hltL2fL1sSingleMu7L1f0L2Filtered20",
        "hltL3fL1sSingleMu3L1f0L2f0L3Filtered3",
        "hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered5",
        "hltL3fL1sSingleMu3OR5L1f0L2f0L3Filtered7",
        "hltL3fL1sSingleMu7L1f0L2f0L3Filtered12",
        "hltL3fL1sSingleMu7L1f0L2f0L3Filtered15",
        "hltL3fL1sSingleMu7L1f0L2f0L3Filtered20",
        "hltL2fL1sSingleMu3L1f0L2NHitQ10L2Filtered3",
        "hltL3fL1sSingleMu3L1f0L2f0L3NHitQ10L3Filtered3",
        "hltL2fL1sSingleMu3OR5L1f0L2NHitQ10L2Filtered5",
        "hltL3fL1sSingleMu3OR5L1f0L2f0L3NHitQ10L3Filtered5",
        )
}

if isMC:
    globalTag = '92X_upgrade2017_realistic_v11'
else:
    globalTag = '94X_dataRun2_ReReco_EOY17_v6'

# load the Geometry and Magnetic Field
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
    
# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')
    
#----------------------------------------------------------------------------

# For OniaTree Analyzer
from HiAnalysis.HiOnia.oniaTreeAnalyzer_cff import oniaTreeAnalyzer
oniaTreeAnalyzer(process,
                 #muonTriggerList=triggerList, HLTProName=HLTProcess, #useL1Stage2=True, 
                 muonSelection=muonSelection, isMC=isMC, outputFileName=options.outputFile, muonlessPV=muonLessPV, doTrimu=doTrimuons)
process.oniaTreeAna = cms.Path(process.oniaTreeAna)

process.onia2MuMuPatGlbGlb.dimuonSelection       = cms.string("2.2 < mass && mass < 4.0 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
process.onia2MuMuPatGlbGlb.trimuonSelection      = cms.string("2.9 < mass && mass < 8.3 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
#process.onia2MuMuPatGlbGlb.lowerPuritySelection  = cms.string("")
#process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("") ## No need to repeat lowerPuritySelection in there, already included 
if applyCuts:
    process.onia2MuMuPatGlbGlb.LateDimuonSel         = cms.string("userFloat(\"vProb\")>0.002")
    process.onia2MuMuPatGlbGlb.LateTrimuonSel        = cms.string("userFloat(\"vProb\")>0.005 && userFloat(\"ppdlPV3D\")>0 && userFloat(\"ppdlPV\")>0 && userFloat(\"cosAlpha\")>0.2")
process.onia2MuMuPatGlbGlb.onlySoftMuons         = cms.bool(OnlySoftMuons)

process.hionia.minimumFlag      = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.useGeTracks      = cms.untracked.bool(keepExtraColl) #for Reco_trk_*
process.hionia.fillRecoTracks   = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.srcTracks        = cms.InputTag("generalTracks")
process.hionia.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")
process.hionia.genParticles     = cms.InputTag("genParticles")
process.hionia.SofterSgMuAcceptance = cms.bool(SofterSgMuAcceptance)
process.hionia.applyCuts        = cms.bool(applyCuts)
process.hionia.AtLeastOneCand   = cms.bool(atLeastOneCand)
process.hionia.OneMatchedHLTMu  = cms.int32(OneMatchedHLTMu)
process.hionia.useSVfinder      = cms.bool(useSVfinder)

process.NoScraping = cms.EDFilter("FilterOutScraping",
                          applyfilter = cms.untracked.bool(True),
                          debugOn = cms.untracked.bool(False),
                          numtrack = cms.untracked.uint32(10),
                          thresh = cms.untracked.double(0.25)
                          )
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                   src = cms.InputTag("offlinePrimaryVertices"),
                                   cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"), 
                                   filter = cms.bool(True),   # otherwise it won't filter the events
                                   )

if applyEventSel:
    process.oniaTreeAna.replace(process.hionia, process.primaryVertexFilter * process.NoScraping * process.hionia )

if atLeastOneCand:
  if doTrimuons:
      process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilterTrimu)
      process.oniaTreeAna.replace(process.patMuonSequence, process.filter3mu * process.pseudoDimuonFilterSequence * process.patMuonSequence)
  elif doDimuonTrk:
      process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilterDimutrk)
      process.oniaTreeAna.replace(process.patMuonSequence, process.pseudoDimuonFilterSequence * process.patMuonSequence)
  else:
      process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter)
      #BEWARE, pseudoDimuonFilterSequence asks for opposite-sign dimuon in given mass range. But saves a lot of time by filtering before running PAT muons
      process.oniaTreeAna.replace(process.patMuonSequence, process.pseudoDimuonFilterSequence * process.patMuonSequence)


if useSVfinder:
    from RecoVertex.AdaptiveVertexFinder.inclusiveVertexFinder_cfi import *
    from RecoVertex.AdaptiveVertexFinder.vertexMerger_cfi import *
    from RecoVertex.AdaptiveVertexFinder.trackVertexArbitrator_cfi import *
    
    process.inclusiveVertexFinderLoose = inclusiveVertexFinder.clone(
        vertexMinDLen2DSig = 1.5,
        vertexMinDLenSig = 1.25,
        vertexMinAngleCosine = 0.001,
        maximumLongitudinalImpactParameter = 0.6, #default = 0.3
        maxNTracks = 10, #default = 30
        minPt = 1.2, #following muon acceptance
        #useVertexReco = False,
        #fitterSigmacut = 3.,
        clusterizer = cms.PSet(
            seedMax3DIPSignificance = cms.double(9999.),#default
            seedMax3DIPValue = cms.double(9999.),#default
            seedMin3DIPSignificance = cms.double(1.6), # default=1.2
            seedMin3DIPValue = cms.double(0.005), # default = 0.005
            clusterMaxDistance = cms.double(0.05),#default = 0.05
            clusterMaxSignificance = cms.double(3.),#default = 4.5
            distanceRatio = cms.double(10.),#default = 20
            clusterMinAngleCosine = cms.double(0.001), # default = 0.5
            maxTimeSignificance = cms.double(3.5),#default
        ),
    )
    
    process.vertexMergerLoose = vertexMerger.clone(
        secondaryVertices = "inclusiveVertexFinderLoose"
    )
    process.trackVertexArbitratorLoose = trackVertexArbitrator.clone(
        secondaryVertices = cms.InputTag("vertexMergerLoose")
    )
    process.inclusiveSecondaryVerticesLoose = vertexMerger.clone(
        secondaryVertices = "trackVertexArbitratorLoose",
        maxFraction = 0.2, # default 0.7 - 0.2
        minSignificance = 3. # default 2 - 10
    )
    process.inclusiveVertexingTask = cms.Task(
        process.inclusiveVertexFinderLoose,
        process.vertexMergerLoose,
        process.trackVertexArbitratorLoose,
        process.inclusiveSecondaryVerticesLoose
    )
    process.inclusiveVertexing = cms.Sequence(process.inclusiveVertexingTask)
    process.oniaTreeAna.replace(process.hionia, process.inclusiveVertexing*process.hionia)
    
    from RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi import *
    inclusiveSecondaryVertexFinderLooseTagInfos = secondaryVertexTagInfos.clone()
    # use external SV collection made from IVF
    inclusiveSecondaryVertexFinderLooseTagInfos.extSVCollection     = cms.InputTag('inclusiveSecondaryVerticesLoose')
    inclusiveSecondaryVertexFinderLooseTagInfos.useExternalSV = cms.bool(True)

#----------------------------------------------------------------------------
'''
# For HLTBitAnalyzer
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.HLTProcessName              = HLTProcess
process.hltbitanalysis.hltresults                  = cms.InputTag("TriggerResults","",HLTProcess)
process.hltbitanalysis.l1tAlgBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.l1tExtBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.gObjectMapRecord            = cms.InputTag("hltGtStage2ObjectMap","",HLTProcess)
process.hltbitanalysis.gmtStage2Digis              = cms.string("hltGtStage2Digis")
process.hltbitanalysis.caloStage2Digis             = cms.string("hltGtStage2Digis")
process.hltbitanalysis.UseL1Stage2                 = cms.untracked.bool(True)
process.hltbitanalysis.getPrescales                = cms.untracked.bool(False)
process.hltbitanalysis.getL1InfoFromEventSetup     = cms.untracked.bool(False)
process.hltbitanalysis.UseTFileService             = cms.untracked.bool(True)
process.hltbitanalysis.RunParameters.HistogramFile = cms.untracked.string(options.outputFile)
process.hltbitanalysis.RunParameters.isData        = cms.untracked.bool(not isMC)
process.hltbitanalysis.RunParameters.Monte         = cms.bool(isMC)
process.hltbitanalysis.RunParameters.GenTracks     = cms.bool(False)
if (HLTProcess == "HLT") :
	process.hltbitanalysis.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.gmtStage2Digis    = cms.string("gtStage2Digis")
	process.hltbitanalysis.caloStage2Digis   = cms.string("gtStage2Digis")
        if saveHLTBit:
          process.hltBitAna = cms.EndPath(process.hltbitanalysis)
#
##----------------------------------------------------------------------------
#
# For HLTObject Analyzer
process.load("HeavyIonsAnalysis.EventAnalysis.hltobject_cfi")
process.hltobject.processName = cms.string(HLTProcess)
process.hltobject.treeName = cms.string(options.outputFile)
process.hltobject.loadTriggersFromHLT = cms.untracked.bool(False)
process.hltobject.triggerNames = triggerList['DoubleMuonTrigger'] + triggerList['SingleMuonTrigger']
process.hltobject.triggerResults = cms.InputTag("TriggerResults","",HLTProcess)
process.hltobject.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","",HLTProcess)
if saveHLTobj:
  process.hltObjectAna = cms.EndPath(process.hltobject)
'''

#----------------------------------------------------------------------------
#Options:
process.source = cms.Source("PoolSource",
#process.source = cms.Source("NewEventStreamFileReader", # for streamer data
		fileNames = cms.untracked.vstring( options.inputFiles ),
		)
process.TFileService = cms.Service("TFileService", 
		fileName = cms.string( options.outputFile )
		)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#if saveHLTobj or saveHLTBit:
#    process.schedule  = cms.Schedule( process.oniaTreeAna , process.hltBitAna , process.hltObjectAna )
#else:
process.schedule  = cms.Schedule( process.oniaTreeAna )
