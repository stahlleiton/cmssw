import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIOnia")

# Conditions
isPbPb = True
isMC   = False
useEventPlane    = False
useGeneralTracks = False
muonSelection  = "Glb" # Single muon selection: Glb(isGlobal), Trk(isTracker) and GlbTrk(isGlobal&&isTracker) are the available options

# Setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = "Jpsi_TEST.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = 'file:/home/llr/cms/stahl/OniaCode/CMSSW_7_5_3_patch1/src/HiAnalysis/HiOnia/test/onia2MuMuPAT_740.root'
#file:/home/llr/cms/chapon/data_CMS/promptskims2015/CMSSW_7_5_4/test/step2_reRECO_740_9_1_e5r.root'
options.maxEvents = 100 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if isMC:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_hi', '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')

#Centrality Tags for CMSSW 7_5_X: 
#            Sample                            Centrality Tag                                             
#        Hydjet 2.76 TeV:  "CentralityTable_HFtowers200_HydjetDrum_run1v750x01_mc"             
#        Hydjet 5 TeV:     "CentralityTable_HFtowers200_HydjetDrum5_v750x02_mc"                   
#        Data 2.76 TeV:    "CentralityTable_HFtowers200_Glauber2010A_eff99_run1v750x01_offline"       
#        Data Run2:        "CentralityTable_HFtowers200_Glauber2015A_v750x01_offline"                    
#Centrality Variables: 
#            Sample             Variable 
#        Hydjet 2.76 TeV:   HFtowersHydjetDrum
#        Hydjet 5 TeV:      HFtowersHydjetDrum5
#        Data 2.76 TeV:     HFtowers
#        Data Run2:         HFtowers
#nonDefaultGlauberModels: 
#            Sample           Name
#        Hydjet 2.76 TeV:   HydjetDrum
#        Hydjet 5 TeV:      HydjetDrum5

if isPbPb:
  process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
  process.centralityBin.Centrality = cms.InputTag("hiCentrality")
  if isMC:
    process.centralityBin.centralityVariable = cms.string("HFtowersHydjetDrum5")
    process.centralityBin.nonDefaultGlauberModel = cms.string("HydjetDrum5")  
    process.GlobalTag.toGet.extend([cms.PSet(
          record  = cms.string("HeavyIonRcd"),
          tag     = cms.string("CentralityTable_HFtowers200_HydjetDrum5_v750x02_mc"),
          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
          label   = cms.untracked.string("HFtowersHydjetDrum5")
          )])
  else:
    process.centralityBin.centralityVariable = cms.string("HFtowers")
    process.centralityBin.nonDefaultGlauberModel = cms.string("")    
    process.GlobalTag.toGet.extend([cms.PSet(
          record  = cms.string("HeavyIonRcd"),
          tag     = cms.string("CentralityTable_HFtowers200_Glauber2010A_eff99_run1v750x01_offline"),
          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
          label   = cms.untracked.string("HFtowers")
          )])


# Event plane (Not working currently)
#process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")

#Options:
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    )
)

'''
#Trigger Filter
process.hltDblMuOpen = cms.EDFilter("HLTHighLevel",
                 TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                 HLTPaths = cms.vstring("HLT_HIL1DoubleMu0_HighQ_v*"),
                 eventSetupPathsKey = cms.string(''),
                 andOr = cms.bool(True),
                 throw = cms.bool(False)
)
'''

process.hionia = cms.EDAnalyzer('HiOniaAnalyzer',
                                #-- Collections
                                srcMuon             = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig       = cms.InputTag("patMuonsWithoutTrigger"),
                                src                 = cms.InputTag("onia2MuMuPatGlbGlb"),
                                srcTracks           = cms.InputTag("hiGeneralTracks"),
                                EvtPlane            = cms.InputTag("hiEvtPlane","recoLevel"),
                                EvtPlaneFlat        = cms.InputTag("hiEvtPlaneFlat",""),

                                triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),

                                #-- Reco Details
                                useBeamSpot = cms.bool(False),
                                useRapidity = cms.bool(True),
                                
                                #--
                                maxAbsZ = cms.double(24.0),
                                
                                pTBinRanges      = cms.vdouble(0.0, 6.0, 8.0, 9.0, 10.0, 12.0, 15.0, 40.0),
                                etaBinRanges     = cms.vdouble(0.0, 2.5),
                                centralityRanges = cms.vdouble(20,40,100),

                                onlyTheBest     = cms.bool(False),		
                                applyCuts       = cms.bool(True),
                                storeEfficiency = cms.bool(False),
                      
                                removeSignalEvents = cms.untracked.bool(False),
                                removeTrueMuons    = cms.untracked.bool(False),
                                storeSameSign      = cms.untracked.bool(True),
                                
                                #-- Gen Details
                                oniaPDG = cms.int32(443),
                                muonSel = cms.string(muonSelection),
                                isHI = cms.untracked.bool(isPbPb),
                                isPA = cms.untracked.bool(False),
                                isMC = cms.untracked.bool(isMC),
                                isPromptMC = cms.untracked.bool(False),
                                useEvtPlane = cms.untracked.bool(useEventPlane),
                                useGeTracks = cms.untracked.bool(useGeneralTracks),
                                runVersionChange = cms.untracked.uint32(182133),

                                #-- Histogram configuration
                                combineCategories = cms.bool(False),
                                fillRooDataSet    = cms.bool(False),
                                fillTree          = cms.bool(True),
                                fillHistos        = cms.bool(False),
                                minimumFlag       = cms.bool(True),
                                fillSingleMuons   = cms.bool(True),
                                fillRecoTracks    = cms.bool(False),
                                histFileName      = cms.string(options.outputFile),		
                                dataSetName       = cms.string(options.secondaryOutputFile),
                                
                                #--
                                # NumberOfTriggers = cms.uint32(8),
                                dblTriggerPathNames = cms.vstring("HLT_HIL1DoubleMu0_HighQ_v2",
                                                                  "HLT_HIL2DoubleMu3_v2",
                                                                  "HLT_HIL3DoubleMuOpen_v2",
                                                                  "HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2"),
                                dblTriggerFilterNames = cms.vstring("hltHIDoubleMuLevel1PathL1HighQFiltered",
                                                                    "hltHIL2DoubleMu3L2Filtered",
                                                                    "hltHIDimuonL3FilteredOpen",
                                                                    "hltHIDimuonL3FilteredMg2OSnoCowboy"),
                                sglTriggerPathNames = cms.vstring("HLT_HIL2Mu3_NHitQ_v2",
                                                                  "HLT_HIL2Mu7_v2",
                                                                  "HLT_HIL2Mu15_v2",
                                                                  "HLT_HIL3Mu3_v2"),
                                sglTriggerFilterNames = cms.vstring("hltHIL2Mu3NHitL2Filtered",
                                                                    "hltHIL2Mu7L2Filtered",
                                                                    "hltHIL2Mu15L2Filtered",
                                                                    "hltHISingleMu3L3Filtered")
                                )

if isPbPb:
  process.hionia.primaryVertexTag = cms.InputTag("hiSelectedVertex")
  process.hionia.genParticles     = cms.InputTag("genParticles")
  process.hionia.muonLessPV       = cms.bool(False)
  process.hionia.CentralitySrc    = cms.InputTag("hiCentrality")
  if isMC:
    process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowersHydjetDrum5")  
  else:
    process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
  
  process.p = cms.Path(process.centralityBin*process.hionia)
else:    
  process.hionia.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")
  process.hionia.genParticles     = cms.InputTag("genParticles")
  process.hionia.muonLessPV       = cms.bool(True)
  process.hionia.CentralitySrc    = cms.InputTag("")
  process.hionia.CentralityBinSrc = cms.InputTag("")

  process.p = cms.Path(process.hionia)

                           

