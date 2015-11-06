import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("HIOnia")

# Conditions
isPbPb = True;
isData = True;
useEventPlane = False;
muonSelection = "GlbTrk" # GlbGlb, GlbTrk, TrkTrk are availale

# Setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = "Jpsi_Histos.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles = 'file:/home/llr/cms/chapon/data_CMS/promptskims2015/CMSSW_7_5_4/test/step2_reRECO_740_9_1_e5r.root'
options.maxEvents = 10 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if isData:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')
else:
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_HIon', '')

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
#        Data Run2:         HFtowers#
#nonDefaultGlauberModels: 
#            Sample           Name
#        Hydjet 2.76 TeV:   HydjetDrum
#        Hydjet 5 TeV:      HydjetDrum5

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi") 
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")   # Only for MC Hydjet  

process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
      tag = cms.string("CentralityTable_HFtowers200_Glauber2010A_eff99_run1v750x01_offline"),
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
      label = cms.untracked.string("HFtowers")
   ),
])

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
                                srcMuon             = cms.InputTag("patMuonsWithTrigger"),
                                srcMuonNoTrig       = cms.InputTag("patMuonsWithoutTrigger"),
                                src                 = cms.InputTag("onia2MuMuPatGlbGlb"),
                                srcTracks           = cms.InputTag("hiGeneralTracks"),
                                genParticles        = cms.InputTag("hiGenParticles"),
                                primaryVertexTag    = cms.InputTag("hiSelectedVertex"),
                                triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
                             
                                CentralitySrc    = cms.InputTag("hiCentrality"),
                                CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
                                EvtPlane         = cms.InputTag("hiEvtPlane","recoLevel"),
                                EvtPlaneFlat     = cms.InputTag("hiEvtPlaneFlat",""),   

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
                                muonLessPV         = cms.bool(False), # this has to be FALSE for HI reco algo
                                
                                #-- Gen Details
                                oniaPDG = cms.int32(443),
                                isHI = cms.untracked.bool(True),
                                isPA = cms.untracked.bool(False),
                                isMC = cms.untracked.bool(False),
                                isPromptMC = cms.untracked.bool(False),
                                useEvtPlane = cms.untracked.bool(useEventPlane),
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

process.p = cms.Path(process.centralityBin*process.hionia)
