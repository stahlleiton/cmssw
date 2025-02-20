from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

config = config()
config.section_('General')
date = '2025_2_10'
config.General.workArea = 'crab_projects/'+date+'/MC'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../forest_miniAOD_ParticleTransformer_run3_MC.py'
config.JobType.numCores = 1
config.section_('Data')
config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/group/cmst3/group/hintt/Run3/HiForest/'+date
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*', 'T2_US_Caltech', 'T2_US_Purdue', 'T2_US_Vanderbilt', 'T1_FR_*', 'T2_FR_*', 'T2_CH_CERN']

dataMap = {}
dataMap["TT_hvq_POWHEG_NONEMB_Official"           ] = { "PD": "/TT_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-NoPU_132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["TT_hvq_POWHEG_Hydjet_Official"           ] = { "PD": "/TT_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["TT012J_NLO_FXFX_MADGRAPH_Hydjet_Official"] = { "PD": "/TT-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }

dataMap["TWminus_POWHEG_Hydjet_Official"   ] = { "PD": "/TWminus_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["TWplus_POWHEG_Hydjet_Official"    ] = { "PD": "/TWplus_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["SingleT_POWHEG_Hydjet_Official"   ] = { "PD": "/singleT_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["SingleTbar_POWHEG_Hydjet_Official"] = { "PD": "/singleTbar_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }

dataMap["DYToEE_M_50_POWHEG_Hydjet_Official"                  ] = { "PD": "/DYto2E_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["DYToMuMu_M_50_POWHEG_Hydjet_Official"                ] = { "PD": "/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["DYToTauTau_M_50_POWHEG_Hydjet_Official"              ] = { "PD": "/DYto2Tau_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["DY012JToLL_M_50_NLO_FXFX_MADGRAPH_Hydjet_Official"   ] = { "PD": "/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["DY012JToLL_M_10_50_NLO_FXFX_MADGRAPH_Hydjet_Official"] = { "PD": "/DYto2L-2Jets_MLL-10to50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }

dataMap["WpToE_POWHEG_Hydjet_Official"                ] = { "PD": "/WplusToEplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime":  2749 }
dataMap["WpToMu_POWHEG_Hydjet_Official"               ] = { "PD": "/WplusToMuplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime":  2749 }
dataMap["WpToTau_POWHEG_Hydjet_Official"              ] = { "PD": "/WplusToTauplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime":  2749 }
dataMap["WmToE_POWHEG_Hydjet_Official"                ] = { "PD": "/WminusToEminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime":  2749 }
dataMap["WmToMu_POWHEG_Hydjet_Official"               ] = { "PD": "/WminusToMuminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime":  2749 }
dataMap["WmToTau_POWHEG_Hydjet_Official"              ] = { "PD": "/WminusToTauminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime":  2749 }
dataMap["WWto2L2Nu_POWHEG_Hydjet_Official"            ] = { "PD": "/WWto2L2Nu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["WWtoLNu2Q_POWHEG_Hydjet_Official"            ] = { "PD": "/WWtoLNu2Q_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM", "Units": 10, "Memory": 3900, "RunTime": 2749 }
dataMap["W012JToLNu_NLO_FXFX_MADGRAPH_Hydjet_Official"] = { "PD": "/WtoLNu-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3900,    "RunTime": 2749 }

dataMap["QCDToMu_PYTHIA8_Hydjet_Official"] = { "PD": "/QCDtoMuons_Pthat20_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v1/MINIAODSIM", "Units": 100, "Memory": 3900, "RunTime": 2749 }
dataMap["QCDToE_PYTHIA8_Hydjet_Official" ] = { "PD": "/QCDtoElectrons_Pthat20_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v1/MINIAODSIM", "Units": 20, "Memory": 3900, "RunTime": 2749 }

dataMap["TT_hvq_POWHEG_HERWIG_NONEMB" ] = { "PD": "/POWHEG_5p36TeV_2023Run3/cmst3-TT_hvq_POWHEG_HERWIG_NONEMB_5p36TeV_TuneCP5_2023Run3_MINIAOD_2024_04_19-b75d3b5331198c44808663383a6aa555/USER", "Units": 100, "Memory": 3900, "RunTime": 2749 }
dataMap["DiJet_pTHat15_PYTHIA8_Hydjet"] = { "PD": "/PYTHIA8_2023Run3/cmst3-DiJet_pTHat15_PYTHIA8_Hydjet_5p36TeV_TuneCP5_2023Run3_MINIAOD_2024_04_19-b75d3b5331198c44808663383a6aa555/USER", "Units": 20, "Memory": 3900, "RunTime": 1440 }
dataMap["BJet_pTHat15_PYTHIA8_Hydjet" ] = { "PD": "/PYTHIA8_2023Run3/cmst3-BJet_pTHat15_PYTHIA8_Hydjet_5p36TeV_TuneCP5_2023Run3_MINIAOD_2024_04_19-b75d3b5331198c44808663383a6aa555/USER", "Units": 20, "Memory": 3900, "RunTime": 1440 }
dataMap["BJet_pTHat80_PYTHIA8_Hydjet" ] = { "PD": "/PYTHIA8_2023Run3/cmst3-BJet_pTHat80_PYTHIA8_Hydjet_5p36TeV_TuneCP5_2023Run3_MINIAOD_2024_04_19-b75d3b5331198c44808663383a6aa555/USER", "Units": 20, "Memory": 3900, "RunTime": 1440 }
dataMap["Z012JToNuNu_Zpt100_NLO_FXFX_MADGRAPH_Hydjet"] = { "PD": "/MADGRAPH_5p36TeV_2023Run3/cmst3-Z012JToNuNu_Zpt100_NLO_FXFX_MADGRAPH_Hydjet_5p36TeV_TuneCP5_2023Run3_MINIAOD_2024_04_19-b75d3b5331198c44808663383a6aa555/USER", "Units": 12, "Memory": 3900, "RunTime": 1440 }

## Submit PDs
for key, val in dataMap.items():
    config.General.requestName = f'HiForest_{key}_5p36TeV_TuneCP5_2023Run3_'+date
    config.Data.inputDataset = val["PD"]
    config.Data.inputDBS = 'global' if ("HINPbPbSpring23MiniAOD" in val["PD"]) else 'phys03'
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))
