from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

config = config()
config.section_('General')
date = '2026_07_13'
config.General.workArea = 'crab_projects/'+date+'/MC'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../forest_miniAOD_ParticleTransformer_run3_MC.py'
config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/hintt/Run3_2024_ppRef/HiForest/'+date+'/MC'
config.Data.publication = False
config.Data.inputDBS = 'global'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*', 'T1_IT_*', 'T1_FR_*', 'T2_US_*', 'T2_IT_*', 'T2_FR_*', 'T2_DE_*', 'T2_CH_*']

dataMap = {}

dataMap["TT_hvq_POWHEG_HERWIG_NOPU"] = { "PD": "/TTbar_TuneCH3_5p36TeV_powheg-herwig7/RunIIIpp5p36Winter24MiniAOD-NoPU_141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
dataMap["TT_hvq_POWHEG_NOPU"       ] = { "PD": "/TT_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-NoPU_141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
dataMap["TT_hvq_POWHEG"            ] = { "PD": "/TT_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
dataMap["TT012J_NLO_FXFX_MADGRAPH" ] = { "PD": "/TT-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }

#dataMap["TWminus_POWHEG"   ] = { "PD": "/T_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["TWplus_POWHEG"    ] = { "PD": "/Tbar_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["SingleT_POWHEG"   ] = { "PD": "/T-tChannel_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["SingleTbar_POWHEG"] = { "PD": "/Tbar-tChannel_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }

#dataMap["DYToEE_M_50_POWHEG"                  ] = { "PD": "/DYToEE_M-50_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["DYToMuMu_M_50_POWHEG"                ] = { "PD": "/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["DYToTauTau_M_50_POWHEG"              ] = { "PD": "/DYToTauTau_M-50_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["DYToEE_M_10_50_POWHEG"               ] = { "PD": "/DYToEE_M-10to50_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["DYToMuMu_M_10_50_POWHEG"             ] = { "PD": "/DYToMuMu_M-10to50_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["DYToTauTau_M_10_50_POWHEG"           ] = { "PD": "/DYToTauTau_M-10to50_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
dataMap["DY012JToLL_M_50_NLO_FXFX_MADGRAPH"   ] = { "PD": "/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["DY012JToLL_M_10_50_NLO_FXFX_MADGRAPH"] = { "PD": "/DYto2L-2Jets_MLL-10to50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }

#dataMap["WpToE_POWHEG"                ] = { "PD": "/WplusToEplusNu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WpToMu_POWHEG"               ] = { "PD": "/WplusToMuplusNu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WpToTau_POWHEG"              ] = { "PD": "/WplusToTauplusNu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WmToE_POWHEG"                ] = { "PD": "/WminusToEminusNu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WmToMu_POWHEG"               ] = { "PD": "/WminusToMuminusNu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WmToTau_POWHEG"              ] = { "PD": "/WminusToTauminusNu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WWto2L2Nu_POWHEG"            ] = { "PD": "/WWTo2L2Nu_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v1/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["WWtoLNu2Q_POWHEG"            ] = { "PD": "/WWToLNu2Q_TuneCP5_5p36TeV_powheg-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
dataMap["W012JToLNu_NLO_FXFX_MADGRAPH"] = { "PD": "/WtoLNu-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["W01234JToLNu_LO_MLM_MADGRAPH"] = { "PD": "/WtoLNu-4Jets_TuneCP5_5p36TeV_madgraphMLM-pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
#dataMap["W01234JToLNu_4J_LO_MLM_MADGRAPH"] = { "PD": "", "Units": 1, "Split": "FileBased", "Memory": 2999, "RunTime": 720 }

dataMap["QCDToMu_PYTHIA8"] = { "PD": "/QCD-Mu_pThat-20_TuneCP5_5p36TeV_pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }
dataMap["QCDToE_PYTHIA8" ] = { "PD": "/QCD-E_pThat-20_TuneCP5_5p36TeV_pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 720 }

##dataMap["DiJet_pTHat15_PYTHIA8"] = { "PD": "/QCD_pThat-15to1200_TuneCP5_5p36TeV_pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 5, "Memory": 2999, "RunTime": 720, "MaxUnits": 10000 }
#dataMap["BJet_pTHat15_PYTHIA8" ] = { "PD": "/QCD_BEnriched_pThat-15to500_TuneCP5_5p36TeV_pythia8/RunIIIpp5p36Winter24MiniAOD-141X_mcRun3_2024_realistic_ppRef5TeV_v7-v2/MINIAODSIM", "Units": 5, "Memory": 2999, "RunTime": 720, "MaxUnits": 10000 }

## Submit PDs
for key, val in dataMap.items():
    config.General.requestName = f'HiForest_{key}_ppRef_5p36TeV_2024Run3_'+date
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.Data.splitting = val['Split'] if "Split" in val else 'LumiBased'
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.allowNonValidInputDataset = val["PRODUCTION"] if "PRODUCTION" in val else False
    config.Data.totalUnits = val["MaxUnits"] if "MaxUnits" in val else 10000000000
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))
