from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

config = config()
config.section_('General')
date = '2026_06_26'
config.General.workArea = 'crab_projects/'+date+'/MC'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../forest_miniAOD_ParticleTransformer_run3_MC.py'
config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/hintt/Run3_2025_OO/HiForest/'+date+'/MC'
config.Data.publication = False
config.Data.inputDBS = 'global'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*', 'T1_IT_*', 'T1_FR_*', 'T2_US_*', 'T2_FR_*', 'T2_CH_CERN']

dataMap = {}

dataMap["DYToEE_M_10_50_POWHEG"] = { "PD": "/DYto2E_MLL-10to50_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToMuMu_M_10_50_POWHEG"] = { "PD": "/DYto2Mu_MLL-10to50_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToTauTau_M_10_50_POWHEG"] = { "PD": "/DYto2Tau_MLL-10to50_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }

dataMap["DYToEE_M_50_POWHEG"] = { "PD": "/DYto2E_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToMuMu_M_50_POWHEG"] = { "PD": "/DYto2Mu_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToTauTau_M_50_POWHEG"] = { "PD": "/DYto2Tau_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }

dataMap["WmToENu_POWHEG"] = { "PD": "/WminusToEminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v1/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["WpToENu_POWHEG"] = { "PD": "/WplusToEplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["WmToMuNu_POWHEG"] = { "PD": "/WminusToMuminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["WpToMuNu_POWHEG"] = { "PD": "/WplusToMuplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v3/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["WmToTauNu_POWHEG"] = { "PD": "/WminusToTauminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v1/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }
dataMap["WpToTauNu_POWHEG"] = { "PD": "/WplusToTauplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINOOSpring25MiniAOD-150X_mcRun3_2025_forOO_realistic_v9-v1/MINIAODSIM", "Units": 10, "Memory": 3000,   "RunTime": 2749 }

## Submit PDs
for key, val in dataMap.items():
    config.General.requestName = f'HiForest_{key}_OO_5p36TeV_2025Run3_'+date
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.Data.splitting = val['Split'] if "Split" in val else 'LumiBased'
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))
