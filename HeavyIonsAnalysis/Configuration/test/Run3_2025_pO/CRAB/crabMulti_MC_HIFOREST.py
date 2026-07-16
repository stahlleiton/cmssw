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
config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/CERN/pO2025/HiForest/'+date+'/MC'
config.Data.publication = False
config.Data.splitting = 'FileBased'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

def submit(config, dryrun):
    try:
        crabCommand('submit', config = config, dryrun=dryrun)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

dataMap = {}
dataMap["DYToEE_M_10_50_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-DYToEE_M_10_50_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToMuMu_M_10_50_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-DYToMuMu_M_10_50_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToTauTau_M_10_50_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-DYToTauTau_M_10_50_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000,   "RunTime": 2749 }

dataMap["DYToEE_M_50_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-DYToEE_M_50_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }
dataMap["DYToMuMu_M_50_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-DYToMuMu_M_50_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000,   "RunTime": 2749 }
dataMap["DYToTauTau_M_50_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-DYToTauTau_M_50_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000,   "RunTime": 2749 }

dataMap["WmToENu_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-WmToENu_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }
dataMap["WpToENu_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-WpToENu_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }
dataMap["WmToMuNu_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-WmToMuNu_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }
dataMap["WpToMuNu_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-WpToMuNu_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }
dataMap["WmToTauNu_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-WmToTauNu_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }
dataMap["WpToTauNu_POWHEG"] = { "PD": "/POWHEG_9p62TeV_2025Run3/phys_heavyions-WpToTauNu_POWHEG_pO_9p62TeV_TuneCP5_2025Run3_RECO_2025_10_10-a387798f8dc4c0c5c1970bf99cec353f/USER", "Units": 6, "Memory": 3000, "RunTime": 2749 }

## Submit PDs
for key, val in dataMap.items():
    config.General.requestName = f'HiForest_{key}_pO_9p62TeV_2025Run3_'+date
    config.Data.inputDataset = val["PD"]
    config.Data.inputDBS = 'global' if ("HINPbPbSpring23MiniAOD" in val["PD"]) else 'phys03'
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    submit(config = config, dryrun=False)
