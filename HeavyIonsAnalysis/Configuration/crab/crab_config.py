#Crab config for RECO AOD
from CRABClient.UserUtilities import config

config = config()

config.section_('General')
config.General.workArea = 'MC_projects'
config.General.requestName ='HiForest_HYDJET'
config.General.transferOutputs = True
config.General.transferLogs = True


config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyIonsAnalysis/Configuration/test/forest_miniAOD_run3_MC_HYDJET.py'
config.JobType.allowUndistributedCMSSW =True

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.inputDataset = '/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbAutumn22DR-NoPU_125X_mcRun3_2022_realistic_HI_v13-v2/AODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = 'HiForest_HYDJET'

config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/GO2023/HiForest_HYDJET'
config.Site.storageSite = 'T2_CH_CERN'
