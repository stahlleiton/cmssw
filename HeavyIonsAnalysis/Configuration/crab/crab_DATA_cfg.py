#Crab config for RECO AOD
from CRABClient.UserUtilities import config

config = config()

config.section_('General')
config.General.workArea = 'crab_projects_new'
config.General.requestName ='HiForest_TestRun_data_ZB1_GT_20230811'
config.Data.outputDatasetTag = 'HiForest_TestRun_data_ZB1_GT'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyIonsAnalysis/Configuration/test/forest_AODPAT_run3_DATA.py'
config.JobType.allowUndistributedCMSSW =True
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDataset = '/HITestRaw1/HIRun2022A-PromptReco-v1/AOD'
config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Collisions2022HISpecial/Cert_Collisions2022HISpecial_362293_362323_Golden.json'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 4

#config.Data.allowNonValidInputDataset = True

config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/GO2023/HiForest_TestRun_data_originalGT'
config.Site.storageSite = 'T2_CH_CERN'
