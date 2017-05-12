from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'PARun2016C-v1_Run_285479_286504_EWQANA_pA_20170430'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'lo.py'
config.JobType.maxMemoryMB = 2400

config.section_('Data')
config.Data.inputDataset ='/PASingleMuon/anstahll-HIEWQ_pPb8TeV_AOD_20170207-ea47420fa7d4dfad582bb824ce988a1b/USER'
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/EWQAnalysis2017/Tree/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T2_FR_GRIF_*']
config.Site.storageSite = 'T2_FR_GRIF_LLR'
