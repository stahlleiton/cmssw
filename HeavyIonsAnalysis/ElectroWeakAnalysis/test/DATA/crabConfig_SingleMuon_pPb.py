from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'PARun2016C-v1_Run_285410_285951_EWQANA_pPb_20170204'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'lo.py'
config.JobType.maxMemoryMB = 2400

config.section_('Data')
config.Data.inputDataset ='/PASingleMuon/PARun2016C-PromptReco-v1/AOD'  
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 10
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
config.Data.runRange = '285410-285951'
config.Data.outLFNDirBase = '/store/user/%s/EWQAnalysis2017/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False
#config.Data.outputDatasetTag = config.General.requestName

config.section_('Site')
config.Site.storageSite = 'T2_FR_GRIF_LLR'
