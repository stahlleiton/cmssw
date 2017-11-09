from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'HIEWQ_pPb8TeV_AOD_20170913'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'skim_pPb.py'
config.JobType.maxMemoryMB = 2400

config.section_('Data')
config.Data.inputDataset ='/PASingleMuon/PARun2016C-PromptReco-v1/AOD'  
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 10
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/home/llr/cms/stahl/ElectroWeakAnalysis/CMSSW_8_0_29/src/HeavyIonsAnalysis/ElectroWeakAnalysis/test/DATA/Cert_285479-286496_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
config.Data.runRange = '285479-286496'
config.Data.outLFNDirBase = '/store/user/%s/EWQAnalysis2017/SKIM/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName

config.section_('Site')
config.Site.storageSite = 'T2_FR_GRIF_LLR'
