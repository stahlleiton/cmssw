from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.section_('General')
config.General.requestName = 'QCD_5p5TeV_TuneCP5_MTD_ANA_20190306'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'TrackAnalysis/test/test.py'
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 1400
config.section_('Data')
config.Data.inputDataset = '/QCD_Pt_120_170_5p5TeV_TuneCP5_Pythia8/PhaseIIMTDTDRAutumn18DR-NoPU_103X_upgrade2023_realistic_v2-v1/FEVT'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/MTD/MC/Official/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_FR_*' , 'T2_CH_*' , 'T2_US_*' , 'T2_FR_*']
config.Site.storageSite = 'T2_CH_CERN'
