from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
config = config()
config.section_('General')
tryMulti = True
date = '2025_2_10'
config.General.workArea = 'crab_projects/'+date+'/DATA'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../forest_miniAOD_ParticleTransformer_run3_SKIM_DATA.py'
config.JobType.maxMemoryMB = 3500
config.JobType.maxJobRuntimeMin = 2749
config.section_('Data')
config.Data.outLFNDirBase = '/store/group/cmst3/group/hintt/Run3/HiForest/'+date+'/DATA'
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True # to fix issue of missing blocks
#config.Site.blacklist = ['T2_EE_Estonia']

if tryMulti:
    config.Data.splitting = 'LumiBased'
    config.Data.unitsPerJob = 100
    config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23HI/Cert_Collisions2023HI_374288_375823_Golden.json'
    config.Data.inputDBS = 'global'
    ## Submit the muon PDs
    for i in range(0, 32, 1):
        config.General.requestName = f'HiForest_HIPhysicsRawPrime{i}_HIRun2023A_PromptReco_v2_SKIM_'+date
        config.Data.inputDataset = f'/HIPhysicsRawPrime{i}/HIRun2023A-PromptReco-v2/MINIAOD'
        config.Data.outputDatasetTag = config.General.requestName
        crabCommand('submit', config = config, dryrun=False)
else:
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 40
    infile = 'minbias_prompt_v2.txt'
    config.Data.userInputFiles = open(infile).readlines()
    config.Data.totalUnits = len(config.Data.userInputFiles)
    config.Data.outputPrimaryDataset = 'HIPhysicsRawPrime'
    config.Site.whitelist = ['T2_US_Vanderbilt']
    config.General.requestName = 'HiForest_HIPhysicsRawPrime_HIRun2023A_PromptReco_v2_SKIM_'+date
    config.Data.outputDatasetTag = config.General.requestName
    crabCommand('submit', config = config, dryrun=False)
