from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

config = config()
config.section_('General')
date = '2026_06_26'
config.General.workArea = 'crab_projects/'+date+'/DATA/PHOLEP'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../forest_miniAOD_ParticleTransformer_run3_DATA.py'
config.JobType.inputFiles = ['../phoEleReg_Run3_2025_PbPb.db']
config.JobType.maxMemoryMB = 2999
config.JobType.maxJobRuntimeMin = 1749
config.section_('Data')
config.Data.outLFNDirBase = '/store/group/cmst3/group/hintt/Run3_2025_PbPb/HiForest/'+date+'/DATA/PHOLEP'
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True # to fix issue of missing blocks
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*', 'T1_FR_*', 'T2_US_MIT', 'T2_FR_*', 'T2_ES_*', 'T2_UK_*', 'T2_US_Vanderbilt', 'T2_CH_CERN']
config.Site.blacklist = ['T2_CN_*', 'T2_TW_*', 'T2_DE_*', 'T2_EE_*']

def submit(config, dryrun):
    try:
        crabCommand('submit', config = config, dryrun=dryrun)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 25
config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions25HI/Cert_Collisions2025_HI_399465_400426_Golden.json'
config.Data.inputDBS = 'global'
## Submit the muon PDs
for i in range(0, 60, 1):
    config.General.requestName = f'HiForest_HIPhysicsRawPrime{i}_HIRun2025A_PromptReco_v1_PHOLEPSKIM_'+date
    config.Data.inputDataset = f'/HIPhysicsRawPrime{i}/HIRun2025A-PromptReco-v1/MINIAOD'
    config.Data.outputDatasetTag = config.General.requestName
    submit(config = config, dryrun=False)
