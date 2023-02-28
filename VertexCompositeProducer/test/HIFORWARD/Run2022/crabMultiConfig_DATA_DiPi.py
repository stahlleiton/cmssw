from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

date = '20230207'

config.section_("General")
config.General.workArea = 'crab_projects/'+date
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['lumiData.csv']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Collisions2022HISpecial/Cert_Collisions2022HISpecial_362293_362323_Muon.json'
config.Data.runRange = '362293-362323'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {}
for i in range(24):
    dataMap["HITestRaw%s"%i] = { "PD": "/HITestRaw%s/HIRun2022A-PromptReco-v1/AOD"%i, "Units": 25, "Memory": 1700, "RunTime": 1400, "PSet": "PbPbSkimAndTree2022_DiPi_ZDC_ParticleAnalyzer_cfg.py" }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'ParticleTree_DiPi_'+key+'_HIRun2022A_PromptReco_'+date
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = val["PSet"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/CERN/PbPb2022/TREE/%s/%s' % (date, config.General.requestName)
    #config.Data.outLFNDirBase = '/store/user/anstahll/CERN/PbPb2018/TREE/%s' % (config.General.requestName)

    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
