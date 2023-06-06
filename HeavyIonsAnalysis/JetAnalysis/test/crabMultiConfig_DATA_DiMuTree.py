from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
config.Data.runRange = '326381-327564'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = False

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=True)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {
            "HIForward": { "PD": "/HIForward/HIRun2018A-04Apr2019-v1/AOD", "Units": 45, "Memory": 2400, "RunTime": 1400, "PSet": "runForestAOD_CMSDAS_DATA_103X.py" },
            }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'HiForest_'+key+'_HIRun2018_04Apr2019_CMSDAS_20230606'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = val["PSet"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/CMSDAS/HiForest/%s' % (config.General.requestName)

    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
