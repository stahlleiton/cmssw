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
config.Data.splitting = 'FileBased'
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
            "HIForward": { "PD": "/JPsiToMuMu_5p02TeV_PbPb_Coherent_STARlight/HINPbPbAutumn18GS-103X_upgrade2018_realistic_HI_v11-v3/GEN-SIM", "Units": 1, "Memory": 1800, "RunTime": 1400, "PSet": "runForestAOD_CMSDAS_MC_GEN_103X.py" },
          }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'HiForest_'+key+'_STARlight_CohJpsi2MuMu_PbPb5TeV_GENONLY_CMSDAS_20230526'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = val["PSet"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/CMSDAS/HiForest/%s' % (config.General.requestName)

    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
