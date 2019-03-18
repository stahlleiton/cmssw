
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hioniaanalyzer_PbPbPrompt_103X_DATA_vnmerge_cfg.py'
config.JobType.maxMemoryMB = 2900
config.JobType.maxJobRuntimeMin = 600
config.JobType.inputFiles = ['HeavyIonRPRcd_PbPb2018_offline.db']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 50
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
config.Data.runRange = '326381-327564'
config.Data.publication = False

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

## Submit the v1 HIMinimumBias samples
for x in range(0,20):
    config.General.requestName = 'OniaVNTrees_HIMinimumBias'+str(x)+'_v1_HIRun2018_Upsilon_20190318'
    config.Data.inputDataset = '/HIMinimumBias'+str(x)+'/HIRun2018A-PromptReco-v1/AOD'
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/UpsilonV2/PbPb2018/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)

## Submit the v2 HIMinimumBias samples
for x in range(0,20):
    config.General.requestName = 'OniaVNTrees_HIMinimumBias'+str(x)+'_v2_HIRun2018_Upsilon_20190318'
    config.Data.inputDataset = '/HIMinimumBias'+str(x)+'/HIRun2018A-PromptReco-v2/AOD'
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/UpsilonV2/PbPb2018/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)
