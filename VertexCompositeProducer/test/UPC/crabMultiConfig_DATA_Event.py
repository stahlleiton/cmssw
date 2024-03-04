from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

date = '2024_02_04'
config.section_("General")
config.General.workArea = 'crab_projects/Event/'+date
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.scriptExe = 'submitScript.sh'
config.JobType.inputFiles = ['emap_2023_newZDC_v3.txt', 'CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db']

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23HI/Cert_Collisions2023HI_374288_375823_Golden.json'
#config.Data.runRange = '374288-375823'
config.Data.publication = False

config.section_('Site')
config.Site.whitelist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T2_CH_CERN'

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
inputDir = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/SKIM"
dataMap["HIForward_TRK_UPC"] = { "PD": inputDir+f"/hiforward_skim_track_upcreco.txt" ,  "Units": 10, "Memory": 2000, "RunTime": 350 }
dataMap["HIForward_EGM_UPC"] = { "PD": inputDir+"/hiforward_skim_egamma_upcreco.txt" , "Units": 10, "Memory": 2000, "RunTime": 350 }
dataMap["HIForward_TRK_PP"] = { "PD": inputDir+f"/hiforward_skim_track_ppreco.txt" ,  "Units": 10, "Memory": 2000, "RunTime": 350 }
dataMap["HIForward_EGM_PP"] = { "PD": inputDir+"/hiforward_skim_egamma_ppreco.txt" , "Units": 10, "Memory": 2000, "RunTime": 350 }
dataMap["HIEmptyBX_UPC"] =   { "PD": inputDir+"/hiemptybx_upcreco.txt",              "Units": 1,  "Memory": 2000, "RunTime": 350 }
dataMap["HIEmptyBX_PP"] =    { "PD": inputDir+"/hiemptybx_ppreco.txt",               "Units": 1,  "Memory": 2000, "RunTime": 350 }
dataMap["HIZeroBias_UPC"] =  { "PD": inputDir+"/hizerobias_upcreco.txt",              "Units": 1,  "Memory": 2000, "RunTime": 350 }
dataMap["HIZeroBias_PP"] =   { "PD": inputDir+"/hizerobias_ppreco.txt",               "Units": 1,  "Memory": 2000, "RunTime": 350 }

## Submit the PDs
for key, val in dataMap.items():
    config.General.requestName = 'ParticleAnalyzer_Event_'+key+'_HIRun2023A_'+date
    config.Data.userInputFiles = open(val["PD"]).readlines()
    config.Data.totalUnits = len(config.Data.userInputFiles)
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = "ParticleAnalyzer_EventInfo_cfg.py"
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/Event/'+date

    print("Submitting CRAB job for: "+val["PD"])
    submit(config)
