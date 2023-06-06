### HiForest Configuration
# Collisions: PbPb
# Type: Data
# Input: AOD

import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('HiForest',eras.Run2_2018_pp_on_AA)

###############################################################################
# HiForest labelling info
###############################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest 103X")
import subprocess, os
version = subprocess.check_output(['git',
    '-C', os.path.expandvars('$CMSSW_BASE/src'), 'describe', '--tags'])
if version == '':
    version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

###############################################################################
# Input source
###############################################################################

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #"file:/eos/user/c/cmsdas/2023/long-ex-hin/HIN-HINPbPbAutumn18wmLHEGS-00004.root"
        "file:/afs/cern.ch/user/a/anstahll/work/CMG/CMSDAS/Simulation/CMSSW_10_3_5/src/HIN-HINPbPbAutumn18GS-00055.root"
),
    )

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

###############################################################################
# Load Global Tag, Geometry, etc.
###############################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2018_realistic_HI_v12', '')
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################
# Define tree output
###############################################################################

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestAOD.root"))

###############################################################################
# Additional Reconstruction and Analysis: Main Body
###############################################################################

process.muonTree = cms.EDAnalyzer("MuonTree",
    muons = cms.InputTag(""),
    vertices = cms.InputTag(""),
    genparticles = cms.InputTag("genParticles"),
)

###############################################################################

#########################
# Main analysis list
#########################

process.ana_step = cms.Path(
    process.muonTree
)

###############################################################################

#########################
# Event Selection
#########################

###############################################################################

# Customization
