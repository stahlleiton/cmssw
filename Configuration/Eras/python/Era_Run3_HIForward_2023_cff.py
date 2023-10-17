import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_HIForward_cff import Run3_HIForward
from Configuration.Eras.Modifier_run3_egamma_2023_cff import run3_egamma_2023

Run3_HIForward_2023 = cms.ModifierChain(Run3_HIForward, run3_egamma_2023)
