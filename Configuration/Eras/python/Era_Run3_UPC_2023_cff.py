import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_UPC_cff import Run3_UPC
from Configuration.Eras.Modifier_run3_egamma_2023_cff import run3_egamma_2023

Run3_UPC_2023 = cms.ModifierChain(Run3_UPC, run3_egamma_2023)
