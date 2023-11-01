import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.ProcessModifiers.egamma_lowPt_exclusive_cff import egamma_lowPt_exclusive
from Configuration.Eras.Modifier_highBetaStar_2023_cff import highBetaStar_2023
from Configuration.Eras.Modifier_run3_upc_2023_cff import run3_upc_2023

Run3_UPC = cms.ModifierChain(Run3, egamma_lowPt_exclusive, highBetaStar_2023, run3_upc_2023)
