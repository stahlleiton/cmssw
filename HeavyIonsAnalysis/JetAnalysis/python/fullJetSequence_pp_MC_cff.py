import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.rerecoJets_cff import *

from HeavyIonsAnalysis.JetAnalysis.jets.ak3CaloJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pp_mc_cff import *

from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pp_mc_cff import *

jetSequences = cms.Sequence(
    ak3CaloJets +
    ak3PFJets +

    ak4PFJets +
    ak4CaloJets +

    ak3CaloJetSequence +
    ak3PFJetSequence +

    ak4CaloJetSequence +
    ak4PFJetSequence
)
