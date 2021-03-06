import FWCore.ParameterSet.Config as cms

import Geometry.HcalEventSetup.hcalTrigTowerGeometry_cfi

hcalTrigTowerGeometry = Geometry.HcalEventSetup.hcalTrigTowerGeometry_cfi.hcalTrigTowerGeometry.clone()

import Geometry.HcalEventSetup.hcalTopologyConstants_cfi as hcalTopologyConstants_cfi
hcalTrigTowerGeometry.hcalTopologyConstants = cms.PSet(hcalTopologyConstants_cfi.hcalTopologyConstants)
