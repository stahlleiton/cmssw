#!/bin/bash

# Create CMSSW working area
#cmsrel CMSSW_4_4_7
cmsenv
git cms-init

# Add remote repository
git remote add CMS-HIN-dilepton https://github.com/CMS-HIN-dilepton/cmssw.git
git remote update
git fetch CMS-HIN-dilepton

# Bring a working branch from remote repository
git checkout -b CMSSW_4_4_X_Lxyz CMS-HIN-dilepton/CMSSW_4_4_X_Lxyz

# Copy CMSSW modules from remote repository
git cms-addpkg CmsHi/Analysis2010
git cms-addpkg CondFormats/HIObjects
git cms-addpkg DataFormats/HeavyIonEvent
git cms-addpkg HeavyIonsAnalysis/Configuration
git cms-addpkg RecoHI/HiCentralityAlgos
git cms-addpkg RecoHI/HiEvtPlaneAlgos
git cms-addpkg RecoLuminosity/LumiDB
git cms-addpkg RecoVertex/PrimaryVertexProducer
git cms-addpkg MuonAnalysis/MuonAssociators
git cms-addpkg MuonAnalysis/TagAndProbe
git cms-addpkg PhysicsTools/TagAndProbe
git cms-addpkg PhysicsTools/UtilAlgos
git cms-addpkg HiSkim/HiOnia2MuMu
git cms-addpkg HiAnalysis/HiOnia

scram b
