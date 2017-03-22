// -*- C++ -*-
//
// Package:    HiMETAnalyzer
// Class:      HiMETAnalyzer
// 
/**\class HiMETAnalyzer HiHIMETAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Andre, Stahl 
//         Created:  Jan  25 2017
// 
//
//

#include "HeavyIonsAnalysis/ElectroWeakAnalysis/interface/HiMETAnalyzer.h"

//Headers for the core items
#include "FWCore/Framework/interface/MakerMacros.h"

//system include files
#include <iostream>


//
// constructors and destructor
//
HiMETAnalyzer::HiMETAnalyzer(const edm::ParameterSet& iConfig):
  _patMETToken ( consumes< edm::View< pat::MET      > >( iConfig.getParameter< edm::InputTag >("patMETTag" ) )),
  _pfMETToken  ( consumes< edm::View< reco::PFMET   > >( iConfig.getParameter< edm::InputTag >("pfMETTag"  ) )),
  _caloMETToken( consumes< edm::View< reco::CaloMET > >( iConfig.getParameter< edm::InputTag >("caloMETTag") )),
  _genMETToken ( consumes< edm::View< reco::GenMET  > >( iConfig.getParameter< edm::InputTag >("genMETTag" ) )),
  _filterResultToken( consumes< edm::TriggerResults >( edm::InputTag("TriggerResults") )),
  _metFilters( iConfig.getParameter< std::vector< std::string > >( "metFilters") ),
  _corrNames( iConfig.getParameter< std::vector< std::string > >("corrections") ),
  _doAll( iConfig.getParameter< bool >( "doAll" ) ),
  _eventFilter( iConfig.getParameter< std::string >( "eventFilter" ) )
{
}

HiMETAnalyzer::~HiMETAnalyzer()
{
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
HiMETAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::map< std::string , bool > doMET;

  // Fill Event Tree always
  doMET["Event"] = true;

  // Extract the results of the MET Filters
  std::map< std::string , Bool_t >  metFilter;
  edm::Handle< edm::TriggerResults > filterResult;
  iEvent.getByToken( _filterResultToken, filterResult );
  if (filterResult.isValid()) {
    edm::TriggerNames const& filterNames = iEvent.triggerNames(*filterResult);
    // Check if we want to keep the event
    uint evtFilterIndex = filterNames.triggerIndex(_eventFilter);
    if ( (evtFilterIndex >= filterNames.size()) || !filterResult->accept(evtFilterIndex)) return;
    // Continue to check the MET filters
    for ( size_t iFilter = 0; iFilter < _metFilters.size(); ++iFilter ) {
      std::string filterName = _metFilters.at(iFilter);
      if (filterName.find("Flag_")==std::string::npos) { filterName = std::string("Flag_") + filterName; }
      uint index = filterNames.triggerIndex(filterName);
      if ( index < filterNames.size() ) {
        metFilter[_metFilters.at(iFilter)] = filterResult->accept(index);
      }
    }
  }
  if (metFilter.size()>0) { doMET["Filter"] = true; }
  else { 
    edm::LogWarning("METAnalyzer_FilterNotFound") << "No filters were found!" << std::endl;
    doMET["Filter"] = false; 
  }

  // If Gen MET tree exist, fill it
  edm::Handle<  edm::View<reco::GenMET> > genMETHandle;
  getCollection(iEvent, _genMETToken, genMETHandle);
  if (genMETHandle.isValid() && genMETHandle->size()!=1) {
    edm::LogError("METAnalyzer_GenMETCollection_corrupt") << "Number of GEN METs in event is different from 1!" << std::endl;
    return;
  }
  if (genMETHandle.isValid() && genMETHandle->size()==1) { doMET["Gen"] = true; }
  else  { doMET["Gen"] = false; }

  // If Calo MET tree exist, fill it
  edm::Handle<  edm::View<reco::CaloMET> > caloMETHandle;
  getCollection(iEvent, _caloMETToken, caloMETHandle);
  if (caloMETHandle.isValid() && caloMETHandle->size()!=1) {
    edm::LogError("METAnalyzer_CaloMETCollection_corrupt") << "Number of Calo METs in event is different from 1!" << std::endl;
    return;
  }
  if (caloMETHandle.isValid() && caloMETHandle->size()==1) { doMET["Calo"] = true; }
  else  { doMET["Calo"] = false; }

  // If PF MET tree exist, fill it
  edm::Handle<  edm::View<reco::PFMET> > pfMETHandle;
  getCollection(iEvent, _pfMETToken, pfMETHandle);
  if (pfMETHandle.isValid() && pfMETHandle->size()!=1) {
    edm::LogError("METAnalyzer_PFMETCollection_corrupt") << "Number of PF METs in event is different from 1!" << std::endl;
    return;
  }
  if (pfMETHandle.isValid() && pfMETHandle->size()==1) { doMET["PF"] = true; }
  else  { doMET["PF"] = false; }

  // If PAT MET tree exist, fill it
  edm::Handle<  edm::View<pat::MET> > patMETHandle;
  getCollection(iEvent, _patMETToken, patMETHandle);
  if (patMETHandle.isValid() && patMETHandle->size()!=1) {
    edm::LogError("METAnalyzer_PATMETCollection_corrupt") << "Number of PAT METs in event is different from 1!" << std::endl;
    return;
  }
  if (patMETHandle.isValid() && patMETHandle->size()==1) { doMET["Pat"] = true; }
  else  { doMET["Pat"] = false; }

  // Check if at least one of the collections was found
  if (!doMET["Pat"] && !doMET["Calo"] && !doMET["PF"] && !doMET["Gen"]) { 
    edm::LogError("METAnalyzer_ProductNotFound") << "No collections were found!" << std::endl;
    return;
  }

  // Fill Reco only if there is at least one reconstructed collection
  if (doMET["Pat"] || doMET["PF"] || doMET["Calo"]) { doMET["Reco"] = true; }
  else { doMET["Reco"] = false; }
  
  std::vector< std::string > treeNames;
  if (_doAll) { treeNames.push_back("All"); }
  else {
    if (doMET["Event"]) treeNames.push_back("Event");
    if (doMET["Reco"])  treeNames.push_back("Reco");
    if (doMET["PF"])    treeNames.push_back("PF");
    if (doMET["Calo"])  treeNames.push_back("Calo");
    if (doMET["Pat"])   for (uint i = 0; i < _corrNames.size(); i++) { treeNames.push_back(_corrNames[i]); }
    if (doMET["Gen"])   treeNames.push_back("Gen");
    if (doMET["Filter"]) treeNames.push_back("Filter");
  }

  for (uint i = 0; i < treeNames.size(); i++) {
    std::string name = treeNames[i];
    metEvt_[name].Clear();
    if (doMET["Event"] && (name=="Event" || name=="All")) { metEvt_[name].Fill( iEvent ); }
    if (doMET["Reco"]  && (name=="Reco" || name=="All")) { 
      if (doMET["Pat"])        { metEvt_[name].FillReco( patMETHandle->front()  ); }
      else if (doMET["PF"])    { metEvt_[name].FillReco( pfMETHandle->front()   ); }
      else if (doMET["Calo"])  { metEvt_[name].FillReco( caloMETHandle->front() ); }
    }
    if (doMET["PF"]    && (name=="PF" || name=="All")) { metEvt_[name].Fill( pfMETHandle->front() ); }
    if (doMET["Calo"]  && (name=="Calo" || name=="All")) { metEvt_[name].Fill( caloMETHandle->front() ); }
    if (doMET["Pat"]) {
      if (name=="All" || ( name!="Event" && name!="Reco" && name!="Gen" && name!="Filter" )) { metEvt_[name].Fill( patMETHandle->front() , name ); }
    }
    if (doMET["Gen"]     && (name=="Gen" || name=="All"))    { metEvt_[name].Fill( genMETHandle->front() ); }
    if (doMET["Filter"]  && (name=="Filter" || name=="All")) { metEvt_[name].Fill( metFilter ); }
    if (firstEvent_) {
      metTree_[name] = fs_->make<TTree>(Form("MET_%s", name.c_str()), "");
      metEvt_[name].SetTree( metTree_[name] );
      metEvt_[name].SetBranches( name , doMET);
    }
  }
  if (firstEvent_) { firstEvent_ = false; }

  std::map< std::string, TTree* >::iterator iter;
  for ( iter = metTree_.begin(); iter != metTree_.end(); iter++ ) {
    (iter->second)->Fill();
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiMETAnalyzer::beginJob()
{
  firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ---------------------
void
HiMETAnalyzer::endJob() 
{
}

//
// constructors and destructor
//
HiMETEvent::HiMETEvent()
{
}

HiMETEvent::~HiMETEvent()
{
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Fill(const edm::Event& iEvent)
{
  this->Event_nRun   = iEvent.id().run();
  this->Event_nLumi  = iEvent.luminosityBlock();
  this->Event_nBX    = std::max(iEvent.bunchCrossing(),0);
  this->Event_nEvent = iEvent.id().event();
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Fill(const pat::MET& patMET, std::string corrLevel)
{
  std::map< std::string, pat::MET::METCorrectionLevel>::iterator iCorrection;
  std::map< std::string, pat::MET::METUncertainty>::iterator    iUncertainty;
  for ( iCorrection = METCorrectionLevelMap.begin(); iCorrection != METCorrectionLevelMap.end(); iCorrection++ ) {
    if ( (corrLevel == iCorrection->first) || (corrLevel == "All") ) {
      std::map< std::string , float > tmpSumEt;
      std::map< std::string , TVector2 > tmpP2;
      for ( iUncertainty = METUncertaintyMap.begin(); iUncertainty != METUncertaintyMap.end(); iUncertainty++ ) {
        tmpSumEt[iUncertainty->first] = patMET.shiftedSumEt(iUncertainty->second, iCorrection->second);
        reco::Candidate::LorentzVector p4 = patMET.shiftedP4(iUncertainty->second, iCorrection->second);
        tmpP2[iUncertainty->first].Set(p4.px(), p4.py());
      }
      this->shiftedSumEt[iCorrection->first] = tmpSumEt;
      this->shiftedP2[iCorrection->first] = tmpP2;
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Fill(const reco::PFMET& pfMET)
{
  this->PF_P2.Set(pfMET.px(), pfMET.py());
  this->PF_MuonEt               = pfMET.muonEt();
  this->PF_MuonEtFraction       = pfMET.muonEtFraction();
  this->PF_ChargedEMEt          = pfMET.electronEt();
  this->PF_ChargedEMEtFraction  = pfMET.electronEtFraction();
  this->PF_NeutralEMEt          = pfMET.photonEt();
  this->PF_NeutralEMEtFraction  = pfMET.photonEtFraction();
  this->PF_ChargedHadEt         = pfMET.chargedHadronEt();
  this->PF_ChargedHadEtFraction = pfMET.chargedHadronEtFraction();
  this->PF_NeutralHadEt         = pfMET.neutralHadronEt();
  this->PF_NeutralHadEtFraction = pfMET.neutralHadronEtFraction();
  this->PF_HFEMEt               = pfMET.HFEMEt();
  this->PF_HFEMEtFraction       = pfMET.HFEMEtFraction();
  this->PF_HFHadronEt           = pfMET.HFHadronEt();
  this->PF_HFHadronEtFraction   = pfMET.HFHadronEtFraction();
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Fill(const reco::CaloMET& caloMET)
{
  this->Calo_P2.Set(caloMET.px(), caloMET.py());
  this->Calo_METInmHF           = caloMET.CaloMETInmHF();
  this->Calo_METPhiInmHF        = caloMET.CaloMETPhiInmHF();
  this->Calo_SETInmHF           = caloMET.CaloSETInmHF();
  this->Calo_METInpHF           = caloMET.CaloMETInpHF();
  this->Calo_METPhiInpHF        = caloMET.CaloMETPhiInpHF();
  this->Calo_SETInpHF           = caloMET.CaloSETInpHF();
  this->Calo_MaxEtInEmTowers    = caloMET.maxEtInEmTowers();
  this->Calo_MaxEtInHadTowers   = caloMET.maxEtInHadTowers();
  this->Calo_EtFractionHadronic = caloMET.etFractionHadronic();
  this->Calo_EMEtFraction       = caloMET.emEtFraction();
  this->Calo_HadEtInHB          = caloMET.hadEtInHB();
  this->Calo_HadEtInHE          = caloMET.hadEtInHE();
  this->Calo_HadEtInHF          = caloMET.hadEtInHF();
  this->Calo_HadEtInHO          = caloMET.hadEtInHO();
  this->Calo_EMEtInEB           = caloMET.emEtInEB();
  this->Calo_EMEtInEE           = caloMET.emEtInEE();
  this->Calo_EMEtInHF           = caloMET.emEtInHF();
  this->Calo_MetSignificance    = caloMET.metSignificance();
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Fill(const reco::GenMET& genMET)
{
  this->Gen_P2.Set(genMET.px(), genMET.py());
  this->Gen_InvisibleEt          = genMET.InvisibleEt();
  this->Gen_InvisibleEtFraction  = genMET.InvisibleEtFraction();
  this->Gen_MuonEt               = genMET.MuonEt();
  this->Gen_MuonEtFraction       = genMET.MuonEtFraction();
  this->Gen_ChargedEMEt          = genMET.ChargedEMEt();
  this->Gen_ChargedEMEtFraction  = genMET.ChargedEMEtFraction();
  this->Gen_NeutralEMEt          = genMET.NeutralEMEt();
  this->Gen_NeutralEMEtFraction  = genMET.NeutralEMEtFraction();
  this->Gen_ChargedHadEt         = genMET.ChargedHadEt();
  this->Gen_ChargedHadEtFraction = genMET.ChargedHadEtFraction();
  this->Gen_NeutralHadEt         = genMET.NeutralHadEt();
  this->Gen_NeutralHadEtFraction = genMET.NeutralHadEtFraction();
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Fill(const std::map< std::string, Bool_t >& metFilter)
{
  std::map< std::string, Bool_t >::const_iterator iFilter;
  for ( iFilter = metFilter.begin(); iFilter != metFilter.end(); iFilter++ ) {
    this->Filter[iFilter->first] = iFilter->second;
  }
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::FillReco(const reco::MET& recoMET)
{
  this->Reco_P2.Set(recoMET.px(), recoMET.py());
  if (this->Reco_SigMatrix.GetNoElements()!=4) { this->Reco_SigMatrix.ResizeTo(2,2); }
  for (int i=0; i<2; i++) { for (int j=0; j<2; j++) this->Reco_SigMatrix(i,j) = recoMET.getSignificanceMatrix().At(i,j); }
  this->Reco_SumEt        = recoMET.sumEt();
  this->Reco_Significance = recoMET.significance();
  this->Reco_mEtSig       = recoMET.mEtSig();
}

//--------------------------------------------------------------------------------------------------
void 
HiMETEvent::SetBranches(const std::string name, const StringBoolMap& doMET)
{
  if ( doMET.at("Event") && (name == "Event" || name == "All") ) {
    this->tree->Branch("Event_Run",          &(this->Event_nRun),               "Event_Run/i");
    this->tree->Branch("Event_Lumi",         &(this->Event_nLumi),              "Event_Lumi/s");
    this->tree->Branch("Event_Bx",           &(this->Event_nBX),                "Event_Bx/i");
    this->tree->Branch("Event_Number",       &(this->Event_nEvent),             "Event_Number/i");
  }
  if ( doMET.at("Reco") && (name == "Reco" || name == "All") ) {
    this->tree->Branch("Reco_MET_Mom",       "TVector2",                        &(this->Reco_P2));
    this->tree->Branch("Reco_MET_SigM",      "TMatrixD",                        &(this->Reco_SigMatrix));
    this->tree->Branch("Reco_MET_Sig",       &(this->Reco_Significance),        "Reco_MET_Sig/F");
    this->tree->Branch("Reco_MET_sumEt",     &(this->Reco_SumEt),               "Reco_MET_sumEt/F");
    this->tree->Branch("Reco_MET_mEtSig",    &(this->Reco_mEtSig),              "Reco_MET_mEtSig/F");
  }
  if ( doMET.at("PF") && (name == "PF" || name == "All") ) {
    this->tree->Branch("PF_MET_Mom",             "TVector2",                         &(this->PF_P2));
    this->tree->Branch("PF_MET_Muon_Et",         &(this->PF_MuonEt),                "PF_MET_Muon_Et/F");
    this->tree->Branch("PF_MET_Muon_EtFrac",     &(this->PF_MuonEtFraction),        "PF_MET_Muon_EtFrac/F");
    this->tree->Branch("PF_MET_EM_Chg_Et",       &(this->PF_ChargedEMEt),           "PF_MET_EM_Chg_Et/F");
    this->tree->Branch("PF_MET_EM_Chg_EtFrac",   &(this->PF_ChargedEMEtFraction),   "PF_MET_EM_Chg_EtFrac/F");
    this->tree->Branch("PF_MET_EM_Neu_Et",       &(this->PF_NeutralEMEt),           "PF_MET_EM_Neu_Et/F");
    this->tree->Branch("PF_MET_EM_Neu_EtFrac",   &(this->PF_NeutralEMEtFraction),   "PF_MET_EM_Neu_EtFrac/F");
    this->tree->Branch("PF_MET_EM_HF_Et",        &(this->PF_HFEMEt),                "PF_MET_EM_HF_Et/F");
    this->tree->Branch("PF_MET_EM_HF_EtFrac",    &(this->PF_HFEMEtFraction),        "PF_MET_EM_HF_EtFrac/F");
    this->tree->Branch("PF_MET_Had_Chg_Et",      &(this->PF_ChargedHadEt),          "PF_MET_Had_Chg_Et/F");
    this->tree->Branch("PF_MET_Had_Chg_EtFrac",  &(this->PF_ChargedHadEtFraction),  "PF_MET_Had_Chg_EtFrac/F");
    this->tree->Branch("PF_MET_Had_Neu_Et",      &(this->PF_NeutralHadEt),          "PF_MET_Had_Neu_Et/F");
    this->tree->Branch("PF_MET_Had_Neu_EtFrac",  &(this->PF_NeutralHadEtFraction),  "PF_MET_Had_Neu_EtFrac/F");
    this->tree->Branch("PF_MET_Had_HF_Et",       &(this->PF_HFHadronEt),            "PF_MET_Had_HF_Et/F");
    this->tree->Branch("PF_MET_Had_HF_EtFrac",   &(this->PF_HFHadronEtFraction),    "PF_MET_Had_HF_EtFrac/F");
  }
  if ( doMET.at("Calo") && (name == "Calo" || name == "All") ) {
    this->tree->Branch("Calo_MET_Mom",           "TVector2",                        &(this->Calo_P2));
    this->tree->Branch("Calo_MET_Sig",           &(this->Calo_MetSignificance),     "Calo_MET_Sig/F");
    this->tree->Branch("Calo_MET_mHF_Et",        &(this->Calo_METInmHF),            "Calo_MET_mHF_Et/F");
    this->tree->Branch("Calo_MET_mHF_Phi",       &(this->Calo_METPhiInmHF),         "Calo_MET_mHF_Phi/F");
    this->tree->Branch("Calo_MET_mHF_sumEt",     &(this->Calo_SETInmHF),            "Calo_MET_mHF_sumEt/F");
    this->tree->Branch("Calo_MET_pHF_Et",        &(this->Calo_METInpHF),            "Calo_MET_pHF_Et/F");
    this->tree->Branch("Calo_MET_pHF_Phi",       &(this->Calo_METPhiInpHF),         "Calo_MET_pHF_Phi/F");
    this->tree->Branch("Calo_MET_pHF_sumEt",     &(this->Calo_SETInpHF),            "Calo_MET_pHF_sumEt/F");
    this->tree->Branch("Calo_MET_EM_EtFrac",     &(this->Calo_EMEtFraction),        "Calo_MET_EM_EtFrac/F");
    this->tree->Branch("Calo_MET_EM_EB_Et",      &(this->Calo_EMEtInEB),            "Calo_MET_EM_EB_Et/F");
    this->tree->Branch("Calo_MET_EM_EE_Et",      &(this->Calo_EMEtInEE),            "Calo_MET_EM_EE_Et/F");
    this->tree->Branch("Calo_MET_EM_HF_Et",      &(this->Calo_EMEtInHF),            "Calo_MET_EM_HF_Et/F");
    this->tree->Branch("Calo_MET_EM_Tow_maxEt",  &(this->Calo_MaxEtInEmTowers),     "Calo_MET_EM_Tow_maxEt/F");
    this->tree->Branch("Calo_MET_Had_EtFrac",    &(this->Calo_EtFractionHadronic),  "Calo_MET_Had_EtFrac/F");
    this->tree->Branch("Calo_MET_Had_HB_Et",     &(this->Calo_HadEtInHB),           "Calo_MET_Had_HB_Et/F");
    this->tree->Branch("Calo_MET_Had_HE_Et",     &(this->Calo_HadEtInHE),           "Calo_MET_Had_HE_Et/F");
    this->tree->Branch("Calo_MET_Had_HF_Et",     &(this->Calo_HadEtInHF),           "Calo_MET_Had_HF_Et/F");
    this->tree->Branch("Calo_MET_Had_HO_Et",     &(this->Calo_HadEtInHO),           "Calo_MET_Had_HO_Et/F");
    this->tree->Branch("Calo_MET_Had_Tow_maxEt", &(this->Calo_MaxEtInHadTowers),    "Calo_MET_Had_Tow_maxEt/F");
  }
  if ((doMET.at("Pat") && (name != "Event" && name != "Reco" && name != "Gen"))) {
    std::map< std::string, pat::MET::METCorrectionLevel >::iterator iCorrection;
    std::map< std::string, pat::MET::METUncertainty >::iterator    iUncertainty;
    for ( iCorrection = METCorrectionLevelMap.begin(); iCorrection != METCorrectionLevelMap.end(); iCorrection++ ) {
      std::string nameC = iCorrection->first;
      if ( (name == nameC) || (name == "All") ) {
        this->tree->Branch(Form("%s_MET_NoShift_Mom", nameC.c_str()),   "TVector2",                              &(this->shiftedP2[nameC]["NoShift"]));
        this->tree->Branch(Form("%s_MET_NoShift_sumEt", nameC.c_str()), &(this->shiftedSumEt[nameC]["NoShift"]), Form("%s_MET_NoShift_sumEt/F", nameC.c_str()));
        for ( iUncertainty = METUncertaintyMap.begin(); iUncertainty != METUncertaintyMap.end(); iUncertainty++ ) {
          if (iUncertainty->first == "NoShift") continue;
          this->tree->Branch(Form("%s_MET_%s_Mom", nameC.c_str(), (iUncertainty->first).c_str()),   "TVector2",                                        &(this->shiftedP2[nameC][iUncertainty->first]));
          this->tree->Branch(Form("%s_MET_%s_sumEt", nameC.c_str(), (iUncertainty->first).c_str()), &(this->shiftedSumEt[nameC][iUncertainty->first]), Form("%s_MET_%s_sumEt/F", nameC.c_str(), (iUncertainty->first).c_str()));
        }
      }
    }
  }
  if ( doMET.at("Gen") && (name == "Gen" || name == "All") ) {
    this->tree->Branch("Gen_MET_Mom",            "TVector2",                        &(this->Gen_P2));
    this->tree->Branch("Gen_MET_Inv_Et",         &(this->Gen_InvisibleEt),          "Gen_MET_Inv_Et/F");
    this->tree->Branch("Gen_MET_Inv_EtFrac",     &(this->Gen_InvisibleEtFraction),  "Gen_MET_Inv_EtFrac/F");
    this->tree->Branch("Gen_MET_Muon_Et",        &(this->Gen_MuonEt),               "Gen_MET_Muon_Et/F");
    this->tree->Branch("Gen_MET_Muon_EtFrac",    &(this->Gen_MuonEtFraction),       "Gen_MET_Muon_EtFrac/F");
    this->tree->Branch("Gen_MET_EM_Chg_Et",      &(this->Gen_ChargedEMEt),          "Gen_MET_EM_Chg_Et/F");
    this->tree->Branch("Gen_MET_EM_Chg_EtFrac",  &(this->Gen_ChargedEMEtFraction),  "Gen_MET_EM_Chg_EtFrac/F");
    this->tree->Branch("Gen_MET_EM_Neu_Et",      &(this->Gen_NeutralEMEt),          "Gen_MET_EM_Neu_Et/F");
    this->tree->Branch("Gen_MET_EM_Neu_EtFrac",  &(this->Gen_NeutralEMEtFraction),  "Gen_MET_EM_Neu_EtFrac/F");
    this->tree->Branch("Gen_MET_Had_Chg_Et",     &(this->Gen_ChargedHadEt),         "Gen_MET_Had_Chg_Et/F");
    this->tree->Branch("Gen_MET_Had_Chg_EtFrac", &(this->Gen_ChargedHadEtFraction), "Gen_MET_Had_Chg_EtFrac/F");
    this->tree->Branch("Gen_MET_Had_Neu_Et",     &(this->Gen_NeutralHadEt),         "Gen_MET_Had_Neu_Et/F");
    this->tree->Branch("Gen_MET_Had_Neu_EtFrac", &(this->Gen_NeutralHadEtFraction), "Gen_MET_Had_Neu_EtFrac/F");
  }
  if ( doMET.at("Filter") && (name == "Filter" || name == "All") ) {
    std::map< std::string, Bool_t >::iterator iFilter;
    for ( iFilter = this->Filter.begin(); iFilter != this->Filter.end(); iFilter++ ) {
      std::string nameF = iFilter->first;
      this->tree->Branch(Form("Flag_%s", nameF.c_str()),   &(this->Filter[nameF]),  Form("Flag_%s/O", nameF.c_str()));
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
HiMETEvent::Clear(void)
{
  this->Reco_SigMatrix.Clear();
  this->Reco_P2                  = TVector2();
  this->PF_P2                    = TVector2();
  this->Calo_P2                  = TVector2();
  this->Gen_P2                   = TVector2();
  this->Event_nRun               = 0;
  this->Event_nLumi              = 0;
  this->Event_nBX                = 0;
  this->Event_nEvent             = 0;
  this->Reco_SumEt               = -1.;
  this->Reco_Significance        = -99.;
  this->Reco_mEtSig              = -1.;
  this->PF_MuonEt                = -1.;
  this->PF_MuonEtFraction        = -1.;
  this->PF_ChargedEMEt           = -1.;
  this->PF_ChargedEMEtFraction   = -1.;
  this->PF_NeutralEMEt           = -1.;
  this->PF_NeutralEMEtFraction   = -1.;
  this->PF_ChargedHadEt          = -1.;
  this->PF_ChargedHadEtFraction  = -1.;
  this->PF_NeutralHadEt          = -1.;
  this->PF_NeutralHadEtFraction  = -1.;
  this->PF_HFEMEt                = -1.;
  this->PF_HFEMEtFraction        = -1.;
  this->PF_HFHadronEt            = -1.;
  this->PF_HFHadronEtFraction    = -1.;
  this->Calo_METInmHF            = -1.;
  this->Calo_METPhiInmHF         = -99.;
  this->Calo_SETInmHF            = -1.;
  this->Calo_METInpHF            = -1.;
  this->Calo_METPhiInpHF         = -99.;
  this->Calo_SETInpHF            = -1.;
  this->Calo_MaxEtInEmTowers     = -1.;
  this->Calo_MaxEtInHadTowers    = -1.;
  this->Calo_EtFractionHadronic  = -1.;
  this->Calo_EMEtFraction        = -1.;
  this->Calo_HadEtInHB           = -1.;
  this->Calo_HadEtInHE           = -1.;
  this->Calo_HadEtInHF           = -1.;
  this->Calo_HadEtInHO           = -1.;
  this->Calo_EMEtInEB            = -1.;
  this->Calo_EMEtInEE            = -1.;
  this->Calo_EMEtInHF            = -1.;
  this->Calo_MetSignificance     = -1.;
  this->Gen_InvisibleEt          = -1.;
  this->Gen_InvisibleEtFraction  = -1.;
  this->Gen_MuonEt               = -1.;
  this->Gen_MuonEtFraction       = -1.;
  this->Gen_ChargedEMEt          = -1.;
  this->Gen_ChargedEMEtFraction  = -1.;
  this->Gen_NeutralEMEt          = -1.;
  this->Gen_NeutralEMEtFraction  = -1.;
  this->Gen_ChargedHadEt         = -1.;
  this->Gen_ChargedHadEtFraction = -1.;
  this->Gen_NeutralHadEt         = -1.;
  this->Gen_NeutralHadEtFraction = -1.;
  std::map< std::string, Bool_t >::iterator iFilter;
  for ( iFilter = this->Filter.begin(); iFilter != this->Filter.end(); iFilter++ ) {
    iFilter->second = false;
  }
}


DEFINE_FWK_MODULE(HiMETAnalyzer);

