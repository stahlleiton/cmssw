// -*- C++ -*-
//
// Package:    HiMuonAnalyzer
// Class:      HiMuonAnalyzer
// 
/**\class HiMuonAnalyzer HIMuonAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Andre, Stahl 
//         Created:  Feb  08 2017
// 
//
//

#include "HeavyIonsAnalysis/ConversionAnalysis/interface/HiMuonAnalyzer.h"

//Headers for the core items
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"

//Headers for the data items
#include "DataFormats/TrackReco/interface/TrackBase.h"

//Headers for services and tools
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"

//system include files
#include <iostream>


//
// constructors and destructor
//
HiMuonAnalyzer::HiMuonAnalyzer(const edm::ParameterSet& iConfig):
  _patMuonsToken       ( consumes< edm::View< pat::Muon          > >( iConfig.getParameter< edm::InputTag >( "patMuonsTag"       ) )),
  _recoMuonsToken      ( consumes< edm::View< reco::Muon         > >( iConfig.getParameter< edm::InputTag >( "recoMuonsTag"      ) )),
  _genParticlesToken   ( consumes< reco::GenParticleCollection     >( iConfig.getParameter< edm::InputTag >( "genParticlesTag"   ) )),
  _pfCandidatesToken   ( consumes< edm::View< reco::PFCandidate  > >( iConfig.getParameter< edm::InputTag >( "pfCandidatesTag"   ) )),
  _primaryVertexToken  ( consumes< edm::View< reco::Vertex       > >( iConfig.getParameter< edm::InputTag >( "primaryVertexTag"  ) )),
  _beamSpotToken       ( consumes< reco::BeamSpot                  >( iConfig.getParameter< edm::InputTag >( "beamSpotTag"       ) )),
  _triggerResultsToken ( consumes< edm::TriggerResults             >( iConfig.getParameter< edm::InputTag >( "triggerResultsTag" ) )),
  _triggerPathNames    ( iConfig.getParameter< std::vector< std::string > >( "triggerPathNames"   ) ),
  _triggerFilterNames  ( iConfig.getParameter< std::vector< std::string > >( "triggerFilterNames" ) ),
  _hltPrescaleProvider ( iConfig, consumesCollector(), *this     ),
  _doAll               ( iConfig.getParameter< bool >( "doAll" ) )
{
}

HiMuonAnalyzer::~HiMuonAnalyzer()
{
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
HiMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  hStats_->Fill(0); // Fill "All" bin for each event
  StringBoolMap doMuon;

  // Fill Event Tree always
  doMuon["Event"] = true;

  // Extract the Trigger Results
  edm::Handle< edm::TriggerResults > triggerResultsHandle;
  getCollection(iEvent, _triggerResultsToken, triggerResultsHandle);

  // Extract the BeamSpot
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  getCollection(iEvent, _beamSpotToken, beamSpotHandle);
  reco::Vertex beamSpot = reco::Vertex();
  if (beamSpotHandle.isValid()) beamSpot = reco::Vertex( beamSpotHandle->position(),  beamSpotHandle->covariance3D() );

  // Extract the Offline Primary Vertices
  edm::Handle< edm::View< reco::Vertex > > primaryVertexHandle;
  getCollection(iEvent, _primaryVertexToken, primaryVertexHandle);
  reco::VertexCollection primaryVertexCollection;
  if (primaryVertexHandle.isValid()) {
    for (uint i = 0; i < primaryVertexHandle->size(); i++ ) {
      const reco::Vertex& primaryVertex = primaryVertexHandle->at(i);
      if ( (primaryVertex.isFake() == false) && (primaryVertex.tracksSize() >= 2) && 
           (fabs(primaryVertex.z()) <= 25) && (primaryVertex.position().rho() <= 2) ) {
        primaryVertexCollection.push_back( primaryVertex );
      }
    }
  }
  const reco::Vertex& priVtx = ( ( primaryVertexCollection.size() > 0 ) ? primaryVertexCollection.at(0) : beamSpot );

  // If Gen Particle collection exist, fill the Gen Muon tree
  edm::Handle< reco::GenParticleCollection > genParticleHandle;
  getCollection(iEvent, _genParticlesToken, genParticleHandle);
  reco::GenParticleRefVector genParticleRefCollection;
  if (genParticleHandle.isValid()) {
    for (auto it = genParticleHandle->begin(); it != genParticleHandle->end(); ++it) {
      const reco::GenParticle& genParticle = *it;
      if ( ( (abs(genParticle.pdgId()) >  10) && (abs(genParticle.pdgId()) < 19) ) || // Keep all leptons
           ( (abs(genParticle.pdgId()) >  21) && (abs(genParticle.pdgId()) < 25) )    // Keep all EWK bosons
           ){
        genParticleRefCollection.push_back( reco::GenParticleRef(genParticleHandle, it - genParticleHandle->begin()) );
      }
    }
    doMuon["Gen"] = true;
  }
  else { doMuon["Gen"] = false; }
  reco::GenParticleCollection genMuonCollection;
  for (uint i = 0; i < genParticleRefCollection.size(); i++ ) {
    if ( (abs(genParticleRefCollection.at(i)->pdgId()) == 13) && (genParticleRefCollection.at(i)->status()==1) ) genMuonCollection.push_back( *(genParticleRefCollection.at(i)) );
  }

  // If Reco Muon collection exist, fill the Reco Muon tree
  edm::Handle< edm::View< reco::Muon > > recoMuonHandle;
  getCollection(iEvent, _recoMuonsToken, recoMuonHandle);
  reco::MuonCollection recoMuonCollection;
  if (recoMuonHandle.isValid()) {
    for (uint i = 0; i < recoMuonHandle->size(); i++ ) {
      const reco::Muon& recoMuon = recoMuonHandle->at(i);
      if ( recoMuon.isGlobalMuon() || recoMuon.isTrackerMuon() || recoMuon.isPFMuon() ) {
        recoMuonCollection.push_back( recoMuon );
      }
    }
    doMuon["Reco"] = true;
  }
  else { doMuon["Reco"] = false; }
  std::sort(recoMuonCollection.begin(), recoMuonCollection.end(), recoMuonPTComparator_);

  // If PF Particle collection exist, fill the PF Muon tree
  edm::Handle<  edm::View<reco::PFCandidate> > pfCandidateHandle;
  getCollection(iEvent, _pfCandidatesToken, pfCandidateHandle);
  reco::PFCandidateCollection pfCandidateCollection;
  if (pfCandidateHandle.isValid()) {
    for (uint i = 0; i < pfCandidateHandle->size(); i++ ) {
      const reco::PFCandidate& pfCandidate = pfCandidateHandle->at(i);
      pfCandidateCollection.push_back( pfCandidate );
    }
    doMuon["PF"] = true; 
  }
  else  { doMuon["PF"] = false; }
  std::sort(pfCandidateCollection.begin(), pfCandidateCollection.end(), pfCandidatePTComparator_);
  reco::PFCandidateCollection pfMuonCollection;
  for (uint i = 0; i < pfCandidateCollection.size(); i++ ) {
    if (pfCandidateCollection.at(i).particleId() == reco::PFCandidate::mu) pfMuonCollection.push_back( pfCandidateCollection.at(i) );
  }

  // If PAT Muon collection exist, fill the PAT Muon tree
  edm::Handle< edm::View< pat::Muon > > patMuonHandle;
  getCollection(iEvent, _patMuonsToken, patMuonHandle);
  pat::MuonCollection patMuonCollection;
  if (patMuonHandle.isValid() ) {
    for (uint i = 0; i < patMuonHandle->size(); i++ ) {
      const pat::Muon& patMuon = patMuonHandle->at(i);
      if ( patMuon.isGlobalMuon() || patMuon.isTrackerMuon() || patMuon.isPFMuon() ) {
        patMuonCollection.push_back( patMuon );
      }
    }
    doMuon["Pat"]  = true;
    doMuon["Reco"] = true;
  }
  else  { doMuon["Pat"] = false; }
  std::sort(patMuonCollection.begin(), patMuonCollection.end(), patMuonPTComparator_);
  reco::GenParticleCollection patGenMuonCollection;
  for (uint i = 0; i < patMuonCollection.size(); i++ ) {
    reco::GenParticleRef genRef = patMuonCollection.at(i).genParticleById(13, 1, true);
    if (genRef.isNonnull() && genRef.isAvailable()) patGenMuonCollection.push_back( *genRef );
    else patGenMuonCollection.push_back( reco::GenParticle() );
  }
  
  // Check if at least one of the collections was found
  if (!doMuon["Pat"] && !doMuon["Reco"] && !doMuon["PF"] && !doMuon["Gen"]) { 
    edm::LogError("HiMuonAnalyzer_ProductNotFound") << "No collections were found!" << std::endl;
    return;
  }

  // Do matching between muon collections
  IndexMap indexMap_GEN, indexMap_PF, indexMap_RECO;
  std::vector< std::vector< Char_t > > matchRECOGEN , matchRECOPF , matchPFGEN;
  if (doMuon["Pat"]) matchRECOGEN = doMatching(patGenMuonCollection, genMuonCollection, 0.00001, 0.00001);
  else               matchRECOGEN = doMatching(recoMuonCollection,   genMuonCollection, 0.5,     0.5    );
  if (doMuon["Pat"]) matchRECOPF  = doMatching(patMuonCollection,    pfMuonCollection,  0.01,    0.01   );
  else               matchRECOPF  = doMatching(recoMuonCollection,   pfMuonCollection,  0.01,    0.01   );
  matchPFGEN   = doMatching(pfMuonCollection,     genMuonCollection, 0.5,     0.5    );
  indexMap_RECO["GEN"] = matchRECOGEN[0];
  indexMap_GEN["RECO"] = matchRECOGEN[1];
  indexMap_RECO["PF"]  = matchRECOPF[0];
  indexMap_PF["RECO"]  = matchRECOPF[1];
  indexMap_PF["GEN"]   = matchPFGEN[0];
  indexMap_GEN["PF"]   = matchPFGEN[1];

  std::vector< std::string > treeNames;
  if (_doAll) { treeNames.push_back("All"); }
  else {
    if (doMuon["Event"]) treeNames.push_back("Event");
    if (doMuon["Reco"])  treeNames.push_back("Reco");
    if (doMuon["PF"])    treeNames.push_back("PF");
    if (doMuon["Gen"])   treeNames.push_back("Gen");
  }

  for (uint i = 0; i < treeNames.size(); i++) {
    std::string name = treeNames[i];
    if (firstEvent_) muonEvt_[name].IniArrays();
    muonEvt_[name].Clear();
    if (doMuon["Event"] && (name=="Event" || name=="All")) { 
      muonEvt_[name].SetPrimaryVertex       ( priVtx );
      muonEvt_[name].SetHLTTriggerNames     ( _triggerPathNames     );
      muonEvt_[name].SetHLTPrescaleProvider ( &_hltPrescaleProvider );
      muonEvt_[name].SetHLTTriggerResults   ( triggerResultsHandle  );
      muonEvt_[name].Fill( iEvent, iSetup, primaryVertexCollection  ); 
    }
    if (doMuon["Reco"]  && (name=="Reco" || name=="All")) { 
      muonEvt_[name].IniTrackBuilder   ( iSetup );
      muonEvt_[name].SetHLTFilterNames ( _triggerFilterNames );
      muonEvt_[name].SetPrimaryVertex  ( priVtx );
      if (doMuon["Pat"]) { muonEvt_[name].Fill ( patMuonCollection,  indexMap_RECO ); }
      else               { muonEvt_[name].Fill ( recoMuonCollection, indexMap_RECO ); }
    }
    if (doMuon["PF"]   && (name=="PF" || name=="All")) { 
      muonEvt_[name].IniTrackBuilder  ( iSetup );
      muonEvt_[name].SetPrimaryVertex ( priVtx );
      muonEvt_[name].Fill( pfCandidateCollection, indexMap_PF, primaryVertexCollection ); 
    }
    if (doMuon["Gen"]  && (name=="Gen" || name=="All")) { 
      muonEvt_[name].Fill( genParticleRefCollection, indexMap_GEN ); 
    }
    if (firstEvent_) {
      muonTree_[name] = fs_->make<TTree>(Form("Muon_%s", name.c_str()), "");
      muonEvt_[name].SetTree( muonTree_[name] );
      muonEvt_[name].SetBranches( name , doMuon);
    }
  }
  if (firstEvent_) { firstEvent_ = false; }

  std::map< std::string, TTree* >::iterator iter;
  for ( iter = muonTree_.begin(); iter != muonTree_.end(); iter++ ) {
    (iter->second)->Fill();
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiMuonAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(_triggerResultsToken, triggerResultsLabel);
  _hltPrescaleProvider.init(iRun, iSetup, triggerResultsLabel.process, changed);
}

//--------------------------------------------------------------------------------------------------
void 
HiMuonAnalyzer::beginJob()
{
  firstEvent_ = true;
  uint NTRIGGERS = _triggerPathNames.size();
  hStats_ = fs_->make<TH1F>("hStats","hStats;;Number of Events",2*NTRIGGERS+2,0,2*NTRIGGERS+2);
  hStats_->GetXaxis()->SetBinLabel( 1 , "All_Events" );
  hStats_->GetXaxis()->SetBinLabel( 2+NTRIGGERS, "All_Muons" );
  for (uint i = 0; i < NTRIGGERS; ++i ) {
    hStats_->GetXaxis()->SetBinLabel( (i+2) , _triggerPathNames.at(i).c_str() ); // event trigger firing
    hStats_->GetXaxis()->SetBinLabel( (i+3+NTRIGGERS) , _triggerPathNames.at(i).c_str() ); // muon trigger matching
  }
  hStats_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ---------------------
void
HiMuonAnalyzer::endJob() 
{
}

//
// constructors and destructor
//
HiMuonEvent::HiMuonEvent()
{
}

HiMuonEvent::~HiMuonEvent()
{
}

//--------------------------------------------------------------------------------------------------
void
HiMuonEvent::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::VertexCollection& vtxCol)
{
  this->Event_nRun   = iEvent.id().run();
  this->Event_nLumi  = iEvent.luminosityBlock();
  this->Event_nBX    = std::max(iEvent.bunchCrossing(),0);
  this->Event_nOrbit = std::max(iEvent.orbitNumber(),0);
  this->Event_nEvent = iEvent.id().event();
  // Primary Vertex Information
  this->Event_PriVtx_Position = TVector3( privtx_.x()      , privtx_.y()      , privtx_.z()      );
  this->Event_PriVtx_Error    = TVector3( privtx_.xError() , privtx_.yError() , privtx_.zError() );
  this->Event_nPV = vtxCol.size();
  // Trigger Information
  if (triggerResultsHandle_.isValid()) {
    edm::TriggerNames const& triggerNames = iEvent.triggerNames(*triggerResultsHandle_);
    for (ushort iTr = 0; iTr < triggerNames_.size(); iTr++) {
      bool isTriggerFired = false;
      uint triggerIndex = triggerNames.triggerIndex( triggerNames_.at(iTr) );
      if ( triggerIndex >= triggerNames.size() ) continue;
      if ( triggerResultsHandle_->accept(triggerIndex) ) isTriggerFired = true;
      int prescaleValue = -1;
      if ( hltPrescaleProvider_->hltConfigProvider().inited() && hltPrescaleProvider_->prescaleSet(iEvent,iSetup)>=0 ) {
        std::pair<std::vector<std::pair<std::string,int> >,int> presInfo = hltPrescaleProvider_->prescaleValuesInDetail(iEvent, iSetup, triggerNames_.at(iTr));
        const int hltPres = presInfo.second;
        int l1Pres = ( (presInfo.first.size()<1) ? -1 : 1 );
        if (presInfo.first.size()==1) l1Pres = presInfo.first.at(0).second;
        prescaleValue = hltPres*l1Pres;
      }
      this->Event_Trigger_isFired.push_back  ( isTriggerFired       );
      this->Event_Trigger_Prescale.push_back ( prescaleValue        );
      if (isTriggerFired) hStats_->Fill(iTr+1);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
HiMuonEvent::Fill(const pat::MuonCollection& patMuonCollection, const IndexMap& indexMap)
{
  reco::MuonCollection recoMuonCollection;
  for (uint imuon = 0; imuon < patMuonCollection.size(); imuon++) {
    recoMuonCollection.push_back( patMuonCollection.at(imuon) );
  }
  for (ushort imuon = 0; imuon < patMuonCollection.size(); imuon++) { 
    const pat::Muon& patMuon = patMuonCollection.at(imuon);
    // Pat Muon Information
    std::vector<unsigned char> trigBits;
    for (ushort iTr = 0; iTr < filterNames_.size(); iTr++) {
      const pat::TriggerObjectStandAloneCollection muHLTMatchesFilter = patMuon.triggerObjectMatchesByFilter( filterNames_.at(iTr) );
      trigBits.push_back( muHLTMatchesFilter.size() > 0 );
      if (muHLTMatchesFilter.size() > 0) { hStats_->Fill(2+filterNames_.size()+iTr); }
    }
    hStats_->Fill(1+filterNames_.size());
    this->Pat_Muon_dB.push_back     ( patMuon.dB  ( pat::Muon::PV2D ) );
    this->Pat_Muon_dBErr.push_back  ( patMuon.edB ( pat::Muon::PV2D ) );
    this->Pat_Muon_TriggerMatched.push_back( trigBits );
    // Fill Reco Muon
    this->Fill(recoMuonCollection, indexMap, imuon);
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiMuonEvent::IniTrackBuilder(const edm::EventSetup& iSetup)
{
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", _theTTBuilder);
}

//--------------------------------------------------------------------------------------------------
void
HiMuonEvent::Fill(const reco::MuonCollection& recoMuonCollection, const IndexMap& indexMap, const short muonIndex)
{
  this->Reco_Muon_N = recoMuonCollection.size();
  for (ushort imuon = 0; imuon < recoMuonCollection.size(); imuon++) {
    if ( (muonIndex >= 0) && (imuon != muonIndex) ) continue;
    const reco::Muon& recoMuon = recoMuonCollection.at(imuon);
    // Reco Muon Information
    TLorentzVector muP4 = TLorentzVector();
    muP4.SetPtEtaPhiM( recoMuon.pt() , recoMuon.eta() , recoMuon.phi() , recoMuon.mass() );
    new ((*this->Reco_Muon_P4)[imuon]) TLorentzVector( muP4 );
    this->Reco_Muon_Charge.push_back    ( recoMuon.charge()         );
    this->Reco_Muon_Gen_Index.push_back ( indexMap.at("GEN")[imuon] );
    this->Reco_Muon_PF_Index.push_back  ( indexMap.at("PF")[imuon]  );
    // Reco Muon Flags
    this->Reco_Muon_isPFMuon.push_back                 ( recoMuon.isPFMuon()                       );
    this->Reco_Muon_isGlobalMuon.push_back             ( recoMuon.isGlobalMuon()                   );
    this->Reco_Muon_isTrackerMuon.push_back            ( recoMuon.isTrackerMuon()                  );
    this->Reco_Muon_isStandAloneMuon.push_back         ( recoMuon.isStandAloneMuon()               );
    this->Reco_Muon_isLooseMuon.push_back              ( muon::isLooseMuon  ( recoMuon )           );
    this->Reco_Muon_isMediumMuon.push_back             ( muon::isMediumMuon ( recoMuon )           );
    this->Reco_Muon_isHighPtMuon.push_back             ( muon::isHighPtMuon ( recoMuon , privtx_ ) );
    this->Reco_Muon_isSoftMuon.push_back               ( muon::isSoftMuon   ( recoMuon , privtx_ ) );
    this->Reco_Muon_isTightMuon.push_back              ( muon::isTightMuon  ( recoMuon , privtx_ ) );
    this->Reco_Muon_isArbitrated.push_back             ( muon::isGoodMuon( recoMuon , muon::AllArbitrated          ) );
    this->Reco_Muon_TrackerMuonArbitrated.push_back    ( muon::isGoodMuon( recoMuon , muon::TrackerMuonArbitrated  ) );
    this->Reco_Muon_TMLastStationLoose.push_back       ( muon::isGoodMuon( recoMuon , muon::TMLastStationLoose     ) );
    this->Reco_Muon_TMLastStationTight.push_back       ( muon::isGoodMuon( recoMuon , muon::TMLastStationTight     ) );
    this->Reco_Muon_TM2DCompatibilityLoose.push_back   ( muon::isGoodMuon( recoMuon , muon::TM2DCompatibilityLoose ) );
    this->Reco_Muon_TM2DCompatibilityTight.push_back   ( muon::isGoodMuon( recoMuon , muon::TM2DCompatibilityTight ) );
    this->Reco_Muon_TMOneStationLoose.push_back        ( muon::isGoodMuon( recoMuon , muon::TMOneStationLoose      ) );
    this->Reco_Muon_TMOneStationTight.push_back        ( muon::isGoodMuon( recoMuon , muon::TMOneStationTight      ) );
    this->Reco_Muon_GlobalMuonPromptTight.push_back    ( muon::isGoodMuon( recoMuon , muon::GlobalMuonPromptTight  ) );
    this->Reco_Muon_GMTkChiCompatibility.push_back     ( muon::isGoodMuon( recoMuon , muon::GMTkChiCompatibility   ) );
    this->Reco_Muon_GMStaChiCompatibility.push_back    ( muon::isGoodMuon( recoMuon , muon::GMStaChiCompatibility  ) );
    this->Reco_Muon_GMTkKinkTight.push_back            ( muon::isGoodMuon( recoMuon , muon::GMTkKinkTight          ) );
    this->Reco_Muon_TMLastStationAngLoose.push_back    ( muon::isGoodMuon( recoMuon , muon::TMLastStationAngLoose  ) );
    this->Reco_Muon_TMLastStationAngTight.push_back    ( muon::isGoodMuon( recoMuon , muon::TMLastStationAngTight  ) );
    this->Reco_Muon_TMOneStationAngLoose.push_back     ( muon::isGoodMuon( recoMuon , muon::TMOneStationAngLoose   ) );
    this->Reco_Muon_TMOneStationAngTight.push_back     ( muon::isGoodMuon( recoMuon , muon::TMOneStationAngTight   ) );
    // Reco Muon Quality Variables
    this->Reco_Muon_NumberOfMatchedStations.push_back  ( recoMuon.numberOfMatchedStations()   );
    this->Reco_Muon_NumberOfMatches.push_back          ( recoMuon.numberOfMatches()           );
    this->Reco_Muon_SegmentCompatibility.push_back     ( muon::segmentCompatibility(recoMuon) );
    if (recoMuon.isQualityValid()) {
      this->Reco_Muon_Chi2LocalPosition.push_back ( recoMuon.combinedQuality().chi2LocalPosition );
      this->Reco_Muon_TrkKink.push_back           ( recoMuon.combinedQuality().trkKink           );
    }
    else {
      this->Reco_Muon_Chi2LocalPosition.push_back ( -1. );
      this->Reco_Muon_TrkKink.push_back           ( -1. );
    }
    // Reco Muon Track Information
    if (recoMuon.innerTrack().isNonnull() && recoMuon.innerTrack().isAvailable()) {
      TVector3 p3 = TVector3();
      p3.SetPtEtaPhi( recoMuon.innerTrack()->pt() , recoMuon.innerTrack()->eta() , recoMuon.innerTrack()->phi() );
      new ((*this->Reco_Muon_InnerTrack_P3)[imuon]) TVector3( p3 );
      this->Reco_Muon_InnerTrack_PtError.push_back          ( recoMuon.innerTrack()->ptError()                                   );
      this->Reco_Muon_InnerTrack_isHighPurity.push_back     ( recoMuon.innerTrack()->quality(reco::TrackBase::highPurity)        );
      this->Reco_Muon_InnerTrack_NumOfValHits.push_back     ( recoMuon.innerTrack()->numberOfValidHits()                         );
      this->Reco_Muon_InnerTrack_NumOfLostHits.push_back    ( recoMuon.innerTrack()->numberOfLostHits()                          );
      this->Reco_Muon_InnerTrack_NumOfValPixHits.push_back  ( recoMuon.innerTrack()->hitPattern().numberOfValidPixelHits()       );
      this->Reco_Muon_InnerTrack_TrkLayersWithMea.push_back ( recoMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement() );
      this->Reco_Muon_InnerTrack_PixLayersWithMea.push_back ( recoMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement()   );
      this->Reco_Muon_InnerTrack_dXY.push_back              ( recoMuon.innerTrack()->dxy(privtx_.position())                     );
      this->Reco_Muon_InnerTrack_dXYErr.push_back           ( recoMuon.innerTrack()->dxyError()                                  );
      this->Reco_Muon_InnerTrack_dZ.push_back               ( recoMuon.innerTrack()->dz(privtx_.position())                      );
      this->Reco_Muon_InnerTrack_dZErr.push_back            ( recoMuon.innerTrack()->dzError()                                   );
      this->Reco_Muon_InnerTrack_ValFrac.push_back          ( recoMuon.innerTrack()->validFraction()                             );
      this->Reco_Muon_InnerTrack_NormChi2.push_back         ( recoMuon.innerTrack()->normalizedChi2()                            );
    }
    else {
      new ((*this->Reco_Muon_InnerTrack_P3)[imuon]) TVector3();
      this->Reco_Muon_InnerTrack_PtError.push_back          ( -1. );
      this->Reco_Muon_InnerTrack_isHighPurity.push_back     (  0  );
      this->Reco_Muon_InnerTrack_NumOfValHits.push_back     ( -1  );
      this->Reco_Muon_InnerTrack_NumOfLostHits.push_back    ( -1  );
      this->Reco_Muon_InnerTrack_NumOfValPixHits.push_back  ( -1  );
      this->Reco_Muon_InnerTrack_TrkLayersWithMea.push_back ( -1  );
      this->Reco_Muon_InnerTrack_PixLayersWithMea.push_back ( -1  );
      this->Reco_Muon_InnerTrack_dXY.push_back              ( -99.);
      this->Reco_Muon_InnerTrack_dXYErr.push_back           ( -1. );
      this->Reco_Muon_InnerTrack_dZ.push_back               ( -99.);
      this->Reco_Muon_InnerTrack_dZErr.push_back            ( -1. );
      this->Reco_Muon_InnerTrack_ValFrac.push_back          ( -1. );
      this->Reco_Muon_InnerTrack_NormChi2.push_back         ( -1. );
    }
    if (recoMuon.globalTrack().isNonnull() && recoMuon.globalTrack().isAvailable()) {
      TVector3 p3 = TVector3();
      p3.SetPtEtaPhi( recoMuon.globalTrack()->pt() , recoMuon.globalTrack()->eta() , recoMuon.globalTrack()->phi() );
      new ((*this->Reco_Muon_GlobalTrack_P3)[imuon]) TVector3( p3 );
      this->Reco_Muon_GlobalTrack_PtError.push_back          ( recoMuon.globalTrack()->ptError()                            );
      this->Reco_Muon_GlobalTrack_NumOfValMuonHits.push_back ( recoMuon.globalTrack()->hitPattern().numberOfValidMuonHits() );
      this->Reco_Muon_GlobalTrack_NormChi2.push_back         ( recoMuon.globalTrack()->normalizedChi2()                     );
    }
    else {
      new ((*this->Reco_Muon_GlobalTrack_P3)[imuon]) TVector3();
      this->Reco_Muon_GlobalTrack_PtError.push_back          ( -1. );
      this->Reco_Muon_GlobalTrack_NumOfValMuonHits.push_back ( -1  );
      this->Reco_Muon_GlobalTrack_NormChi2.push_back         ( -1. );
    }
    if (recoMuon.muonBestTrack().isNonnull() && recoMuon.muonBestTrack().isAvailable()) {
      TVector3 p3 = TVector3();
      p3.SetPtEtaPhi( recoMuon.muonBestTrack()->pt() , recoMuon.muonBestTrack()->eta() , recoMuon.muonBestTrack()->phi() );
      TVector3 vertex = TVector3(recoMuon.muonBestTrack()->vertex().X(), recoMuon.muonBestTrack()->vertex().Y(), recoMuon.muonBestTrack()->vertex().Z());
      new ((*this->Reco_Muon_BestTrack_P3)[imuon]    ) TVector3( p3     );
      new ((*this->Reco_Muon_BestTrack_Vertex)[imuon]) TVector3( vertex );
      this->Reco_Muon_BestTrack_Type.push_back    ( recoMuon.muonBestTrackType()                            );
      this->Reco_Muon_BestTrack_PtError.push_back ( recoMuon.muonBestTrack()->ptError()                     );
      this->Reco_Muon_BestTrack_dXY.push_back     ( recoMuon.muonBestTrack()->dxy(privtx_.position())       );
      this->Reco_Muon_BestTrack_dXYErr.push_back  ( recoMuon.muonBestTrack()->dxyError()                    );
      this->Reco_Muon_BestTrack_dZ.push_back      ( recoMuon.muonBestTrack()->dz(privtx_.position())        );
      this->Reco_Muon_BestTrack_dZErr.push_back   ( recoMuon.muonBestTrack()->dzError()                     );
    }
    else {
      new ((*this->Reco_Muon_BestTrack_P3)[imuon]    ) TVector3();
      new ((*this->Reco_Muon_BestTrack_Vertex)[imuon]) TVector3();
      this->Reco_Muon_BestTrack_Type.push_back    ( -1  );
      this->Reco_Muon_BestTrack_PtError.push_back ( -1. );
      this->Reco_Muon_BestTrack_dXY.push_back     ( -99.);
      this->Reco_Muon_BestTrack_dXYErr.push_back  ( -1. );
      this->Reco_Muon_BestTrack_dZ.push_back      ( -99.);
      this->Reco_Muon_BestTrack_dZErr.push_back   ( -1. );
    }
    // Reco Muon Isolation
    if (recoMuon.isPFIsolationValid()) {
      Float_t isoR03PUCorr   = recoMuon.pfIsolationR03().sumChargedHadronPt + std::max(0., recoMuon.pfIsolationR03().sumNeutralHadronEt + recoMuon.pfIsolationR03().sumPhotonEt - 0.5*recoMuon.pfIsolationR03().sumPUPt );
      Float_t isoR03NoPUCorr = recoMuon.pfIsolationR03().sumChargedHadronPt + std::max(0.0f, recoMuon.pfIsolationR03().sumNeutralHadronEt + recoMuon.pfIsolationR03().sumPhotonEt );
      this->Reco_Muon_PFIsoR03_ChargedEM_SumPt.push_back    ( recoMuon.pfIsolationR03().sumChargedParticlePt - recoMuon.pfIsolationR03().sumChargedHadronPt );
      this->Reco_Muon_PFIsoR03_NeutralEM_SumEt.push_back    ( recoMuon.pfIsolationR03().sumPhotonEt          );
      this->Reco_Muon_PFIsoR03_ChargedHad_SumPt.push_back   ( recoMuon.pfIsolationR03().sumChargedHadronPt   );
      this->Reco_Muon_PFIsoR03_NeutralHad_SumEt.push_back   ( recoMuon.pfIsolationR03().sumNeutralHadronEt   );
      this->Reco_Muon_PFIsoR03_ChargedHadPU_SumPt.push_back ( recoMuon.pfIsolationR03().sumPUPt              );
      this->Reco_Muon_PFIsoR03_IsoPUCorr.push_back          ( isoR03PUCorr   / (recoMuon.pt()+1E-12)         );
      this->Reco_Muon_PFIsoR03_IsoNoPUCorr.push_back        ( isoR03NoPUCorr / (recoMuon.pt()+1E-12)         );
      Float_t isoR04PUCorr   = recoMuon.pfIsolationR04().sumChargedHadronPt + std::max(0., recoMuon.pfIsolationR04().sumNeutralHadronEt + recoMuon.pfIsolationR04().sumPhotonEt - 0.5*recoMuon.pfIsolationR04().sumPUPt );
      Float_t isoR04NoPUCorr = recoMuon.pfIsolationR04().sumChargedHadronPt + std::max(0.0f, recoMuon.pfIsolationR04().sumNeutralHadronEt + recoMuon.pfIsolationR04().sumPhotonEt );
      this->Reco_Muon_PFIsoR04_ChargedEM_SumPt.push_back    ( recoMuon.pfIsolationR04().sumChargedParticlePt - recoMuon.pfIsolationR04().sumChargedHadronPt );
      this->Reco_Muon_PFIsoR04_ChargedHad_SumPt.push_back   ( recoMuon.pfIsolationR04().sumChargedHadronPt   );
      this->Reco_Muon_PFIsoR04_NeutralHad_SumEt.push_back   ( recoMuon.pfIsolationR04().sumNeutralHadronEt   );
      this->Reco_Muon_PFIsoR04_NeutralEM_SumEt.push_back    ( recoMuon.pfIsolationR04().sumPhotonEt          );
      this->Reco_Muon_PFIsoR04_ChargedHadPU_SumPt.push_back ( recoMuon.pfIsolationR04().sumPUPt              );
      this->Reco_Muon_PFIsoR04_IsoPUCorr.push_back          ( isoR04PUCorr   / (recoMuon.pt()+1E-12)         );
      this->Reco_Muon_PFIsoR04_IsoNoPUCorr.push_back        ( isoR04NoPUCorr / (recoMuon.pt()+1E-12)         );
    }
    else {
      this->Reco_Muon_PFIsoR03_ChargedEM_SumPt.push_back    ( -1. );
      this->Reco_Muon_PFIsoR03_NeutralEM_SumEt.push_back    ( -1. );
      this->Reco_Muon_PFIsoR03_ChargedHad_SumPt.push_back   ( -1. );
      this->Reco_Muon_PFIsoR03_NeutralHad_SumEt.push_back   ( -1. );
      this->Reco_Muon_PFIsoR03_ChargedHadPU_SumPt.push_back ( -1. );
      this->Reco_Muon_PFIsoR03_IsoPUCorr.push_back          ( -1. );
      this->Reco_Muon_PFIsoR03_IsoNoPUCorr.push_back        ( -1. );
      this->Reco_Muon_PFIsoR04_ChargedEM_SumPt.push_back    ( -1. );
      this->Reco_Muon_PFIsoR04_NeutralEM_SumEt.push_back    ( -1. );
      this->Reco_Muon_PFIsoR04_ChargedHad_SumPt.push_back   ( -1. );
      this->Reco_Muon_PFIsoR04_NeutralHad_SumEt.push_back   ( -1. );
      this->Reco_Muon_PFIsoR04_ChargedHadPU_SumPt.push_back ( -1. );
      this->Reco_Muon_PFIsoR04_IsoPUCorr.push_back          ( -1. );
      this->Reco_Muon_PFIsoR04_IsoNoPUCorr.push_back        ( -1. );
    }
    if (recoMuon.isIsolationValid()) {
      this->Reco_Muon_IsoR03_SumPt.push_back ( recoMuon.isolationR03().sumPt );
      this->Reco_Muon_IsoR05_SumPt.push_back ( recoMuon.isolationR05().sumPt );
      this->Reco_Muon_IsoR03_Iso.push_back   ( recoMuon.isolationR03().sumPt / (recoMuon.pt()+1E-12) );
      this->Reco_Muon_IsoR05_Iso.push_back   ( recoMuon.isolationR05().sumPt / (recoMuon.pt()+1E-12) );
    }
    else {
      this->Reco_Muon_IsoR03_SumPt.push_back ( -1. );
      this->Reco_Muon_IsoR05_SumPt.push_back ( -1. );
      this->Reco_Muon_IsoR03_Iso.push_back   ( -1. );
      this->Reco_Muon_IsoR05_Iso.push_back   ( -1. );
    }
    // Reco Dimuon
    for (ushort imuon2 = imuon+1; imuon2 < recoMuonCollection.size(); imuon2++) {
      const reco::Muon& recoMuon2 = recoMuonCollection.at(imuon2);
      // Reco Dimuon Information
      TLorentzVector mu2P4 = TLorentzVector();
      mu2P4.SetPtEtaPhiM( recoMuon2.pt() , recoMuon2.eta() , recoMuon2.phi() , recoMuon2.mass() );
      TLorentzVector diMuonP4 = TLorentzVector();
      diMuonP4 = muP4 + mu2P4;
      Char_t charge = recoMuon.charge() + recoMuon2.charge();
      Bool_t isCowBoy = ( ( recoMuon.charge() * TVector2::Phi_mpi_pi( recoMuon.phi() - recoMuon2.phi() ) ) > 0. );
      TVector3 diMuonVtx = TVector3();
      Float_t vProb    = -1. ;
      Float_t dca      = -1.;
      Float_t massWErr = -1.;
      if ( ( recoMuon.muonBestTrack().isNonnull()  && recoMuon.muonBestTrack().isAvailable()  ) && 
           ( recoMuon2.muonBestTrack().isNonnull() && recoMuon2.muonBestTrack().isAvailable() ) ) {
        std::vector<reco::TransientTrack> t_tks;
        t_tks.push_back( _theTTBuilder->build( *recoMuon.muonBestTrack()  ) );
        t_tks.push_back( _theTTBuilder->build( *recoMuon2.muonBestTrack() ) );
        KalmanVertexFitter  vtxFitter = KalmanVertexFitter(true);
        TransientVertex diMuonVertex = vtxFitter.vertex( t_tks );
        if (diMuonVertex.isValid()) {
          Float_t vChi2 = diMuonVertex.totalChiSquared();
          Float_t vNDF  = diMuonVertex.degreesOfFreedom();
          vProb = TMath::Prob(vChi2,(int)vNDF);
          diMuonVtx.SetXYZ( diMuonVertex.position().x() , diMuonVertex.position().y() , diMuonVertex.position().z() );
        }
        TrajectoryStateClosestToPoint mu1TS = t_tks[0].impactPointTSCP();
        TrajectoryStateClosestToPoint mu2TS = t_tks[1].impactPointTSCP();
        if (mu1TS.isValid() && mu2TS.isValid()) {
          ClosestApproachInRPhi cApp;
          cApp.calculate(mu1TS.theState(), mu2TS.theState());
          if ( cApp.status() ) dca = cApp.distance();
        }
        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );
        InvariantMassFromVertex massCalculator;
        massWErr = massCalculator.invariantMass( VtxForInvMass, _muMasses ).error();
      }
      ushort idimuon = this->Reco_DiMuon_Charge.size();
      this->Reco_DiMuon_N = idimuon + 1;
      new ((*this->Reco_DiMuon_P4)[idimuon]) TLorentzVector( diMuonP4 );
      this->Reco_DiMuon_Charge.push_back      ( charge    );
      this->Reco_DiMuon_isCowBoy.push_back    ( isCowBoy  );
      new ((*this->Reco_DiMuon_Vertex)[idimuon]) TVector3( diMuonVtx );
      this->Reco_DiMuon_VertexProb.push_back  ( vProb     );
      this->Reco_DiMuon_DCA.push_back         ( dca       );
      this->Reco_DiMuon_MassError.push_back   ( massWErr  );
      this->Reco_DiMuon_Muon1_Index.push_back ( imuon     );
      this->Reco_DiMuon_Muon2_Index.push_back ( imuon2    );
    }
  }
}

//--------------------------------------------------------------------------------------------------
bool 
HiMuonEvent::isMatched(const reco::Candidate& cand1, const reco::Candidate& cand2, double maxDeltaR, double maxDPtRel)
{
  double deltaR = reco::deltaR( cand1.eta(), cand2.phi(), cand1.eta(), cand2.phi());
  double dPtRel = abs( cand1.pt() - cand2.pt() ) / ( cand2.pt() + 1E-12 );
  return ( (deltaR < maxDeltaR) && (dPtRel < maxDPtRel) && (cand1.charge() == cand2.charge()) );
}

//--------------------------------------------------------------------------------------------------
double 
HiMuonEvent::pfIsolation(const reco::Candidate& muon, const reco::PFCandidateCollection& pfCandColl, double coneSize, double vetoArea)
{
  double sumPt = 0.0;
  for (uint icand = 0; icand < pfCandColl.size(); icand++) {
    const reco::PFCandidate& pfCand = pfCandColl.at(icand);
    if ( (pfCand.particleId() == reco::PFCandidate::mu) && isMatched(muon, pfCand, 0.0000001, 0.0000001) ) continue;
    double deltaR = reco::deltaR(muon.eta(), muon.phi(), pfCand.eta(), pfCand.phi());
    if ( (deltaR < coneSize) && (deltaR >= vetoArea) ) ( (pfCand.charge()!=0) ? (sumPt += pfCand.pt()) : (sumPt += pfCand.et()) );
  }
  return sumPt;
}

//--------------------------------------------------------------------------------------------------
void
HiMuonEvent::Fill(const reco::PFCandidateCollection& pfCandidateCollection, const IndexMap& indexMap, const reco::VertexCollection& priVtxCollection)
{
  // Fill PF Muon collection
  reco::PFCandidateCollection pfMuonCollection;
  for (uint i = 0; i < pfCandidateCollection.size(); i++ ) {
    if ( pfCandidateCollection.at(i).particleId() == reco::PFCandidate::mu ) pfMuonCollection.push_back( pfCandidateCollection.at(i) );
  }
  // Make PF collections and compute PF MET
  reco::PFCandidateCollection pfChargedEMCollection;
  reco::PFCandidateCollection pfChargedHadronCollection;
  reco::PFCandidateCollection pfChargedHadronPUCollection;
  reco::PFCandidateCollection pfNeutralEMCollection;
  reco::PFCandidateCollection pfNeutralHadronCollection;
  TLorentzVector pfMETP4 = TLorentzVector();
  for (uint i = 0; i < pfCandidateCollection.size(); i++ ) {
    const reco::PFCandidate pfCandidate = pfCandidateCollection.at(i);
    // Compute PF MET
    TLorentzVector pfP4 = TLorentzVector();
    pfP4.SetPtEtaPhiM( pfCandidate.pt() , pfCandidate.eta(), pfCandidate.phi(), pfCandidate.mass() );
    pfMETP4 -= pfP4;
    // Determine if Pile-Up
    bool isPFPU = false;
    PFPileUpAlgo puAlgo;
    puAlgo.setCheckClosestZVertex(true);
    int ivertex = puAlgo.chargedHadronVertex( priVtxCollection, pfCandidate );
    if( ivertex!=-1 && ivertex!=0 ) { isPFPU = true; }
    bool keepPFCandidate = false;
    // Fill PF Charged EM collection
    if ( pfCandidate.particleId()==reco::PFCandidate::e || pfCandidate.particleId()==reco::PFCandidate::mu ) {
      pfChargedEMCollection.push_back( pfCandidate ); keepPFCandidate = true;
    }
    // Fill PF Charged Hadron collection
    if ( pfCandidate.particleId()==reco::PFCandidate::h && isPFPU==false ) {
      pfChargedHadronCollection.push_back( pfCandidate ); keepPFCandidate = true;
    }
    // Fill PF Charged Hadron from PU collection
    if ( pfCandidate.particleId()==reco::PFCandidate::h && isPFPU==true ) {
      if (pfCandidate.pt() > 0.5 ) { pfChargedHadronPUCollection.push_back( pfCandidate ); keepPFCandidate = true; }
    }
    // Fill PF Neutral EM collection
    if ( pfCandidate.particleId()==reco::PFCandidate::gamma || pfCandidate.particleId()==reco::PFCandidate::egamma_HF ) {
      if (pfCandidate.pt() > 0.5 ) { pfNeutralEMCollection.push_back( pfCandidate ); keepPFCandidate = true; }
    }
    // Fill PF Neutral Hadron collection
    if ( pfCandidate.particleId()==reco::PFCandidate::h0 || pfCandidate.particleId()==reco::PFCandidate::h_HF ) {
      if (pfCandidate.pt() > 0.5 ) { pfNeutralHadronCollection.push_back( pfCandidate );  keepPFCandidate = true; }
    }
    if (keepPFCandidate && pfMuonCollection.size()>0) {
      this->PF_Candidate_isPU.push_back ( isPFPU );
      this->PF_Candidate_Id.push_back   ( pfCandidate.particleId() );
      this->PF_Candidate_Eta.push_back  ( pfCandidate.eta()        );
      this->PF_Candidate_Phi.push_back  ( pfCandidate.phi()        );
      this->PF_Candidate_Pt.push_back   ( (pfCandidate.charge()!=0) ? pfCandidate.pt() : pfCandidate.et() );
    }
  }
  pfMETP4.SetPtEtaPhiM( pfMETP4.Pt(), 0.0, pfMETP4.Phi(), 0.0 );
  this->PF_MET_P2.Set( pfMETP4.Px() , pfMETP4.Py() );
  // PF Muon
  this->PF_Muon_N = pfMuonCollection.size();
  for (ushort imuon = 0, idimuon = 0; imuon < pfMuonCollection.size(); imuon++) { 
    const reco::PFCandidate& pfMuon = pfMuonCollection.at(imuon);
    // PF Muon Information
    TLorentzVector muP4 = TLorentzVector();
    muP4.SetPtEtaPhiM( pfMuon.pt() , pfMuon.eta() , pfMuon.phi() , pfMuon.mass() );
    new ((*this->PF_Muon_P4)[imuon]) TLorentzVector( muP4 );
    this->PF_Muon_Charge.push_back     ( pfMuon.charge()            );
    this->PF_Muon_Gen_Index.push_back  ( indexMap.at("GEN")[imuon]  );
    this->PF_Muon_Reco_Index.push_back ( indexMap.at("RECO")[imuon] );
    // Compute PF Isolation
    Float_t sumChargedEMPtR03       = pfIsolation ( pfMuon, pfChargedEMCollection,       0.3, 0.0001 );
    Float_t sumChargedHadronPtR03   = pfIsolation ( pfMuon, pfChargedHadronCollection,   0.3, 0.0001 );
    Float_t sumChargedHadronPUPtR03 = pfIsolation ( pfMuon, pfChargedHadronPUCollection, 0.3, 0.01   );
    Float_t sumNeutralEMEtR03       = pfIsolation ( pfMuon, pfNeutralEMCollection,       0.3, 0.01   );
    Float_t sumNeutralHadronEtR03   = pfIsolation ( pfMuon, pfNeutralHadronCollection,   0.3, 0.01   );
    Float_t isoR03PUCorr   = sumChargedHadronPtR03 + std::max(0., sumNeutralHadronEtR03 + sumNeutralEMEtR03 - 0.5*sumChargedHadronPUPtR03 );
    Float_t isoR03NoPUCorr = sumChargedHadronPtR03 + std::max(0.0f, sumNeutralHadronEtR03 + sumNeutralEMEtR03 );
    this->PF_Muon_PFIsoR03_ChargedEM_SumPt.push_back    ( sumChargedEMPtR03                    );
    this->PF_Muon_PFIsoR03_NeutralEM_SumEt.push_back    ( sumNeutralEMEtR03                    );
    this->PF_Muon_PFIsoR03_ChargedHad_SumPt.push_back   ( sumChargedHadronPtR03                );
    this->PF_Muon_PFIsoR03_NeutralHad_SumEt.push_back   ( sumNeutralHadronEtR03                );
    this->PF_Muon_PFIsoR03_ChargedHadPU_SumPt.push_back ( sumChargedHadronPUPtR03              );
    this->PF_Muon_PFIsoR03_IsoPUCorr.push_back          ( isoR03PUCorr   / (pfMuon.pt()+1E-12) );
    this->PF_Muon_PFIsoR03_IsoNoPUCorr.push_back        ( isoR03NoPUCorr / (pfMuon.pt()+1E-12) );
    Float_t sumChargedEMPtR04       = pfIsolation ( pfMuon, pfChargedEMCollection,       0.4, 0.0001 );
    Float_t sumChargedHadronPtR04   = pfIsolation ( pfMuon, pfChargedHadronCollection,   0.4, 0.0001 );
    Float_t sumChargedHadronPUPtR04 = pfIsolation ( pfMuon, pfChargedHadronPUCollection, 0.4, 0.01   );
    Float_t sumNeutralEMEtR04       = pfIsolation ( pfMuon, pfNeutralEMCollection,       0.4, 0.01   );
    Float_t sumNeutralHadronEtR04   = pfIsolation ( pfMuon, pfNeutralHadronCollection,   0.4, 0.01   );
    Float_t isoR04PUCorr   = sumChargedHadronPtR04 + std::max(0., sumNeutralHadronEtR04 + sumNeutralEMEtR04 - 0.5*sumChargedHadronPUPtR04 );
    Float_t isoR04NoPUCorr = sumChargedHadronPtR04 + std::max(0.0f, sumNeutralHadronEtR04 + sumNeutralEMEtR04 );
    this->PF_Muon_PFIsoR04_ChargedEM_SumPt.push_back    ( sumChargedEMPtR04                    );
    this->PF_Muon_PFIsoR04_NeutralEM_SumEt.push_back    ( sumNeutralEMEtR04                    );
    this->PF_Muon_PFIsoR04_ChargedHad_SumPt.push_back   ( sumChargedHadronPtR04                );
    this->PF_Muon_PFIsoR04_NeutralHad_SumEt.push_back   ( sumNeutralHadronEtR04                );
    this->PF_Muon_PFIsoR04_ChargedHadPU_SumPt.push_back ( sumChargedHadronPUPtR04              );
    this->PF_Muon_PFIsoR04_IsoPUCorr.push_back          ( isoR04PUCorr   / (pfMuon.pt()+1E-12) );
    this->PF_Muon_PFIsoR04_IsoNoPUCorr.push_back        ( isoR04NoPUCorr / (pfMuon.pt()+1E-12) );
    // Comput Muon-MET transverse momenta
    TLorentzVector pfMuonP4T = TLorentzVector();
    pfMuonP4T.SetPtEtaPhiM(pfMuon.pt(), 0.0, pfMuon.phi(), pfMuon.mass());
    new ((*this->PF_MuonMET_P4T)[imuon]) TLorentzVector( pfMuonP4T + pfMETP4 );
    // PF Dimuon
    for (ushort imuon2 = imuon+1; imuon2 < pfMuonCollection.size(); imuon2++) {
      const reco::PFCandidate& pfMuon2 = pfMuonCollection.at(imuon2);
      // PF Dimuon Information
      TLorentzVector mu2P4 = TLorentzVector();
      mu2P4.SetPtEtaPhiM( pfMuon2.pt() , pfMuon2.eta() , pfMuon2.phi() , pfMuon2.mass() );
      TLorentzVector diMuonP4 = muP4 + mu2P4;
      Char_t charge = pfMuon.charge() + pfMuon2.charge();
      TVector3 diMuonVtx = TVector3();
      Float_t vProb    = -1. ;
      Float_t dca      = -1.;
      Float_t massWErr = -1.;
      if ( pfMuon.trackRef().isNonnull() && pfMuon2.trackRef().isNonnull() ) {
        std::vector<reco::TransientTrack> t_tks;
        t_tks.push_back( _theTTBuilder->build( *pfMuon.trackRef()  ) );
        t_tks.push_back( _theTTBuilder->build( *pfMuon2.trackRef() ) );
        KalmanVertexFitter  vtxFitter = KalmanVertexFitter(true);
        TransientVertex diMuonVertex = vtxFitter.vertex( t_tks );
        if (diMuonVertex.isValid()) {
          Float_t vChi2 = diMuonVertex.totalChiSquared();
          Float_t vNDF  = diMuonVertex.degreesOfFreedom();
          vProb = TMath::Prob(vChi2,(int)vNDF);
          diMuonVtx.SetXYZ( diMuonVertex.position().x() , diMuonVertex.position().y() , diMuonVertex.position().z() );
        }
        TrajectoryStateClosestToPoint mu1TS = t_tks[0].impactPointTSCP();
        TrajectoryStateClosestToPoint mu2TS = t_tks[1].impactPointTSCP();
        if (mu1TS.isValid() && mu2TS.isValid()) {
          ClosestApproachInRPhi cApp;
          cApp.calculate(mu1TS.theState(), mu2TS.theState());
          if ( cApp.status() ) dca = cApp.distance();
        }
        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );
        InvariantMassFromVertex massCalculator;
        massWErr = massCalculator.invariantMass( VtxForInvMass, _muMasses ).error();
      }
      this->PF_DiMuon_N = idimuon + 1;
      new ((*this->PF_DiMuon_P4)[idimuon]) TLorentzVector( diMuonP4 );
      this->PF_DiMuon_Charge.push_back      ( charge    );
      new ((*this->PF_DiMuon_Vertex)[idimuon]) TVector3( diMuonVtx );
      this->PF_DiMuon_VertexProb.push_back  ( vProb     );
      this->PF_DiMuon_DCA.push_back         ( dca       );
      this->PF_DiMuon_MassError.push_back   ( massWErr  );
      this->PF_DiMuon_Muon1_Index.push_back ( imuon     );
      this->PF_DiMuon_Muon2_Index.push_back ( imuon2    );
      idimuon++;
    }
  }
}

//--------------------------------------------------------------------------------------------------
short 
HiMuonEvent::findGenIndex(const reco::GenParticleRef& genMatchRef, const reco::GenParticleRefVector& genRefCollection)
{
  if ( genMatchRef.isNonnull() && genMatchRef.isAvailable() ) {
    for (ushort ipar = 0; ipar < genRefCollection.size(); ipar++ ) {
      const reco::GenParticleRef& genParticleRef = genRefCollection.at(ipar);
      if (genMatchRef == genParticleRef) return ipar;
    }
  }
  return -1;
}

//--------------------------------------------------------------------------------------------------
void
HiMuonEvent::Fill(const reco::GenParticleRefVector& genRefCollection, const IndexMap& indexMap)
{
  // Search for missing daughters and add them to the gen collection
  reco::GenParticleRefVector genParticleRefCollection(genRefCollection);
  for (uint ipar = 0; ipar < genParticleRefCollection.size(); ipar++ ) {
    const reco::GenParticleRef& genParticleRef = genParticleRefCollection.at(ipar);
    uint numberOfDaughters = genParticleRef->numberOfDaughters();
    for (uint idau = 0; idau < numberOfDaughters; idau++ ) {
      reco::GenParticleRef genDaughterRef = genParticleRef->daughterRef(idau);
      if ( genDaughterRef.isNonnull() && genDaughterRef.isAvailable() && (findGenIndex(genDaughterRef, genParticleRefCollection) == -1) ) { 
        genParticleRefCollection.push_back( genDaughterRef ); 
      }
    }
  }
  // Loop over the Gen Particle Collection
  std::vector< UChar_t > genMuonIndex;
  reco::GenParticleCollection genMuonCollection;
  for (uint ipar = 0; ipar < genParticleRefCollection.size(); ipar++ ) {
    const reco::GenParticle& genParticle = *(genParticleRefCollection.at(ipar));
    // Fill the Gen Muon Collection
    if ( (abs(genParticle.pdgId()) == 13) && (genParticle.status() == 1) ) { genMuonCollection.push_back( genParticle ); genMuonIndex.push_back( ipar ); }
    // Fill the Gen Particle Variables
    TLorentzVector p4 = TLorentzVector();
    p4.SetPtEtaPhiM( genParticle.pt() , genParticle.eta() , genParticle.phi() , genParticle.mass() );
    new ((*this->Gen_Particle_P4)[ipar]) TLorentzVector( p4 );
    this->Gen_Particle_Status.push_back ( genParticle.status() );
    this->Gen_Particle_PdgId.push_back  ( genParticle.pdgId()  );
    std::vector< ushort > motherIndex;
    for (ushort imom = 0; imom < genParticle.numberOfMothers(); imom++ ) {
      short index = findGenIndex(genParticle.motherRef(imom), genParticleRefCollection);
      if (index >= 0) motherIndex.push_back( index ); 
    }
    this->Gen_Particle_Mother_Index.push_back( motherIndex );
    std::vector< ushort > daughterIndex;
    for (ushort idau = 0; idau < genParticle.numberOfDaughters(); idau++ ) {
      short index = findGenIndex(genParticle.daughterRef(idau), genParticleRefCollection);
      if (index >= 0) daughterIndex.push_back( index );
      if (index <  0)  std::cout << "[Error] Daughter not found: " << genParticle.daughterRef(idau)->pdgId() << " pt: " << genParticle.daughterRef(idau)->pt() << std::endl;
    }
    this->Gen_Particle_Daughter_Index.push_back( daughterIndex );
  }
  // Gen Muon
  this->Gen_Muon_N = genMuonCollection.size();
  for (ushort imuon = 0; imuon < genMuonCollection.size(); imuon++ ) {
    const reco::GenParticle& genMuon = genMuonCollection.at(imuon);
    // Gen Muon Information
    TLorentzVector muP4 = TLorentzVector();
    muP4.SetPtEtaPhiM( genMuon.pt() , genMuon.eta() , genMuon.phi() , genMuon.mass() );
    new ((*this->Gen_Muon_P4)[imuon]) TLorentzVector( muP4 );
    this->Gen_Muon_Charge.push_back     ( genMuon.charge()           );
    this->Gen_Muon_Particle_Index.push_back ( genMuonIndex[imuon]    );
    this->Gen_Muon_Reco_Index.push_back ( indexMap.at("RECO")[imuon] );
    this->Gen_Muon_PF_Index.push_back   ( indexMap.at("PF")[imuon]   );
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiMuonEvent::IniArrays()
{
  // Initialize the TCloneArrays
  this->Reco_Muon_P4               = new TClonesArray("TLorentzVector", 100);
  this->Reco_Muon_InnerTrack_P3    = new TClonesArray("TVector3",       100);
  this->Reco_Muon_GlobalTrack_P3   = new TClonesArray("TVector3",       100);
  this->Reco_Muon_BestTrack_P3     = new TClonesArray("TVector3",       100);
  this->Reco_Muon_BestTrack_Vertex = new TClonesArray("TVector3",       100);
  this->Reco_DiMuon_P4             = new TClonesArray("TLorentzVector", 1000);
  this->Reco_DiMuon_Vertex         = new TClonesArray("TVector3",       1000);
  this->PF_Muon_P4                 = new TClonesArray("TLorentzVector", 100);
  this->PF_DiMuon_P4               = new TClonesArray("TLorentzVector", 1000);
  this->PF_DiMuon_Vertex           = new TClonesArray("TVector3",       1000);
  this->PF_MuonMET_P4T             = new TClonesArray("TLorentzVector", 100);
  this->Gen_Particle_P4            = new TClonesArray("TLorentzVector", 1000);
  this->Gen_Muon_P4                = new TClonesArray("TLorentzVector", 100);
}

//--------------------------------------------------------------------------------------------------
void 
HiMuonEvent::SetBranches(const std::string name, const StringBoolMap& doMuon)
{
  // Set Branches
  if ( doMuon.at("Event") && (name == "Event" || name == "All") ) {
    // Event Info
    tree_->Branch("Event_Run",                        &(this->Event_nRun),           "Event_Run/i"    );
    tree_->Branch("Event_Lumi",                       &(this->Event_nLumi),          "Event_Lumi/s"   );
    tree_->Branch("Event_Bx",                         &(this->Event_nBX),            "Event_Bx/i"     );
    tree_->Branch("Event_Orbit",                      &(this->Event_nOrbit),         "Event_Orbit/l"  );
    tree_->Branch("Event_Number",                     &(this->Event_nEvent),         "Event_Number/l" );
    tree_->Branch("Event_nPV",                        &(this->Event_nPV),            "Event_nPV/b"    );
    tree_->Branch("Event_PriVtx_Pos",                 "TVector3",      &(this->Event_PriVtx_Position) );
    tree_->Branch("Event_PriVtx_Err",                 "TVector3",      &(this->Event_PriVtx_Error   ) );
    // Trigger Info
    tree_->Branch("Event_Trig_Fired",                 &(this->Event_Trigger_isFired)                  );
    tree_->Branch("Event_Trig_Presc",                 &(this->Event_Trigger_Prescale)                 );
  }
  if ( doMuon.at("Reco") && (name == "Reco" || name == "All") ) {
    // Reco Muon Kinematic
    tree_->Branch("Reco_Muon_N",                      &(this->Reco_Muon_N),          "Reco_Muon_N/b"   );
    tree_->Branch("Reco_Muon_Mom",                    "TClonesArray",   &(this->Reco_Muon_P4), 32000, 0 );
    tree_->Branch("Reco_Muon_Charge",                 &(this->Reco_Muon_Charge)                       );
    // Reco Muon Matched Index
    if (doMuon.at("Gen")) tree_->Branch("Reco_Muon_Gen_Idx",             &(this->Reco_Muon_Gen_Index) );
    if (doMuon.at("PF") ) tree_->Branch("Reco_Muon_PF_Idx",              &(this->Reco_Muon_PF_Index)  );
    if (doMuon.at("Pat")) {
      // PAT Muon Info
      tree_->Branch("Pat_Muon_Trig",                  &(this->Pat_Muon_TriggerMatched)                );
      tree_->Branch("Pat_Muon_dB",                    &(this->Pat_Muon_dB)                            );
      tree_->Branch("Pat_Muon_dBErr",                 &(this->Pat_Muon_dBErr)                         );
    }
    // Reco Muon ID Flags
    tree_->Branch("Reco_Muon_isPF",                   &(this->Reco_Muon_isPFMuon)                     );
    tree_->Branch("Reco_Muon_isGlobal",               &(this->Reco_Muon_isGlobalMuon)                 );
    tree_->Branch("Reco_Muon_isTracker",              &(this->Reco_Muon_isTrackerMuon)                );
    tree_->Branch("Reco_Muon_isStandAlone",           &(this->Reco_Muon_isStandAloneMuon)             );
    tree_->Branch("Reco_Muon_isLoose",                &(this->Reco_Muon_isLooseMuon)                  );
    tree_->Branch("Reco_Muon_isMedium",               &(this->Reco_Muon_isMediumMuon)                 );
    tree_->Branch("Reco_Muon_isHighPt",               &(this->Reco_Muon_isHighPtMuon)                 );
    tree_->Branch("Reco_Muon_isSoft",                 &(this->Reco_Muon_isSoftMuon)                   );
    tree_->Branch("Reco_Muon_isTight",                &(this->Reco_Muon_isTightMuon)                  );
    tree_->Branch("Reco_Muon_isArbitrated",           &(this->Reco_Muon_isArbitrated)                 );
    tree_->Branch("Reco_Muon_TrackerArbitrated",      &(this->Reco_Muon_TrackerMuonArbitrated)        );
    tree_->Branch("Reco_Muon_GlobalPromptTight",      &(this->Reco_Muon_GlobalMuonPromptTight)        );
    tree_->Branch("Reco_Muon_TMLastStationLoose",     &(this->Reco_Muon_TMLastStationLoose)           );
    tree_->Branch("Reco_Muon_TMLastStationTight",     &(this->Reco_Muon_TMLastStationTight)           );
    tree_->Branch("Reco_Muon_TM2DCompatibilityLoose", &(this->Reco_Muon_TM2DCompatibilityLoose)       );
    tree_->Branch("Reco_Muon_TM2DCompatibilityTight", &(this->Reco_Muon_TM2DCompatibilityTight)       );
    tree_->Branch("Reco_Muon_TMOneStationLoose",      &(this->Reco_Muon_TMOneStationLoose)            );
    tree_->Branch("Reco_Muon_TMOneStationTight",      &(this->Reco_Muon_TMOneStationTight)            );
    tree_->Branch("Reco_Muon_GMTkChiCompatibility",   &(this->Reco_Muon_GMTkChiCompatibility)         );
    tree_->Branch("Reco_Muon_GMStaChiCompatibility",  &(this->Reco_Muon_GMStaChiCompatibility)        );
    tree_->Branch("Reco_Muon_GMTkKinkTight",          &(this->Reco_Muon_GMTkKinkTight)                );
    tree_->Branch("Reco_Muon_TMLastStationAngLoose",  &(this->Reco_Muon_TMLastStationAngLoose)        );
    tree_->Branch("Reco_Muon_TMLastStationAngTight",  &(this->Reco_Muon_TMLastStationAngTight)        );
    tree_->Branch("Reco_Muon_TMOneStationAngLoose",   &(this->Reco_Muon_TMOneStationAngLoose)         );
    tree_->Branch("Reco_Muon_TMOneStationAngTight",   &(this->Reco_Muon_TMOneStationAngTight)         );
    // Reco Muon ID vars
    tree_->Branch("Reco_Muon_MatchedStations",        &(this->Reco_Muon_NumberOfMatchedStations)      );
    tree_->Branch("Reco_Muon_Matches",                &(this->Reco_Muon_NumberOfMatches)              );
    tree_->Branch("Reco_Muon_SegmentComp",            &(this->Reco_Muon_SegmentCompatibility)         );
    tree_->Branch("Reco_Muon_Chi2Pos",                &(this->Reco_Muon_Chi2LocalPosition)            );
    tree_->Branch("Reco_Muon_TrkKink",                &(this->Reco_Muon_TrkKink)                      );
    // Reco Muon Inner Track
    tree_->Branch("Reco_Muon_InTrk_Mom",          "TClonesArray", &(this->Reco_Muon_InnerTrack_P3), 32000, 0 );
    tree_->Branch("Reco_Muon_InTrk_PtErr",            &(this->Reco_Muon_InnerTrack_PtError)           );
    tree_->Branch("Reco_Muon_InTrk_isHighPurity",     &(this->Reco_Muon_InnerTrack_isHighPurity)      );
    tree_->Branch("Reco_Muon_InTrk_ValidHits",        &(this->Reco_Muon_InnerTrack_NumOfValHits)      );
    tree_->Branch("Reco_Muon_InTrk_LostHits",         &(this->Reco_Muon_InnerTrack_NumOfLostHits)     );
    tree_->Branch("Reco_Muon_InTrk_ValidPixHits",     &(this->Reco_Muon_InnerTrack_NumOfValPixHits)   );
    tree_->Branch("Reco_Muon_InTrk_TrkLayers",        &(this->Reco_Muon_InnerTrack_TrkLayersWithMea)  );
    tree_->Branch("Reco_Muon_InTrk_PixLayers",        &(this->Reco_Muon_InnerTrack_PixLayersWithMea)  );
    tree_->Branch("Reco_Muon_InTrk_dXY",              &(this->Reco_Muon_InnerTrack_dXY)               );
    tree_->Branch("Reco_Muon_InTrk_dXYErr",           &(this->Reco_Muon_InnerTrack_dXYErr)            );
    tree_->Branch("Reco_Muon_InTrk_dZ",               &(this->Reco_Muon_InnerTrack_dZ)                );
    tree_->Branch("Reco_Muon_InTrk_dZErr",            &(this->Reco_Muon_InnerTrack_dZErr)             );
    tree_->Branch("Reco_Muon_InTrk_ValFrac",          &(this->Reco_Muon_InnerTrack_ValFrac)           );
    tree_->Branch("Reco_Muon_InTrk_NormChi2",         &(this->Reco_Muon_InnerTrack_NormChi2)          );
    // Reco Muon Global Track
    tree_->Branch("Reco_Muon_GlbTrk_Mom",         "TClonesArray", &(this->Reco_Muon_GlobalTrack_P3), 32000, 0 );
    tree_->Branch("Reco_Muon_GlbTrk_PtErr",           &(this->Reco_Muon_GlobalTrack_PtError)          );
    tree_->Branch("Reco_Muon_GlbTrk_ValidMuonHits",   &(this->Reco_Muon_GlobalTrack_NumOfValMuonHits) );
    tree_->Branch("Reco_Muon_GlbTrk_NormChi2",        &(this->Reco_Muon_GlobalTrack_NormChi2)         );
    // Reco Muon Best Track
    tree_->Branch("Reco_Muon_BestTrk_Type",           &(this->Reco_Muon_BestTrack_Type)               );
    tree_->Branch("Reco_Muon_BestTrk_Mom",        "TClonesArray", &(this->Reco_Muon_BestTrack_P3), 32000, 0 );
    tree_->Branch("Reco_Muon_BestTrk_Vertex",     "TClonesArray", &(this->Reco_Muon_BestTrack_Vertex), 32000, 0 );
    tree_->Branch("Reco_Muon_BestTrk_PtErr",          &(this->Reco_Muon_BestTrack_PtError)            );
    tree_->Branch("Reco_Muon_BestTrk_dXY",            &(this->Reco_Muon_BestTrack_dXY)                );
    tree_->Branch("Reco_Muon_BestTrk_dXYErr",         &(this->Reco_Muon_BestTrack_dXYErr)             );
    tree_->Branch("Reco_Muon_BestTrk_dZ",             &(this->Reco_Muon_BestTrack_dZ)                 );
    tree_->Branch("Reco_Muon_BestTrk_dZErr",          &(this->Reco_Muon_BestTrack_dZErr)              );
    // Reco Muon Isolation
    tree_->Branch("Reco_Muon_IsoPFR03",               &(this->Reco_Muon_PFIsoR03_IsoPUCorr)           );
    tree_->Branch("Reco_Muon_IsoPFR03NoPUCorr",       &(this->Reco_Muon_PFIsoR03_IsoNoPUCorr)         );
    tree_->Branch("Reco_Muon_IsoPFR04",               &(this->Reco_Muon_PFIsoR04_IsoPUCorr)           );
    tree_->Branch("Reco_Muon_IsoPFR04NoPUCorr",       &(this->Reco_Muon_PFIsoR04_IsoNoPUCorr)         );
    tree_->Branch("Reco_Muon_EM_Chg_sumR03Pt",        &(this->Reco_Muon_PFIsoR03_ChargedEM_SumPt)     );
    tree_->Branch("Reco_Muon_EM_Chg_sumR04Pt",        &(this->Reco_Muon_PFIsoR04_ChargedEM_SumPt)     );
    tree_->Branch("Reco_Muon_EM_Neu_sumR03Et",        &(this->Reco_Muon_PFIsoR03_NeutralEM_SumEt)     );
    tree_->Branch("Reco_Muon_EM_Neu_sumR04Et",        &(this->Reco_Muon_PFIsoR04_NeutralEM_SumEt)     );
    tree_->Branch("Reco_Muon_Had_Chg_sumR03Pt",       &(this->Reco_Muon_PFIsoR03_ChargedHad_SumPt)    );
    tree_->Branch("Reco_Muon_Had_Chg_sumR04Pt",       &(this->Reco_Muon_PFIsoR04_ChargedHad_SumPt)    );
    tree_->Branch("Reco_Muon_Had_Neu_sumR03Et",       &(this->Reco_Muon_PFIsoR03_NeutralHad_SumEt)    );
    tree_->Branch("Reco_Muon_Had_Neu_sumR04Et",       &(this->Reco_Muon_PFIsoR04_NeutralHad_SumEt)    );
    tree_->Branch("Reco_Muon_Had_PU_sumR03Pt",        &(this->Reco_Muon_PFIsoR03_ChargedHadPU_SumPt)  );
    tree_->Branch("Reco_Muon_Had_PU_sumR04Pt",        &(this->Reco_Muon_PFIsoR04_ChargedHadPU_SumPt)  );
    tree_->Branch("Reco_Muon_IsoR03",                 &(this->Reco_Muon_IsoR03_Iso)                   );
    tree_->Branch("Reco_Muon_IsoR05",                 &(this->Reco_Muon_IsoR05_Iso)                   );
    tree_->Branch("Reco_Muon_Trk_sumR03Pt",           &(this->Reco_Muon_IsoR03_SumPt)                 );
    tree_->Branch("Reco_Muon_Trk_sumR05Pt",           &(this->Reco_Muon_IsoR05_SumPt)                 );
    // Reco DiMuon
    tree_->Branch("Reco_DiMuon_N",                    &(this->Reco_DiMuon_N),       "Reco_DiMuon_N/s" );
    tree_->Branch("Reco_DiMuon_Mom",                  "TClonesArray", &(this->Reco_DiMuon_P4), 32000, 0 );
    tree_->Branch("Reco_DiMuon_Charge",               &(this->Reco_DiMuon_Charge)                     );
    tree_->Branch("Reco_DiMuon_Muon1_Idx",            &(this->Reco_DiMuon_Muon1_Index)                );
    tree_->Branch("Reco_DiMuon_Muon2_Idx",            &(this->Reco_DiMuon_Muon2_Index)                );
    tree_->Branch("Reco_DiMuon_isCowBoy",             &(this->Reco_DiMuon_isCowBoy)                   );
    tree_->Branch("Reco_DiMuon_Vertex",               "TClonesArray", &(this->Reco_DiMuon_Vertex), 32000, 0 );
    tree_->Branch("Reco_DiMuon_VtxProb",              &(this->Reco_DiMuon_VertexProb)                 );
    tree_->Branch("Reco_DiMuon_DCA",                  &(this->Reco_DiMuon_DCA)                        );
    tree_->Branch("Reco_DiMuon_MassErr",              &(this->Reco_DiMuon_MassError)                  );
  }  
  if ( doMuon.at("PF") && (name == "PF" || name == "All") ) {
    // PF Candidate
    tree_->Branch("PF_Candidate_isPU",                &(this->PF_Candidate_isPU)                      );
    tree_->Branch("PF_Candidate_Id",                  &(this->PF_Candidate_Id)                        );
    tree_->Branch("PF_Candidate_Eta",                 &(this->PF_Candidate_Eta)                       );
    tree_->Branch("PF_Candidate_Phi",                 &(this->PF_Candidate_Phi)                       );
    tree_->Branch("PF_Candidate_Pt",                  &(this->PF_Candidate_Pt)                        );
    // PF Muon
    tree_->Branch("PF_Muon_N",                        &(this->PF_Muon_N),           "PF_Muon_N/b"     );
    tree_->Branch("PF_Muon_Mom",                      "TClonesArray", &(this->PF_Muon_P4), 32000, 0   );
    tree_->Branch("PF_Muon_Charge",                   &(this->PF_Muon_Charge)                         );
    if (doMuon.at("Gen") ) tree_->Branch("PF_Muon_Gen_Idx",               &(this->PF_Muon_Gen_Index)  );
    if (doMuon.at("Reco")) tree_->Branch("PF_Muon_Reco_Idx",              &(this->PF_Muon_Reco_Index) );
    // Reco Muon Isolation
    tree_->Branch("PF_Muon_IsoPFR03",                 &(this->PF_Muon_PFIsoR03_IsoPUCorr)             );
    tree_->Branch("PF_Muon_IsoPFR03NoPUCorr",         &(this->PF_Muon_PFIsoR03_IsoNoPUCorr)           );
    tree_->Branch("PF_Muon_IsoPFR04",                 &(this->PF_Muon_PFIsoR04_IsoPUCorr)             );
    tree_->Branch("PF_Muon_IsoPFR04NoPUCorr",         &(this->PF_Muon_PFIsoR04_IsoNoPUCorr)           );
    tree_->Branch("PF_Muon_EM_Chg_sumR03Pt",          &(this->PF_Muon_PFIsoR03_ChargedEM_SumPt)       );
    tree_->Branch("PF_Muon_EM_Chg_sumR04Pt",          &(this->PF_Muon_PFIsoR04_ChargedEM_SumPt)       );
    tree_->Branch("PF_Muon_EM_Neu_sumR03Et",          &(this->PF_Muon_PFIsoR03_NeutralEM_SumEt)       );
    tree_->Branch("PF_Muon_EM_Neu_sumR04Et",          &(this->PF_Muon_PFIsoR04_NeutralEM_SumEt)       );
    tree_->Branch("PF_Muon_Had_Chg_sumR03Pt",         &(this->PF_Muon_PFIsoR03_ChargedHad_SumPt)      );
    tree_->Branch("PF_Muon_Had_Chg_sumR04Pt",         &(this->PF_Muon_PFIsoR04_ChargedHad_SumPt)      );
    tree_->Branch("PF_Muon_Had_Neu_sumR03Et",         &(this->PF_Muon_PFIsoR03_NeutralHad_SumEt)      );
    tree_->Branch("PF_Muon_Had_Neu_sumR04Et",         &(this->PF_Muon_PFIsoR04_NeutralHad_SumEt)      );
    tree_->Branch("PF_Muon_Had_PU_sumR03Pt",          &(this->PF_Muon_PFIsoR03_ChargedHadPU_SumPt)    );
    tree_->Branch("PF_Muon_Had_PU_sumR04Pt",          &(this->PF_Muon_PFIsoR04_ChargedHadPU_SumPt)    );
    // PF DiMuon
    tree_->Branch("PF_DiMuon_N",                      &(this->PF_DiMuon_N),         "PF_DiMuon_N/s"   );
    tree_->Branch("PF_DiMuon_Mom",                    "TClonesArray", &(this->PF_DiMuon_P4), 32000, 0 );
    tree_->Branch("PF_DiMuon_Charge",                 &(this->PF_DiMuon_Charge)                       );
    tree_->Branch("PF_DiMuon_Muon1_Idx",              &(this->PF_DiMuon_Muon1_Index)                  );
    tree_->Branch("PF_DiMuon_Muon2_Idx",              &(this->PF_DiMuon_Muon2_Index)                  );
    tree_->Branch("PF_DiMuon_Vertex",                 "TClonesArray", &(this->PF_DiMuon_Vertex), 32000, 0 );
    tree_->Branch("PF_DiMuon_VtxProb",                &(this->PF_DiMuon_VertexProb)                   );
    tree_->Branch("PF_DiMuon_DCA",                    &(this->PF_DiMuon_DCA)                          );
    tree_->Branch("PF_DiMuon_MassErr",                &(this->PF_DiMuon_MassError)                    );
    // PF MET
    tree_->Branch("PF_MET_Mom",                       "TVector2",             &(this->PF_MET_P2)      );
    // PF Muon-MET
    tree_->Branch("PF_MuonMET_TransMom",              "TClonesArray", &(this->PF_MuonMET_P4T), 32000, 0 );
  }
  if ( doMuon.at("Gen") && (name == "Gen" || name == "All") ) {
    // Gen Particle
    tree_->Branch("Gen_Particle_Mom",                 "TClonesArray", &(this->Gen_Particle_P4), 32000, 0 );
    tree_->Branch("Gen_Particle_PdgId",               &(this->Gen_Particle_PdgId)                     );
    tree_->Branch("Gen_Particle_Status",              &(this->Gen_Particle_Status)                    );
    tree_->Branch("Gen_Particle_Mother_Idx",          &(this->Gen_Particle_Mother_Index)              );
    tree_->Branch("Gen_Particle_Daughter_Idx",        &(this->Gen_Particle_Daughter_Index)            );
    // Gen Muon
    tree_->Branch("Gen_Muon_N",                       &(this->Gen_Muon_N),          "Gen_Muon_N/b"    );
    tree_->Branch("Gen_Muon_Mom",                     "TClonesArray", &(this->Gen_Muon_P4), 32000, 0  );
    tree_->Branch("Gen_Muon_Charge",                  &(this->Gen_Muon_Charge)                        );
    tree_->Branch("Gen_Muon_Particle_Idx",            &(this->Gen_Muon_Particle_Index)                );
    if (doMuon.at("Reco")) tree_->Branch("Gen_Muon_Reco_Idx",         &(this->Gen_Muon_Reco_Index)    );
    if (doMuon.at("PF")  ) tree_->Branch("Gen_Muon_PF_Idx",           &(this->Gen_Muon_PF_Index)      );
  }
}

//--------------------------------------------------------------------------------------------------
void
HiMuonEvent::Clear(void)
{
  // Event Info
  this->Event_nRun   = 0;
  this->Event_nLumi  = 0;
  this->Event_nBX    = 0;
  this->Event_nOrbit = 0;
  this->Event_nEvent = 0;
  this->Event_nPV    = 0;
  this->Event_PriVtx_Position = TVector3();
  this->Event_PriVtx_Error    = TVector3();
  // Trigger Info
  this->Event_Trigger_isFired.clear();
  this->Event_Trigger_Prescale.clear();
  // PAT Muon Info
  this->Pat_Muon_TriggerMatched.clear();
  this->Pat_Muon_dB.clear();
  this->Pat_Muon_dBErr.clear();
  // Reco Muon Kinematic
  this->Reco_Muon_N = 0;
  this->Reco_Muon_P4->Clear();
  this->Reco_Muon_Charge.clear();
  // Reco Muon Matched Index
  this->Reco_Muon_Gen_Index.clear();
  this->Reco_Muon_PF_Index.clear();
  // Reco Muon ID Flags
  this->Reco_Muon_isPFMuon.clear();
  this->Reco_Muon_isGlobalMuon.clear();
  this->Reco_Muon_isTrackerMuon.clear();
  this->Reco_Muon_isStandAloneMuon.clear();
  this->Reco_Muon_isLooseMuon.clear();
  this->Reco_Muon_isMediumMuon.clear();
  this->Reco_Muon_isHighPtMuon.clear();
  this->Reco_Muon_isSoftMuon.clear();
  this->Reco_Muon_isTightMuon.clear();
  this->Reco_Muon_isArbitrated.clear();
  this->Reco_Muon_TrackerMuonArbitrated.clear();
  this->Reco_Muon_GlobalMuonPromptTight.clear();
  this->Reco_Muon_TMLastStationLoose.clear();
  this->Reco_Muon_TMLastStationTight.clear();
  this->Reco_Muon_TM2DCompatibilityLoose.clear();
  this->Reco_Muon_TM2DCompatibilityTight.clear();
  this->Reco_Muon_TMOneStationLoose.clear();
  this->Reco_Muon_TMOneStationTight.clear();
  this->Reco_Muon_GMTkChiCompatibility.clear();
  this->Reco_Muon_GMStaChiCompatibility.clear();
  this->Reco_Muon_GMTkKinkTight.clear();
  this->Reco_Muon_TMLastStationAngLoose.clear();
  this->Reco_Muon_TMLastStationAngTight.clear();
  this->Reco_Muon_TMOneStationAngLoose.clear();
  this->Reco_Muon_TMOneStationAngTight.clear();
  // Reco Muon ID vars
  this->Reco_Muon_NumberOfMatchedStations.clear();
  this->Reco_Muon_NumberOfMatches.clear();
  this->Reco_Muon_SegmentCompatibility.clear();
  this->Reco_Muon_Chi2LocalPosition.clear();
  this->Reco_Muon_TrkKink.clear();
  // Reco Muon Inner Track
  this->Reco_Muon_InnerTrack_P3->Clear();
  this->Reco_Muon_InnerTrack_PtError.clear();
  this->Reco_Muon_InnerTrack_isHighPurity.clear();
  this->Reco_Muon_InnerTrack_NumOfValHits.clear();
  this->Reco_Muon_InnerTrack_NumOfLostHits.clear();
  this->Reco_Muon_InnerTrack_NumOfValPixHits.clear();
  this->Reco_Muon_InnerTrack_TrkLayersWithMea.clear();
  this->Reco_Muon_InnerTrack_PixLayersWithMea.clear();
  this->Reco_Muon_InnerTrack_dXY.clear();
  this->Reco_Muon_InnerTrack_dXYErr.clear();
  this->Reco_Muon_InnerTrack_dZ.clear();
  this->Reco_Muon_InnerTrack_dZErr.clear();
  this->Reco_Muon_InnerTrack_ValFrac.clear();
  this->Reco_Muon_InnerTrack_NormChi2.clear();
  // Reco Muon Global Track
  this->Reco_Muon_GlobalTrack_P3->Clear();
  this->Reco_Muon_GlobalTrack_PtError.clear();
  this->Reco_Muon_GlobalTrack_NumOfValMuonHits.clear();
  this->Reco_Muon_GlobalTrack_NormChi2.clear();
  // Reco Muon Best Track
  this->Reco_Muon_BestTrack_P3->Clear();
  this->Reco_Muon_BestTrack_Vertex->Clear();
  this->Reco_Muon_BestTrack_Type.clear();
  this->Reco_Muon_BestTrack_PtError.clear();
  this->Reco_Muon_BestTrack_dXY.clear();
  this->Reco_Muon_BestTrack_dXYErr.clear();
  this->Reco_Muon_BestTrack_dZ.clear();
  this->Reco_Muon_BestTrack_dZErr.clear();
  // Reco Muon Isolation
  this->Reco_Muon_PFIsoR03_ChargedEM_SumPt.clear();
  this->Reco_Muon_PFIsoR03_NeutralEM_SumEt.clear();
  this->Reco_Muon_PFIsoR03_ChargedHad_SumPt.clear();
  this->Reco_Muon_PFIsoR03_NeutralHad_SumEt.clear();
  this->Reco_Muon_PFIsoR03_ChargedHadPU_SumPt.clear();
  this->Reco_Muon_PFIsoR03_IsoPUCorr.clear();
  this->Reco_Muon_PFIsoR03_IsoNoPUCorr.clear();
  this->Reco_Muon_PFIsoR04_ChargedEM_SumPt.clear();
  this->Reco_Muon_PFIsoR04_ChargedHad_SumPt.clear();
  this->Reco_Muon_PFIsoR04_NeutralHad_SumEt.clear();
  this->Reco_Muon_PFIsoR04_NeutralEM_SumEt.clear();
  this->Reco_Muon_PFIsoR04_ChargedHadPU_SumPt.clear();
  this->Reco_Muon_PFIsoR04_IsoPUCorr.clear();
  this->Reco_Muon_PFIsoR04_IsoNoPUCorr.clear();
  this->Reco_Muon_IsoR03_SumPt.clear();
  this->Reco_Muon_IsoR03_Iso.clear();
  this->Reco_Muon_IsoR05_SumPt.clear();
  this->Reco_Muon_IsoR05_Iso.clear();
  // Reco DiMuon
  this->Reco_DiMuon_N = 0;
  this->Reco_DiMuon_P4->Clear();
  this->Reco_DiMuon_Charge.clear();
  this->Reco_DiMuon_Muon1_Index.clear();
  this->Reco_DiMuon_Muon2_Index.clear();
  this->Reco_DiMuon_isCowBoy.clear();
  this->Reco_DiMuon_Vertex->Clear();
  this->Reco_DiMuon_VertexProb.clear();
  this->Reco_DiMuon_DCA.clear();
  this->Reco_DiMuon_MassError.clear();
  // PF Candidate
  this->PF_Candidate_isPU.clear();
  this->PF_Candidate_Id.clear();
  this->PF_Candidate_Eta.clear();
  this->PF_Candidate_Phi.clear();
  this->PF_Candidate_Pt.clear();
  // PF Muon
  this->PF_Muon_N = 0;
  this->PF_Muon_P4->Clear();
  this->PF_Muon_Charge.clear();
  this->PF_Muon_Gen_Index.clear();
  this->PF_Muon_Reco_Index.clear();
  // Reco Muon Isolation
  this->PF_Muon_PFIsoR03_ChargedEM_SumPt.clear();
  this->PF_Muon_PFIsoR03_NeutralEM_SumEt.clear();
  this->PF_Muon_PFIsoR03_ChargedHad_SumPt.clear();
  this->PF_Muon_PFIsoR03_NeutralHad_SumEt.clear();
  this->PF_Muon_PFIsoR03_ChargedHadPU_SumPt.clear();
  this->PF_Muon_PFIsoR03_IsoPUCorr.clear();
  this->PF_Muon_PFIsoR03_IsoNoPUCorr.clear();
  this->PF_Muon_PFIsoR04_ChargedEM_SumPt.clear();
  this->PF_Muon_PFIsoR04_ChargedHad_SumPt.clear();
  this->PF_Muon_PFIsoR04_NeutralHad_SumEt.clear();
  this->PF_Muon_PFIsoR04_NeutralEM_SumEt.clear();
  this->PF_Muon_PFIsoR04_ChargedHadPU_SumPt.clear();
  this->PF_Muon_PFIsoR04_IsoPUCorr.clear();
  this->PF_Muon_PFIsoR04_IsoNoPUCorr.clear();
  // PF DiMuon
  this->PF_DiMuon_N = 0;
  this->PF_DiMuon_P4->Clear();
  this->PF_DiMuon_Charge.clear();
  this->PF_DiMuon_Muon1_Index.clear();
  this->PF_DiMuon_Muon2_Index.clear();
  this->PF_DiMuon_Vertex->Clear();
  this->PF_DiMuon_VertexProb.clear();
  this->PF_DiMuon_DCA.clear();
  this->PF_DiMuon_MassError.clear();
  // PF MET
  this->PF_MET_P2 = TVector2();
  // PF Muon-MET
  this->PF_MuonMET_P4T->Clear();
  // Gen Particle
  this->Gen_Particle_P4->Clear();
  this->Gen_Particle_PdgId.clear();
  this->Gen_Particle_Status.clear();
  this->Gen_Particle_Mother_Index.clear();
  this->Gen_Particle_Daughter_Index.clear();
  // Gen Muon
  this->Gen_Muon_N = 0;
  this->Gen_Muon_P4->Clear();
  this->Gen_Muon_Charge.clear();
  this->Gen_Muon_Particle_Index.clear();
  this->Gen_Muon_Reco_Index.clear();
  this->Gen_Muon_PF_Index.clear();
}


DEFINE_FWK_MODULE(HiMuonAnalyzer);

