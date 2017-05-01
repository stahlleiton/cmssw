// -*- C++ -*-
//
// Package:    HiBcAnalyzer
// Class:      HiBcAnalyzer
// 
/**\class HiBcAnalyzer HiBcAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Andre, Stahl 
//         Created:  Apr  08 2017
// 
//
//

#include "HeavyIonsAnalysis/ConversionAnalysis/interface/HiBcAnalyzer.h"

//Headers for the core items
#include "FWCore/Framework/interface/MakerMacros.h"

//Headers for the data items
#include "DataFormats/TrackReco/interface/TrackBase.h"

//Headers for services and tools
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

//system include files
#include <iostream>


//
// constructors and destructor
//
HiBcAnalyzer::HiBcAnalyzer(const edm::ParameterSet& iConfig):
  _recoMuonsToken        ( consumes< edm::View< reco::Muon         > >( iConfig.getParameter< edm::InputTag >( "recoMuonsTag"       ) )),
  _recoTracksToken       ( consumes< edm::View< reco::Track        > >( iConfig.getParameter< edm::InputTag >( "recoTracksTag"      ) )),
  _pfMETToken            ( consumes< edm::View< reco::PFMET        > >( iConfig.getParameter< edm::InputTag >( "pfMETTag"           ) )),
  _primaryVertexToken    ( consumes< edm::View< reco::Vertex       > >( iConfig.getParameter< edm::InputTag >( "primaryVertexTag"   ) )),
  _beamSpotToken         ( consumes< reco::BeamSpot                  >( iConfig.getParameter< edm::InputTag >( "beamSpotTag"        ) )),
  _doAll                 ( iConfig.getParameter< bool >( "doAll" ) )
{
}

HiBcAnalyzer::~HiBcAnalyzer()
{
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
HiBcAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  StringBoolMap doBc;
  
  // Fill Event Tree always
  doBc["Event"] = true;
  
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
      if ( (iEvent.isRealData() == false) || 
           ( (primaryVertex.isFake() == false) && (primaryVertex.tracksSize() >= 1) && (fabs(primaryVertex.z()) <= 25) && (primaryVertex.position().rho() <= 2) )
           ) {
        primaryVertexCollection.push_back( primaryVertex );
      }
    }
  }
  const reco::Vertex& priVtx = ( ( primaryVertexCollection.size() > 0 ) ? primaryVertexCollection.at(0) : beamSpot );

  // If Reco Muon collection exist, use them
  edm::Handle< edm::View< reco::Muon > > recoMuonHandle;
  getCollection(iEvent, _recoMuonsToken, recoMuonHandle);
  reco::MuonCollection recoMuonCollection;
  if (recoMuonHandle.isValid()) {
    for (uint i = 0; i < recoMuonHandle->size(); i++ ) {
      const reco::Muon& recoMuon = recoMuonHandle->at(i);
      if (
          ( fabs(recoMuon.Eta())<2.4 && recoMuon.Pt() > 0.8 ) &&
          ( recoMuon.isGlobalMuon() || recoMuon.isTrackerMuon() || recoMuon.isPFMuon() )
          ) {
        recoMuonCollection.push_back( recoMuon );
      }
    }
    doBc["BcToJpsiMuNu"] = true; doBc["BcToJpsiPi"] = true; doBc["BcToJpsiTriPi"] = true;
  }
  else { doBc["BcToJpsiMuNu"] = false; doBc["BcToJpsiPi"] = false; doBc["BcToJpsiTriPi"] = false; }

  // If PF MET tree exist, fill it
  edm::Handle<  edm::View<reco::PFMET> > pfMETHandle;
  getCollection(iEvent, _pfMETToken, pfMETHandle);
  if (pfMETHandle.isValid() && pfMETHandle->size()!=1) {
    edm::LogError("METAnalyzer_PFMETCollection_corrupt") << "Number of PF METs in event is different from 1!" << std::endl;
    return;
  }
  if (pfMETHandle.isValid() && pfMETHandle->size()==1) { doBc["BcToJpsiMuNu"] = true; }
  else { doBc["BcToJpsiMuNu"] = false; }

  // If Reco Track collection exist, use them
  edm::Handle< edm::View< reco::Track > > recoTrackHandle;
  getCollection(iEvent, _recoTracksToken, recoTrackHandle);
  reco::TrackCollection recoTrackCollection;
  if (recoTrackHandle.isValid()) {
    for (uint i = 0; i < recoTrackHandle->size(); i++ ) {
      const reco::Track& recoTrack = recoTrackHandle->at(i);
      if (
          ( fabs(recoTrack.Eta())<2.4 && recoTrack.Pt() > 0.5 ) &&
          ( recoTrack.quality(reco::TrackBase::highPurity)  ) 
          ){
        recoTrackCollection.push_back( recoTrack );
      }
    }
    doBc["BcToJpsiPi"] = true; doBc["BcToJpsiTriPi"] = true;
  }
  else { doBc["BcToJpsiPi"] = false; doBc["BcToJpsiTriPi"] = false; }

  // Check if at least one of the collections was found
  if (!doBc.at("BcToJpsiMuNu") && !doBc.at("BcToJpsiPi") && !doBc.at("BcToJpsiTriPi")) { 
    edm::LogError("HiBcAnalyzer_ProductNotFound") << "No collections were found!" << std::endl;
    return;
  }

  std::vector< std::string > treeNames;
  if (_doAll) { treeNames.push_back("All"); }
  else {
    if (doBc["Event"])         treeNames.push_back("Event");
    if (doBc["BcToJpsiMuNu"])  treeNames.push_back("JpsiMuNu");
    if (doBc["BcToJpsiPi"])    treeNames.push_back("JpsiPi");
    if (doBc["BcToJpsiTriPi"]) treeNames.push_back("JpsiTriPi");
  }

  for (uint i = 0; i < treeNames.size(); i++) {
    std::string name = treeNames[i];
    if (firstEvent_) BcEvt_[name].IniArrays();
    BcEvt_[name].Clear();
    if (doBc.at("Event") && (name=="Event" || name=="All")) {
      BcEvt_[name].SetPrimaryVertex    ( priVtx );
      BcEvt_[name].SetVertexCollection ( primaryVertexCollection );
      BcEvt_[name].Fill( iEvent, iSetup ); 
    }
    if (doBc.at("BcToJpsiMuNu")  && (name=="JpsiMuNu" || name=="All")) {
      BcEvt_[name].IniTrackBuilder     ( iSetup );
      BcEvt_[name].SetPrimaryVertex    ( priVtx );
      BcEvt_[name].SetVertexCollection ( primaryVertexCollection );
      BcEvt_[name].Fill ( recoMuonCollection, pfMET );
    }
    if (doBc.at("BcToJpsiPi")  && (name=="JpsiPi" || name=="All")) {
      BcEvt_[name].IniTrackBuilder     ( iSetup );
      BcEvt_[name].SetPrimaryVertex    ( priVtx );
      BcEvt_[name].SetVertexCollection ( primaryVertexCollection );
      BcEvt_[name].Fill ( recoMuonCollection, recoTrackCollection, "JpsiPi" );
    }
    if (doBc.at("BcToJpsiTriPi")  && (name=="JpsiTriPi" || name=="All")) {
      BcEvt_[name].IniTrackBuilder     ( iSetup );
      BcEvt_[name].SetPrimaryVertex    ( priVtx );
      BcEvt_[name].SetVertexCollection ( primaryVertexCollection );
      BcEvt_[name].Fill ( recoMuonCollection, recoTrackCollection, "JpsiTriPi" );
    }
    if (firstEvent_) {
      BcTree_[name] = fs_->make<TTree>(Form("Bc_%s", name.c_str()), "");
      BcEvt_[name].SetTree( BcTree_[name] );
      BcEvt_[name].SetBranches( name , doBc);
    }
  }
  if (firstEvent_) { firstEvent_ = false; }

  std::map< std::string, TTree* >::iterator iter;
  for ( iter = BcTree_.begin(); iter != BcTree_.end(); iter++ ) {
    (iter->second)->Fill();
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiBcAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

//--------------------------------------------------------------------------------------------------
void 
HiBcAnalyzer::beginJob()
{
  firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ---------------------
void
HiBcAnalyzer::endJob() 
{
}

//
// constructors and destructor
//
HiBcEvent::HiConversionEvent()
{
}

HiBcEvent::~HiConversionEvent()
{
}

//--------------------------------------------------------------------------------------------------
void
HiBcEvent::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  this->Event_nRun   = iEvent.id().run();
  this->Event_nLumi  = iEvent.luminosityBlock();
  this->Event_nBX    = std::max(iEvent.bunchCrossing(),0);
  this->Event_nOrbit = std::max(iEvent.orbitNumber(),0);
  this->Event_nEvent = iEvent.id().event();
  // Primary Vertex Information
  this->Event_PriVtx_Position = TVector3( privtx_.x()      , privtx_.y()      , privtx_.z()      );
  this->Event_PriVtx_Error    = TVector3( privtx_.xError() , privtx_.yError() , privtx_.zError() );
  this->Event_nPV = vtxCol_.size();
}
 
//--------------------------------------------------------------------------------------------------
void 
HiBcEvent::IniTrackBuilder(const edm::EventSetup& iSetup)
{
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", _theTTBuilder);
}
 
 //--------------------------------------------------------------------------------------------------
void
HiBcEvent::Fill(const reco::MuonCollection& recoMuonCollection, const reco:PFMET& pfMET)
{
  // Check MET
  TVector2 pfMETP2 = TVector2(pfMET.px(), pfMET.py());
  if (pfMETP2.Mod() < 20) return false;
  // Check Isolated Muon
  reco::MuonCollection isoMuonCollection;
  for (ushort imuonIso = 0; imuonIso < recoMuonCollection.size(); imuonIso++) {
    const reco::Muon& recoMuonIso = recoMuonCollection.at(imuonIso);
    TLorentzVector muIsoP4 = TLorentzVector();
    muIsoPo4.SetPtEtaPhiM( recoMuonIso.pt() , recoMuonIso.eta() , recoMuonIso.phi() , recoMuonIso.mass() );
    Float_t isoR03NoPUCorr = recoMuonIso.pfIsolationR03().sumChargedHadronPt + std::max(0.0f, recoMuonIso.pfIsolationR03().sumNeutralHadronEt + recoMuonIso.pfIsolationR03().sumPhotonEt );
    if (
        ( fabs(muIsoP4.Eta())<2.4 && muIsoP4.Pt()>0.8 ) && // Muon within CMS acceptance
        ( muon::isMediumMuon ( recoMuon ) ) &&             // Muon passing ID
        ( recoMuon.isPFIsolationValid() &&  ((isoR03NoPUCorr/recoMuonIso.pt()) < 0.15) ) // Muon is Isolated
        )
      {
        isMuonCollection.push_back(recoMuonIso);
        this->Reco_IsoMuon_Muon_Index.push_back(imuonIso);
        this->Reco_IsoMuon_Charge.push_back(imuonIso.charge());
        new ((*this->Reco_IsoMuon_P4)[imuonIso]) TLorentzVector( muIsoP4 );
      }
  }
  // Check Dimuons
  for (ushort imuon = 0, idimuon = 0; imuon < recoMuonCollection.size(); imuon++) {
    const reco::Muon& recoMuon = recoMuonCollection.at(imuon);
    TLorentzVector muP4 = TLorentzVector();
    muP4.SetPtEtaPhiM( recoMuon.pt() , recoMuon.eta() , recoMuon.phi() , recoMuon.mass() );
    for (ushort imuon2 = imuon+1; imuon2 < recoMuonCollection.size(); imuon2++) {
      const reco::Muon& recoMuon2 = recoMuonCollection.at(imuon2);
      TLorentzVector mu2P4 = TLorentzVector();
      mu2P4.SetPtEtaPhiM( recoMuon2.pt() , recoMuon2.eta() , recoMuon2.phi() , recoMuon2.mass() );
      TLorentzVector diMuonP4 = TLorentzVector();
      diMuonP4 = muP4 + mu2P4;
      Char_t charge = recoMuon.charge() + recoMuon2.charge();
      Float_t vProb    = -1. ;
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
        }
      }
      idimuon++;
      if ( 
          ( charge == 0 && vProb > 0.01 ) && // Opposite sign dimuons and Vertex probability larger than 1%
          ( (fabs(muP4.Eta())<2.4 && muP4.Pt()>0.8) && (fabs(mu2P4.Eta())<2.4 && mu2P4.Pt()>0.8) ) && // Muons within CMS acceptance
          ( muon::isSoftMuon(recoMuon , privtx_) && muon::isSoftMuon(recoMuon2 , privtx_) ) &&           // Muons passing ID
          ( fabs(diMuonP4.M()-3.096916)<0.2 || fabs(diMuonP4.M()-9.46030)<0.2 ) // Jpsi + Upsilon
           ) {
        new ((*this->Reco_DiMuon_P4)[imuonIso]) TLorentzVector( diMuonP4 );
        // Loop over the isolated muons
        for (ushort imuonIso = 0; imuonIso < recoMuonIsoCollection.size(); imuonIso++) {
          const reco::Muon& recoMuonIso = recoMuonIsoCollection.at(imuonIso);
          TLorentzVector muIsoP4 = TLorentzVector();
          muIsoP4.SetPtEtaPhiM( recoMuonIso.pt() , recoMuonIso.eta() , recoMuonIso.phi() , recoMuonIso.mass() );
          if (
              ( (recoMuonIso.charge() == recoMuon.charge() ) && (muIsoP4 != muP4)  ) ||
              ( (recoMuonIso.charge() == recoMuon2.charge()) && (muIsoP4 != mu2P4) ) 
              )
            {
              int itrimuon = this->Reco_TriMuon_Charge.size();
              TLorentzVector triMuonP4 = muIsoP4 + diMuonP4;
              new ((*this->Reco_TriMuon_P4)[itrimuon]) TLorentzVector( triMuonP4 );
              this->Reco_TriMuon_Charge.push_back(imuonIso.charge());
              this->Reco_TriMuon_DiMuon_Index.push_back(idimuon);
              this->Reco_TriMuon_IsoMuon_Index.push_back(imuonIso);
              this->Reco_Bc_Mass.push_back(-1.);
              this->Reco_Bc_Type.push_back(0);
              if (fabs(diMuonP4.M()-3.096916)<0.2) { 
                this->Reco_Bc_Mass.at(itrimuon) = ( triMuonP4.M() - diMuonP4.M() + 3.096916);
                this->Reco_BC_Type.at(itrimuon) = 1;
              }
              if (fabs(diMuonP4.M()-9.46030)<0.2 ) {
                this->Reco_Bc_Mass.at(itrimuon) = ( triMuonConvP4.M() - diMuonP4.M() + 9.46030 );
                this->Reco_Bc_Type.at(itrimuon) = 2;
              }
            }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
HiConversionEvent::Fill(const reco::PFCandidateCollection& pfCandidateCollection, const IndexMap& indexMap)
{
  // Fill PF Photon collection
  reco::PFCandidateCollection pfPhotonCollection;
  for (uint i = 0; i < pfCandidateCollection.size(); i++ ) {
    if ( pfCandidateCollection.at(i).particleId() == reco::PFCandidate::gamma ) pfPhotonCollection.push_back( pfCandidateCollection.at(i) );
  }
  // PF Photon
  this->PF_Photon_N = pfPhotonCollection.size();
  for (ushort ipho = 0, idipho = 0; ipho < pfPhotonCollection.size(); ipho++) {
    const reco::PFCandidate& pfPhoton = pfPhotonCollection.at(ipho);
    // PF Photon Information
    TLorentzVector phoP4 = TLorentzVector();
    phoP4.SetPtEtaPhiM( pfPhoton.pt() , pfPhoton.eta() , pfPhoton.phi() , pfPhoton.mass() );
    new ((*this->PF_Photon_P4)[ipho]) TLorentzVector( phoP4 );
    this->PF_Photon_Gen_Index.push_back  ( indexMap.at("GEN")[ipho]  );
    this->PF_Photon_Reco_Index.push_back ( indexMap.at("RECO")[ipho] );
    // PF DiPhoton
    for (ushort ipho2 = ipho+1; ipho2 < pfPhotonCollection.size(); ipho2++) {
      const reco::PFCandidate& pfPhoton2 = pfPhotonCollection.at(ipho2);
      // PF DiPhoton Information
      TLorentzVector pho2P4 = TLorentzVector();
      pho2P4.SetPtEtaPhiM( pfPhoton2.pt() , pfPhoton2.eta() , pfPhoton2.phi() , pfPhoton2.mass() );
      TLorentzVector diPhotonP4 = phoP4 + pho2P4;
      TVector3 diPhotonVtx = TVector3();
      this->PF_DiPhoton_N = idipho + 1;
      new ((*this->PF_DiPhoton_P4)[idipho]) TLorentzVector( diPhotonP4 );
      this->PF_DiPhoton_Photon1_Index.push_back ( ipho  );
      this->PF_DiPhoton_Photon2_Index.push_back ( ipho2 );
      idipho++;
    }
  }
}

//--------------------------------------------------------------------------------------------------
reco::GenParticleRef 
HiConversionEvent::findMotherRef(const reco::GenParticle& genParticle)
{
  if (genParticle.numberOfMothers() <= 0)                return reco::GenParticleRef();
  reco::GenParticleRef genMother = genParticle.motherRef();
  if (genMother.isNull() || !genMother.isAvailable())    return reco::GenParticleRef();
  if (genMother->pdgId() != genParticle.pdgId())         return genMother;
  return findMotherRef( *genMother );
}

//--------------------------------------------------------------------------------------------------
void
HiConversionEvent::Fill(const reco::GenParticleRefVector& genRefCollection, const IndexMap& indexMap)
{
  // Loop over the Gen Particle Collection
  std::vector< UChar_t > genPhotonIndex;
  reco::GenParticleCollection genPhotonCollection;
  for (uint ipar = 0; ipar < genRefCollection.size(); ipar++ ) {
    const reco::GenParticle& genParticle = *(genRefCollection.at(ipar));
    // Fill the Gen Photon Collection
    if ( (abs(genParticle.pdgId()) == 22) && (genParticle.status() == 1) ) { genPhotonCollection.push_back( genParticle ); genPhotonIndex.push_back( ipar ); }
  }
  // Gen Photon
  this->Gen_Photon_N = genPhotonCollection.size();
  for (ushort ipho = 0; ipho < genPhotonCollection.size(); ipho++ ) {
    const reco::GenParticle& genPhoton = genPhotonCollection.at(ipho);
    // Gen Photon Information
    TLorentzVector phoP4 = TLorentzVector();
    phoP4.SetPtEtaPhiM( genPhoton.pt() , genPhoton.eta() , genPhoton.phi() , genPhoton.mass() );
    new ((*this->Gen_Photon_P4)[ipho]) TLorentzVector( phoP4 );
    this->Gen_Photon_Particle_Index.push_back ( genPhotonIndex[ipho]  );
    this->Gen_Photon_Reco_Index.push_back ( indexMap.at("RECO")[ipho] );
  }
  // Gen DiPhoton
  uint idipho = 0;
  for (uint  ipho = 0; ipho < genPhotonCollection.size(); ipho++ ) {
    const reco::GenParticle&    genPhoton = genPhotonCollection.at(ipho);
    const reco::GenParticleRef& genMom    = findMotherRef(genPhoton);
    for (uint  ipho2 = ipho+1; ipho2 < genPhotonCollection.size(); ipho2++ ) {
      const reco::GenParticle&    genPhoton2 = genPhotonCollection.at(ipho2);
      const reco::GenParticleRef& genMom2    = findMotherRef(genPhoton2);
      // Gen DiPhoton Information
      if ( genMom.isNonnull() && genMom.isAvailable() && (genMom == genMom2) ) {
        TLorentzVector diPhotonP4 = TLorentzVector();
        diPhotonP4.SetPtEtaPhiM( genMom->pt() , genMom->eta() , genMom->phi() , genMom->mass() );
        UInt_t pdgId = abs(genMom->pdgId());
        this->Gen_DiPhoton_N = idipho + 1;
        new ((*this->Gen_DiPhoton_P4)[idipho]) TLorentzVector( diPhotonP4 );
        this->Gen_DiPhoton_PdgId.push_back         ( pdgId    );
        this->Gen_DiPhoton_Photon1_Index.push_back ( ipho    );
        this->Gen_DiPhoton_Photon2_Index.push_back ( ipho2   );
        idipho++;
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiConversionEvent::IniArrays()
{
  // Initialize the TCloneArrays
  this->Reco_Conversion_P4             = new TClonesArray("TLorentzVector", 1000);
  this->Reco_Conversion_Vertex         = new TClonesArray("TVector3",       1000);
  this->Reco_Conversion_Track1_P3      = new TClonesArray("TVector3",       1000);
  this->Reco_Conversion_Track2_P3      = new TClonesArray("TVector3",       1000);
  this->Reco_DiConversion_P4           = new TClonesArray("TLorentzVector", 10000);
  this->Reco_DiConversion_Vertex       = new TClonesArray("TVector3",       10000);
  this->Reco_DiMuonConv_P4             = new TClonesArray("TLorentzVector", 10000);
  this->PF_Photon_P4                   = new TClonesArray("TLorentzVector", 1000);
  this->PF_DiPhoton_P4                 = new TClonesArray("TLorentzVector", 10000);
  this->Gen_Photon_P4                  = new TClonesArray("TLorentzVector", 100);
  this->Gen_DiPhoton_P4                = new TClonesArray("TLorentzVector", 100);
}

//--------------------------------------------------------------------------------------------------
void 
HiConversionEvent::SetBranches(const std::string name, const StringBoolMap& doConversion)
{
  // Set Branches
  if ( doConversion.at("Event") && (name == "Event" || name == "All") ) {
    // Event Info
    tree_->Branch("Event_Run",                             &(this->Event_nRun),               "Event_Run/i"           );
    tree_->Branch("Event_Lumi",                            &(this->Event_nLumi),              "Event_Lumi/s"          );
    tree_->Branch("Event_Bx",                              &(this->Event_nBX),                "Event_Bx/i"            );
    tree_->Branch("Event_Orbit",                           &(this->Event_nOrbit),             "Event_Orbit/l"         );
    tree_->Branch("Event_Number",                          &(this->Event_nEvent),             "Event_Number/l"        );
    tree_->Branch("Event_nPV",                             &(this->Event_nPV),                "Event_nPV/b"           );
    tree_->Branch("Event_PriVtx_Pos",                      "TVector3",             &(this->Event_PriVtx_Position)     );
    tree_->Branch("Event_PriVtx_Err",                      "TVector3",             &(this->Event_PriVtx_Error   )     );
  }
  if ( doConversion.at("Reco") && (name == "Reco" || name == "All") ) {
    // Reco Conversion Kinematic
    tree_->Branch("Reco_Conversion_N",                     &(this->Reco_Conversion_N),        "Reco_Conversion_N/s"   );
    tree_->Branch("Reco_Conversion_Mom",                  "TClonesArray",   &(this->Reco_Conversion_P4),     32000, 0 );
    tree_->Branch("Reco_Conversion_Vertex",               "TClonesArray",   &(this->Reco_Conversion_Vertex), 32000, 0 );
    tree_->Branch("Reco_Conversion_Charge",                &(this->Reco_Conversion_Charge)                            );
    // Reco Conversion Matched Index
    if (doConversion.at("Gen")) tree_->Branch("Reco_Conversion_Gen_Idx",    &(this->Reco_Conversion_Gen_Index)        );
    if (doConversion.at("PF") ) tree_->Branch("Reco_Conversion_PF_Idx",     &(this->Reco_Conversion_PF_Index)         );
    // Reco Conversion ID Flags
    tree_->Branch("Reco_Conversion_isDuplicate",           &(this->Reco_Conversion_isDuplicate)                       );
    tree_->Branch("Reco_Conversion_isPionLeg",             &(this->Reco_Conversion_isPionLeg)                         );
    tree_->Branch("Reco_Conversion_isPriVtxComp",          &(this->Reco_Conversion_isPriVtxCompatible)                );
    tree_->Branch("Reco_Conversion_hasCompInnerHits",      &(this->Reco_Conversion_hasCompatibleInnerHits)            );
    tree_->Branch("Reco_Conversion_Quality",               &(this->Reco_Conversion_Quality)                           );
    // Reco Conversion ID vars                                                                                           
    tree_->Branch("Reco_Conversion_NSharedHits",           &(this->Reco_Conversion_NSharedHits)                       );
    tree_->Branch("Reco_Conversion_Algo",                  &(this->Reco_Conversion_Algo)                              );
    tree_->Branch("Reco_Conversion_MVA",                   &(this->Reco_Conversion_MVA)                               );
    tree_->Branch("Reco_Conversion_DCA",                   &(this->Reco_Conversion_DCA)                               );
    tree_->Branch("Reco_Conversion_VertexProb",            &(this->Reco_Conversion_VertexProb)                        );
    tree_->Branch("Reco_Conversion_dZ",                    &(this->Reco_Conversion_dXY)                               );
    tree_->Branch("Reco_Conversion_dXY",                   &(this->Reco_Conversion_dZ)                                );
    tree_->Branch("Reco_Conversion_lXY",                   &(this->Reco_Conversion_lXY)                               );
    tree_->Branch("Reco_Conversion_lZ",                    &(this->Reco_Conversion_lZ)                                );
    tree_->Branch("Reco_Conversion_dPhiTracksAtVtx",       &(this->Reco_Conversion_dPhiTracksAtVtx)                   );
    tree_->Branch("Reco_Conversion_pairCotThetaSep",       &(this->Reco_Conversion_pairCotThetaSep)                   );
    tree_->Branch("Reco_Conversion_EoverP",                &(this->Reco_Conversion_EoverP)                            );
    tree_->Branch("Reco_Conversion_NTracks",               &(this->Reco_Conversion_NTracks)                           );
    // Reco Conversion Track 1
    tree_->Branch("Reco_Conversion_Track1_Mom",       "TClonesArray", &(this->Reco_Conversion_Track1_P3), 32000, 0    );
    tree_->Branch("Reco_Conversion_Track1_PtErr",          &(this->Reco_Conversion_Track1_PtError)                    );
    tree_->Branch("Reco_Conversion_Track1_HitsBeforeVtx",  &(this->Reco_Conversion_Track1_NHitsBeforeVtx)             );
    tree_->Branch("Reco_Conversion_Track1_isHighPurity",   &(this->Reco_Conversion_Track1_isHighPurity)               );
    tree_->Branch("Reco_Conversion_Track1_ValidHits",      &(this->Reco_Conversion_Track1_NumOfValHits)               );
    tree_->Branch("Reco_Conversion_Track1_LostHits",       &(this->Reco_Conversion_Track1_NumOfLostHits)              );
    tree_->Branch("Reco_Conversion_Track1_ValidPixHits",   &(this->Reco_Conversion_Track1_NumOfValPixHits)            );
    tree_->Branch("Reco_Conversion_Track1_TrkLayers",      &(this->Reco_Conversion_Track1_TrkLayersWithMea)           );
    tree_->Branch("Reco_Conversion_Track1_PixLayers",      &(this->Reco_Conversion_Track1_PixLayersWithMea)           );
    tree_->Branch("Reco_Conversion_Track1_dXY",            &(this->Reco_Conversion_Track1_dXY)                        );
    tree_->Branch("Reco_Conversion_Track1_dXYErr",         &(this->Reco_Conversion_Track1_dXYErr)                     );
    tree_->Branch("Reco_Conversion_Track1_dZ",             &(this->Reco_Conversion_Track1_dZ)                         );
    tree_->Branch("Reco_Conversion_Track1_dZErr",          &(this->Reco_Conversion_Track1_dZErr)                      );
    tree_->Branch("Reco_Conversion_Track1_ValFrac",        &(this->Reco_Conversion_Track1_ValFrac)                    );
    tree_->Branch("Reco_Conversion_Track1_NormChi2",       &(this->Reco_Conversion_Track1_NormChi2)                   );
    // Reco Conversion Track 2
    tree_->Branch("Reco_Conversion_Track2_Mom",       "TClonesArray", &(this->Reco_Conversion_Track2_P3), 32000, 0    );
    tree_->Branch("Reco_Conversion_Track2_PtErr",          &(this->Reco_Conversion_Track2_PtError)                    );
    tree_->Branch("Reco_Conversion_Track2_HitsBeforeVtx",  &(this->Reco_Conversion_Track2_NHitsBeforeVtx)             );
    tree_->Branch("Reco_Conversion_Track2_isHighPurity",   &(this->Reco_Conversion_Track2_isHighPurity)               );
    tree_->Branch("Reco_Conversion_Track2_ValidHits",      &(this->Reco_Conversion_Track2_NumOfValHits)               );
    tree_->Branch("Reco_Conversion_Track2_LostHits",       &(this->Reco_Conversion_Track2_NumOfLostHits)              );
    tree_->Branch("Reco_Conversion_Track2_ValidPixHits",   &(this->Reco_Conversion_Track2_NumOfValPixHits)            );
    tree_->Branch("Reco_Conversion_Track2_TrkLayers",      &(this->Reco_Conversion_Track2_TrkLayersWithMea)           );
    tree_->Branch("Reco_Conversion_Track2_PixLayers",      &(this->Reco_Conversion_Track2_PixLayersWithMea)           );
    tree_->Branch("Reco_Conversion_Track2_dXY",            &(this->Reco_Conversion_Track2_dXY)                        );
    tree_->Branch("Reco_Conversion_Track2_dXYErr",         &(this->Reco_Conversion_Track2_dXYErr)                     );
    tree_->Branch("Reco_Conversion_Track2_dZ",             &(this->Reco_Conversion_Track2_dZ)                         );
    tree_->Branch("Reco_Conversion_Track2_dZErr",          &(this->Reco_Conversion_Track2_dZErr)                      );
    tree_->Branch("Reco_Conversion_Track2_ValFrac",        &(this->Reco_Conversion_Track2_ValFrac)                    );
    tree_->Branch("Reco_Conversion_Track2_NormChi2",       &(this->Reco_Conversion_Track2_NormChi2)                   );
    // Reco DiConversion
    tree_->Branch("Reco_DiConversion_N",                   &(this->Reco_DiConversion_N), "Reco_DiConversion_N/i"      );
    tree_->Branch("Reco_DiConversion_Mom",                 "TClonesArray", &(this->Reco_DiConversion_P4), 32000, 0    );
    tree_->Branch("Reco_DiConversion_Charge",              &(this->Reco_DiConversion_Charge)                          );
    tree_->Branch("Reco_DiConversion_Conversion1_Idx",     &(this->Reco_DiConversion_Conversion1_Index)               );
    tree_->Branch("Reco_DiConversion_Conversion2_Idx",     &(this->Reco_DiConversion_Conversion2_Index)               );
    tree_->Branch("Reco_DiConversion_Vertex",             "TClonesArray", &(this->Reco_DiConversion_Vertex), 32000, 0 );
    tree_->Branch("Reco_DiConversion_VtxProb",             &(this->Reco_DiConversion_VertexProb)                      );
    tree_->Branch("Reco_DiConversion_MassErr",             &(this->Reco_DiConversion_MassError)                       );
    // JUST FOR FUN, DiMuon+Conversion
    tree_->Branch("Reco_DiMuonConv_Mom",                   "TClonesArray", &(this->Reco_DiMuonConv_P4), 32000, 0      );
    tree_->Branch("Reco_DiMuonConv_Conversion_Idx",        &(this->Reco_DiMuonConv_Conversion_Index)                  );
    tree_->Branch("Reco_DiMuonConv_DiMuon_Idx",            &(this->Reco_DiMuonConv_DiMuon_Index)                      );
    tree_->Branch("Reco_Chi_Mass",                         &(this->Reco_Chi_Mass)                                     );
    tree_->Branch("Reco_Chi_Type",                         &(this->Reco_Chi_Type)                                     );
  }
  if ( doConversion.at("PF") && (name == "PF" || name == "All") ) {
    // PF Photon
    tree_->Branch("PF_Photon_N",                           &(this->PF_Photon_N),           "PF_Photon_N/s"            );
    tree_->Branch("PF_Photon_Mom",                         "TClonesArray", &(this->PF_Photon_P4), 32000, 0            );
    if (doConversion.at("Gen") ) tree_->Branch("PF_Photon_Gen_Idx",         &(this->PF_Photon_Gen_Index)              );
    if (doConversion.at("Reco")) tree_->Branch("PF_Photon_Reco_Idx",        &(this->PF_Photon_Reco_Index)             );
    // PF DiPhoton
    tree_->Branch("PF_DiPhoton_N",                         &(this->PF_DiPhoton_N),         "PF_DiPhoton_N/i"          );
    tree_->Branch("PF_DiPhoton_Mom",                       "TClonesArray", &(this->PF_DiPhoton_P4), 32000, 0          );
    tree_->Branch("PF_DiPhoton_Photon1_Idx",               &(this->PF_DiPhoton_Photon1_Index)                         );
    tree_->Branch("PF_DiPhoton_Photon2_Idx",               &(this->PF_DiPhoton_Photon2_Index)                         );
  }
  if ( doConversion.at("Gen") && (name == "Gen" || name == "All") ) {
    // Gen Photon
    tree_->Branch("Gen_Photon_N",                          &(this->Gen_Photon_N),          "Gen_Photon_N/b"           );
    tree_->Branch("Gen_Photon_Mom",                        "TClonesArray", &(this->Gen_Photon_P4), 32000, 0           );
    tree_->Branch("Gen_Photon_Particle_Idx",               &(this->Gen_Photon_Particle_Index)                         );
    if (doConversion.at("Reco")) tree_->Branch("Gen_Photon_Reco_Idx",   &(this->Gen_Photon_Reco_Index)                );
    // Gen DiPhoton
    tree_->Branch("Gen_DiPhoton_N",                        &(this->Gen_DiPhoton_N),          "Gen_DiPhoton_N/b"       );
    tree_->Branch("Gen_DiPhoton_Mom",                      "TClonesArray", &(this->Gen_DiPhoton_P4), 32000, 0         );
    tree_->Branch("Gen_DiPhoton_PdgId",                    &(this->Gen_DiPhoton_PdgId)                                );
    tree_->Branch("Gen_DiPhoton_Photon1_Idx",              &(this->Gen_DiPhoton_Photon1_Index)                        );
    tree_->Branch("Gen_DiPhoton_Photon2_Idx",              &(this->Gen_DiPhoton_Photon2_Index)                        );
  }
}

//--------------------------------------------------------------------------------------------------
void
HiConversionEvent::Clear(void)
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
  // Reco Conversion Kinematic
  this->Reco_Conversion_N = 0;
  this->Reco_Conversion_P4->Clear();
  this->Reco_Conversion_Vertex->Clear();
  this->Reco_Conversion_Charge.clear();
  // Reco Conversion Matched Index
  this->Reco_Conversion_Gen_Index.clear();
  this->Reco_Conversion_PF_Index.clear();
  // Reco Conversion ID Flags
  this->Reco_Conversion_isDuplicate.clear();
  this->Reco_Conversion_isPionLeg.clear();
  this->Reco_Conversion_isPriVtxCompatible.clear();
  this->Reco_Conversion_hasCompatibleInnerHits.clear();
  this->Reco_Conversion_Quality.clear();
  // Reco Conversion ID Vars
  this->Reco_Conversion_NSharedHits.clear();
  this->Reco_Conversion_Algo.clear();
  this->Reco_Conversion_MVA.clear();
  this->Reco_Conversion_DCA.clear();
  this->Reco_Conversion_VertexProb.clear();
  this->Reco_Conversion_dXY.clear();
  this->Reco_Conversion_dZ.clear();
  this->Reco_Conversion_lXY.clear();
  this->Reco_Conversion_lZ.clear();
  this->Reco_Conversion_dPhiTracksAtVtx.clear();
  this->Reco_Conversion_pairCotThetaSep.clear();
  this->Reco_Conversion_EoverP.clear();
  this->Reco_Conversion_NTracks.clear();
  // Reco Conversion Track 1
  this->Reco_Conversion_Track1_P3->Clear();
  this->Reco_Conversion_Track1_PtError.clear();
  this->Reco_Conversion_Track1_NHitsBeforeVtx.clear();
  this->Reco_Conversion_Track1_isHighPurity.clear();
  this->Reco_Conversion_Track1_NumOfValHits.clear();
  this->Reco_Conversion_Track1_NumOfLostHits.clear();
  this->Reco_Conversion_Track1_NumOfValPixHits.clear();
  this->Reco_Conversion_Track1_TrkLayersWithMea.clear();
  this->Reco_Conversion_Track1_PixLayersWithMea.clear();
  this->Reco_Conversion_Track1_dXY.clear();
  this->Reco_Conversion_Track1_dXYErr.clear();
  this->Reco_Conversion_Track1_dZ.clear();
  this->Reco_Conversion_Track1_dZErr.clear();
  this->Reco_Conversion_Track1_ValFrac.clear();
  this->Reco_Conversion_Track1_NormChi2.clear();
  // Reco Conversion Track 2
  this->Reco_Conversion_Track2_P3->Clear();
  this->Reco_Conversion_Track2_PtError.clear();
  this->Reco_Conversion_Track2_NHitsBeforeVtx.clear();
  this->Reco_Conversion_Track2_isHighPurity.clear();
  this->Reco_Conversion_Track2_NumOfValHits.clear();
  this->Reco_Conversion_Track2_NumOfLostHits.clear();
  this->Reco_Conversion_Track2_NumOfValPixHits.clear();
  this->Reco_Conversion_Track2_TrkLayersWithMea.clear();
  this->Reco_Conversion_Track2_PixLayersWithMea.clear();
  this->Reco_Conversion_Track2_dXY.clear();
  this->Reco_Conversion_Track2_dXYErr.clear();
  this->Reco_Conversion_Track2_dZ.clear();
  this->Reco_Conversion_Track2_dZErr.clear();
  this->Reco_Conversion_Track2_ValFrac.clear();
  this->Reco_Conversion_Track2_NormChi2.clear();
  // Reco DiConversion
  this->Reco_DiConversion_N = 0;
  this->Reco_DiConversion_P4->Clear();
  this->Reco_DiConversion_Charge.clear();
  this->Reco_DiConversion_Conversion1_Index.clear();
  this->Reco_DiConversion_Conversion2_Index.clear();
  this->Reco_DiConversion_Vertex->Clear();
  this->Reco_DiConversion_VertexProb.clear();
  this->Reco_DiConversion_MassError.clear();
  // JUST FOR FUN, THE CHI
  this->Reco_DiMuonConv_P4->Clear();
  this->Reco_DiMuonConv_Conversion_Index.clear();
  this->Reco_DiMuonConv_DiMuon_Index.clear();
  this->Reco_Chi_Mass.clear();
  this->Reco_Chi_Type.clear();
  // PF Photon
  this->PF_Photon_N = 0;
  this->PF_Photon_P4->Clear();
  this->PF_Photon_Gen_Index.clear();
  this->PF_Photon_Reco_Index.clear();
  // PF DiPhoton
  this->PF_DiPhoton_N = 0;
  this->PF_DiPhoton_P4->Clear();
  this->PF_DiPhoton_Photon1_Index.clear();
  this->PF_DiPhoton_Photon2_Index.clear();
  // Gen Photon
  this->Gen_Photon_N = 0;
  this->Gen_Photon_P4->Clear();
  this->Gen_Photon_Particle_Index.clear();
  this->Gen_Photon_Reco_Index.clear();
  this->Gen_Photon_PF_Index.clear();
  // Gen DiPhoton
  this->Gen_DiPhoton_N = 0;
  this->Gen_DiPhoton_P4->Clear();
  this->Gen_DiPhoton_PdgId.clear();
  this->Gen_DiPhoton_Photon1_Index.clear();
  this->Gen_DiPhoton_Photon2_Index.clear();
}


DEFINE_FWK_MODULE(HiConversionAnalyzer);

