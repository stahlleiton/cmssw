// -*- C++ -*-
//
// Package:    HiConversionAnalyzer
// Class:      HiConversionAnalyzer
// 
/**\class HiConversionAnalyzer HIConversionAnalyzer.cc

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

#include "HeavyIonsAnalysis/ConversionAnalysis/interface/HiConversionAnalyzer.h"

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
HiConversionAnalyzer::HiConversionAnalyzer(const edm::ParameterSet& iConfig):
  _recoMuonsToken        ( consumes< edm::View< reco::Muon         > >( iConfig.getParameter< edm::InputTag >( "recoMuonsTag"       ) )),
  _recoConversionsToken  ( consumes< edm::View< reco::Conversion   > >( iConfig.getParameter< edm::InputTag >( "recoConversionsTag" ) )),
  _genParticlesToken     ( consumes< reco::GenParticleCollection     >( iConfig.getParameter< edm::InputTag >( "genParticlesTag"    ) )),
  _pfCandidatesToken     ( consumes< edm::View< reco::PFCandidate  > >( iConfig.getParameter< edm::InputTag >( "pfCandidatesTag"    ) )),
  _primaryVertexToken    ( consumes< edm::View< reco::Vertex       > >( iConfig.getParameter< edm::InputTag >( "primaryVertexTag"   ) )),
  _beamSpotToken         ( consumes< reco::BeamSpot                  >( iConfig.getParameter< edm::InputTag >( "beamSpotTag"        ) )),
  _doAll                 ( iConfig.getParameter< bool >( "doAll" ) )
{
}

HiConversionAnalyzer::~HiConversionAnalyzer()
{
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
HiConversionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  StringBoolMap doConversion;

  // Fill Event Tree always
  doConversion["Event"] = true;

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

  // If Gen Particle collection exist, fill the Gen Conversion tree
  edm::Handle< reco::GenParticleCollection > genParticleHandle;
  getCollection(iEvent, _genParticlesToken, genParticleHandle);
  reco::GenParticleRefVector genParticleRefCollection;
  if (genParticleHandle.isValid()) {
    for (auto it = genParticleHandle->begin(); it != genParticleHandle->end(); ++it) {
      const reco::GenParticle& genParticle = *it;
      if ( ( (abs(genParticle.pdgId()) >  10) && (abs(genParticle.pdgId()) < 19) ) || // Keep all leptons
           ( (abs(genParticle.pdgId()) >  21) && (abs(genParticle.pdgId()) < 25) ) || // Keep all EWQ bosons
           ( (abs(genParticle.pdgId()) == 10441) || (abs(genParticle.pdgId()) == 20443) || (abs(genParticle.pdgId()) == 445) ) || // Keep Chi_C(1P)
           ( (abs(genParticle.pdgId()) == 443) || (abs(genParticle.pdgId()) == 100443) ) || // Keep Charmonia
           ( (abs(genParticle.pdgId()) == 111) ) // Keep pi0
           ){
        genParticleRefCollection.push_back( reco::GenParticleRef(genParticleHandle, it - genParticleHandle->begin()) );
      }
    }
    doConversion["Gen"] = true;
  }
  else { doConversion["Gen"] = false; }
  reco::GenParticleCollection genPhotonCollection;
  for (uint i = 0; i < genParticleRefCollection.size(); i++ ) {
    if ( (abs(genParticleRefCollection.at(i)->pdgId()) == 22) && (genParticleRefCollection.at(i)->status()==1) ) genPhotonCollection.push_back( *(genParticleRefCollection.at(i)) );
  }

  // If Reco Muon collection exist, use them
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
  }
  std::sort(recoMuonCollection.begin(), recoMuonCollection.end(), recoMuonPTComparator_);

  // If Reco Conversion collection exist, fill the Reco Conversion tree
  edm::Handle< edm::View< reco::Conversion > > recoConversionHandle;
  getCollection(iEvent, _recoConversionsToken, recoConversionHandle);
  reco::ConversionCollection recoConversionCollection;
  if (recoConversionHandle.isValid()) {
    for (uint i = 0; i < recoConversionHandle->size(); i++ ) {
      const reco::Conversion& recoConversion = recoConversionHandle->at(i);
      if ( (recoConversion.tracks().size() == 2) && recoConversion.conversionVertex().isValid() ) {
        recoConversionCollection.push_back( recoConversion );
      }
    }
    doConversion["Reco"] = true;
  }
  else { doConversion["Reco"] = false; }
  std::sort(recoConversionCollection.begin(), recoConversionCollection.end(), recoConversionGreaterByChi2_);

  // If PF Particle collection exist, fill the PF Conversion tree
  edm::Handle<  edm::View<reco::PFCandidate> > pfCandidateHandle;
  getCollection(iEvent, _pfCandidatesToken, pfCandidateHandle);
  reco::PFCandidateCollection pfCandidateCollection;
  if (pfCandidateHandle.isValid()) {
    for (uint i = 0; i < pfCandidateHandle->size(); i++ ) {
      const reco::PFCandidate& pfCandidate = pfCandidateHandle->at(i);
      pfCandidateCollection.push_back( pfCandidate );
    }
    doConversion["PF"] = true; 
  }
  else { doConversion["PF"] = false; }
  std::sort(pfCandidateCollection.begin(), pfCandidateCollection.end(), pfCandidatePTComparator_);
  reco::PFCandidateCollection pfPhotonCollection;
  for (uint i = 0; i < pfCandidateCollection.size(); i++ ) {
    if (pfCandidateCollection.at(i).particleId() == reco::PFCandidate::gamma) pfPhotonCollection.push_back( pfCandidateCollection.at(i) );
  }
  doConversion["PF"] = false;

  // Check if at least one of the collections was found
  if (!doConversion["Reco"] && !doConversion["PF"] && !doConversion["Gen"]) { 
    edm::LogError("HiConversionAnalyzer_ProductNotFound") << "No collections were found!" << std::endl;
    return;
  }

  // Do matching between conversion collections
  // Do matching between conversion collections
  IndexMap indexMap_GEN, indexMap_PF, indexMap_RECO;
  std::vector< std::vector< Char_t > > matchRECOGEN , matchRECOPF , matchPFGEN;
  matchRECOGEN = doMatching(recoConversionCollection, genPhotonCollection, 0.5,  0.5  );
  matchRECOPF  = doMatching(recoConversionCollection, pfPhotonCollection,  0.5,  0.5  );
  matchPFGEN   = doMatching(pfPhotonCollection,       genPhotonCollection, 0.5,  0.5  );
  indexMap_RECO["GEN"] = matchRECOGEN[0];
  indexMap_GEN["RECO"] = matchRECOGEN[1];
  indexMap_RECO["PF"]  = matchRECOPF[0];
  indexMap_PF["RECO"]  = matchRECOPF[1];
  indexMap_PF["GEN"]   = matchPFGEN[0];
  indexMap_GEN["PF"] = matchPFGEN[1];

  std::vector< std::string > treeNames;
  if (_doAll) { treeNames.push_back("All"); }
  else {
    if (doConversion["Event"]) treeNames.push_back("Event");
    if (doConversion["Reco"])  treeNames.push_back("Reco");
    if (doConversion["PF"])    treeNames.push_back("PF");
    if (doConversion["Gen"])   treeNames.push_back("Gen");
  }

  for (uint i = 0; i < treeNames.size(); i++) {
    std::string name = treeNames[i];
    if (firstEvent_) convEvt_[name].IniArrays();
    convEvt_[name].Clear();
    if (doConversion["Event"] && (name=="Event" || name=="All")) {
      convEvt_[name].SetPrimaryVertex    ( priVtx );
      convEvt_[name].SetVertexCollection ( primaryVertexCollection );
      convEvt_[name].Fill( iEvent, iSetup ); 
    }
    if (doConversion["Reco"]  && (name=="Reco" || name=="All")) {
      convEvt_[name].IniTrackBuilder     ( iSetup );
      convEvt_[name].SetPrimaryVertex    ( priVtx );
      convEvt_[name].SetVertexCollection ( primaryVertexCollection );
      convEvt_[name].Fill ( recoConversionCollection, indexMap_RECO, pfPhotonCollection, recoMuonCollection );
    }
    if (doConversion["PF"]   && (name=="PF" || name=="All")) {
      convEvt_[name].IniTrackBuilder     ( iSetup );
      convEvt_[name].SetPrimaryVertex    ( priVtx );
      convEvt_[name].SetVertexCollection ( primaryVertexCollection );
      convEvt_[name].Fill( pfCandidateCollection, indexMap_PF ); 
    }
    if (doConversion["Gen"]  && (name=="Gen" || name=="All")) {
      convEvt_[name].Fill( genParticleRefCollection, indexMap_GEN );
    }
    if (firstEvent_) {
      convTree_[name] = fs_->make<TTree>(Form("Conversion_%s", name.c_str()), "");
      convEvt_[name].SetTree( convTree_[name] );
      convEvt_[name].SetBranches( name , doConversion);
    }
  }
  if (firstEvent_) { firstEvent_ = false; }

  std::map< std::string, TTree* >::iterator iter;
  for ( iter = convTree_.begin(); iter != convTree_.end(); iter++ ) {
    (iter->second)->Fill();
  }
}

//--------------------------------------------------------------------------------------------------
void 
HiConversionAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

//--------------------------------------------------------------------------------------------------
void 
HiConversionAnalyzer::beginJob()
{
  firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ---------------------
void
HiConversionAnalyzer::endJob() 
{
}

//
// constructors and destructor
//
HiConversionEvent::HiConversionEvent()
{
}

HiConversionEvent::~HiConversionEvent()
{
}

//--------------------------------------------------------------------------------------------------
void
HiConversionEvent::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
HiConversionEvent::IniTrackBuilder(const edm::EventSetup& iSetup)
{
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", _theTTBuilder);
}

//--------------------------------------------------------------------------------------------------
bool 
HiConversionEvent::CheckTkVtxCompatible(const reco::Conversion& conv)
{
  std::vector< std::pair< double, short> > idx[2];
  short ik=-1;
  for (auto& tk : conv.tracks()) {
    ik++;
    short count=-1;
    for (auto& vtx : vtxCol_) {
      count++;
      double dz = tk->dz(vtx.position());
      double dzError = tk->dzError();
      dzError = sqrt( dzError*dzError + vtx.covariance(2,2) );
      //if(fabs(dz_)/dzError_ > sigmaTkVtxComp_) continue;
      idx[ik].push_back(std::pair<double,short>(fabs(dz) , count));
    }
    if(idx[ik].size()==0) { return false; }
    std::stable_sort(idx[ik].begin() , idx[ik].end() , lowerByFirstElement);
  }
  if ( (idx[0][0].second == idx[1][0].second) || (idx[0][1].second == idx[1][0].second) || (idx[0][0].second == idx[1][1].second) ) return true;
  return false;
}

//--------------------------------------------------------------------------------------------------
bool 
HiConversionEvent::FoundCompInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB)
{
  size_t count=0;
  uint32_t oldSubStr=0;
  for (int i = 0; ( (i < hitPatA.numberOfHits(reco::HitPattern::HitCategory::TRACK_HITS)) && (count < 2) ); i++ ) {
    uint32_t hitA = hitPatA.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS , i);
    if ( !hitPatA.validHitFilter(hitA) || !hitPatA.trackerHitFilter(hitA) ) continue;
    if ( (hitPatA.getSubStructure(hitA) == oldSubStr) && (hitPatA.getLayer(hitA) == oldSubStr) ) continue;
    if ( hitPatB.getTrackerMonoStereo(reco::HitPattern::HitCategory::TRACK_HITS,hitPatA.getSubStructure(hitA),hitPatA.getLayer(hitA)) != 0 ) return true;
    oldSubStr = hitPatA.getSubStructure(hitA);
    count++;
  }
  return false;
}

//--------------------------------------------------------------------------------------------------
bool 
HiConversionEvent::ConversionEqualByTrack(const reco::Conversion& c1, const reco::Conversion& c2)
{
  for (auto& tk1 : c1.tracks() ) { for (auto& tk2 : c2.tracks() ) { if(tk1 == tk2) { return true; } } }
  return false;
}

//--------------------------------------------------------------------------------------------------
void
HiConversionEvent::Fill(const reco::ConversionCollection& recoConversionCollection, const IndexMap& indexMap, const reco::PFCandidateCollection& pfPhotonCollection, const reco::MuonCollection& recoMuonCollection)
{
  this->Reco_Conversion_N = recoConversionCollection.size();
  for (ushort iconv = 0; iconv < recoConversionCollection.size(); iconv++) {
    const reco::Conversion& recoConversion = recoConversionCollection.at(iconv);
    // Reco Conversion Information
    TLorentzVector convP4 = TLorentzVector();
    convP4.SetPtEtaPhiM( recoConversion.refittedPair4Momentum().pt() , recoConversion.refittedPair4Momentum().eta() , recoConversion.refittedPair4Momentum().phi() , 0.0 );
    TVector3 convVtx( recoConversion.conversionVertex().position().x() , recoConversion.conversionVertex().position().y() , recoConversion.conversionVertex().position().z() );
    int charge1 = 0; for (auto const& track: recoConversion.tracks()) { charge1 += track->charge(); }
    new ((*this->Reco_Conversion_P4)[iconv]) TLorentzVector( convP4 );
    new ((*this->Reco_Conversion_Vertex)[iconv]) TVector3( convVtx );
    this->Reco_Conversion_Charge.push_back    ( charge1 );
    this->Reco_Conversion_Gen_Index.push_back ( indexMap.at("GEN")[iconv] );
    this->Reco_Conversion_PF_Index.push_back  ( indexMap.at("PF")[iconv]  );
    // Reco Conversion Flags
    this->Reco_Conversion_isDuplicate.push_back ( false );
    this->Reco_Conversion_isPionLeg.push_back   ( false );
    for (auto& pfPhoton : pfPhotonCollection) {
      TLorentzVector phoP4 = TLorentzVector(); phoP4.SetPtEtaPhiM( pfPhoton.pt() , pfPhoton.eta() , pfPhoton.phi() , 0.0 );
      if (fabs((phoP4+convP4).M() - 0.1349766) < 0.05) { this->Reco_Conversion_isPionLeg.at(iconv) = true; break; }
    }    
    this->Reco_Conversion_isPriVtxCompatible.push_back ( CheckTkVtxCompatible(recoConversion) );
    bool flagCompatibleInnerHits = false;
    if (recoConversion.tracks().size()==2) {
      reco::HitPattern hitPatA = recoConversion.tracks().at(0)->hitPattern();
      reco::HitPattern hitPatB = recoConversion.tracks().at(1)->hitPattern();
      if ( FoundCompInnerHits(hitPatA,hitPatB) && FoundCompInnerHits(hitPatB,hitPatA) ) flagCompatibleInnerHits = true;
    }
    this->Reco_Conversion_hasCompatibleInnerHits.push_back ( flagCompatibleInnerHits );
    UInt_t quality = 0;
    for (uint i=0; i<12; i++) { if (recoConversion.quality(reco::Conversion::ConversionQuality(i))) { quality += (1 << i); } }
    this->Reco_Conversion_Quality.push_back ( quality );
    // Reco Conversion Quality Variables
    this->Reco_Conversion_NSharedHits.push_back      ( recoConversion.nSharedHits()            );
    this->Reco_Conversion_Algo.push_back             ( UChar_t(recoConversion.algo())          );
    this->Reco_Conversion_MVA.push_back              ( recoConversion.MVAout()                 );
    this->Reco_Conversion_DCA.push_back              ( recoConversion.distOfMinimumApproach()  );
    this->Reco_Conversion_VertexProb.push_back       ( TMath::Prob( recoConversion.conversionVertex().chi2(), recoConversion.conversionVertex().ndof() ) );
    this->Reco_Conversion_dXY.push_back              ( recoConversion.dxy (privtx_.position()) );
    this->Reco_Conversion_dZ.push_back               ( recoConversion.dz  (privtx_.position()) );
    this->Reco_Conversion_lXY.push_back              ( recoConversion.lxy (privtx_.position()) );
    this->Reco_Conversion_lZ.push_back               ( recoConversion.lz  (privtx_.position()) );
    this->Reco_Conversion_dPhiTracksAtVtx.push_back  ( recoConversion.dPhiTracksAtVtx()        );
    this->Reco_Conversion_pairCotThetaSep.push_back  ( recoConversion.pairCotThetaSeparation() );
    this->Reco_Conversion_EoverP.push_back           ( recoConversion.EoverPrefittedTracks()   );
    // Reco Conversion Tracks Information
    this->Reco_Conversion_NTracks.push_back ( recoConversion.nTracks() );
    if (recoConversion.tracks()[0].isNonnull() && recoConversion.tracks()[0].isAvailable()) {
      TVector3 p3 = TVector3();
      p3.SetPtEtaPhi( recoConversion.tracks()[0]->pt() , recoConversion.tracks()[0]->eta() , recoConversion.tracks()[0]->phi() );
      new ((*this->Reco_Conversion_Track1_P3)[iconv]) TVector3( p3 );
      this->Reco_Conversion_Track1_PtError.push_back          ( recoConversion.tracks()[0]->ptError()                                   );
      this->Reco_Conversion_Track1_NHitsBeforeVtx.push_back   ( recoConversion.nHitsBeforeVtx().at(0)                                   );
      this->Reco_Conversion_Track1_isHighPurity.push_back     ( recoConversion.tracks()[0]->quality(reco::TrackBase::highPurity)        );
      this->Reco_Conversion_Track1_NumOfValHits.push_back     ( recoConversion.tracks()[0]->numberOfValidHits()                         );
      this->Reco_Conversion_Track1_NumOfLostHits.push_back    ( recoConversion.tracks()[0]->numberOfLostHits()                          );
      this->Reco_Conversion_Track1_NumOfValPixHits.push_back  ( recoConversion.tracks()[0]->hitPattern().numberOfValidPixelHits()       );
      this->Reco_Conversion_Track1_TrkLayersWithMea.push_back ( recoConversion.tracks()[0]->hitPattern().trackerLayersWithMeasurement() );
      this->Reco_Conversion_Track1_PixLayersWithMea.push_back ( recoConversion.tracks()[0]->hitPattern().pixelLayersWithMeasurement()   );
      this->Reco_Conversion_Track1_dXY.push_back              ( recoConversion.tracks()[0]->dxy(privtx_.position())                     );
      this->Reco_Conversion_Track1_dXYErr.push_back           ( recoConversion.tracks()[0]->dxyError()                                  );
      this->Reco_Conversion_Track1_dZ.push_back               ( recoConversion.tracks()[0]->dz(privtx_.position())                      );
      this->Reco_Conversion_Track1_dZErr.push_back            ( recoConversion.tracks()[0]->dzError()                                   );
      this->Reco_Conversion_Track1_ValFrac.push_back          ( recoConversion.tracks()[0]->validFraction()                             );
      this->Reco_Conversion_Track1_NormChi2.push_back         ( recoConversion.tracks()[0]->normalizedChi2()                            );
    }
    else {
      new ((*this->Reco_Conversion_Track1_P3)[iconv]) TVector3();
      this->Reco_Conversion_Track1_PtError.push_back          ( -1. );
      this->Reco_Conversion_Track1_NHitsBeforeVtx.push_back   ( -1. );
      this->Reco_Conversion_Track1_isHighPurity.push_back     (  0  );
      this->Reco_Conversion_Track1_NumOfValHits.push_back     ( -1  );
      this->Reco_Conversion_Track1_NumOfLostHits.push_back    ( -1  );
      this->Reco_Conversion_Track1_NumOfValPixHits.push_back  ( -1  );
      this->Reco_Conversion_Track1_TrkLayersWithMea.push_back ( -1  );
      this->Reco_Conversion_Track1_PixLayersWithMea.push_back ( -1  );
      this->Reco_Conversion_Track1_dXY.push_back              ( -99.);
      this->Reco_Conversion_Track1_dXYErr.push_back           ( -1. );
      this->Reco_Conversion_Track1_dZ.push_back               ( -99.);
      this->Reco_Conversion_Track1_dZErr.push_back            ( -1. );
      this->Reco_Conversion_Track1_ValFrac.push_back          ( -1. );
      this->Reco_Conversion_Track1_NormChi2.push_back         ( -1. );
    }
    if (recoConversion.tracks()[1].isNonnull() && recoConversion.tracks()[1].isAvailable()) {
      TVector3 p3 = TVector3();
      p3.SetPtEtaPhi( recoConversion.tracks()[1]->pt() , recoConversion.tracks()[1]->eta() , recoConversion.tracks()[1]->phi() );
      new ((*this->Reco_Conversion_Track2_P3)[iconv]) TVector3( p3 );
      this->Reco_Conversion_Track2_PtError.push_back          ( recoConversion.tracks()[1]->ptError()                                   );
      this->Reco_Conversion_Track2_NHitsBeforeVtx.push_back   ( recoConversion.nHitsBeforeVtx().at(1)                                   );
      this->Reco_Conversion_Track2_isHighPurity.push_back     ( recoConversion.tracks()[1]->quality(reco::TrackBase::highPurity)        );
      this->Reco_Conversion_Track2_NumOfValHits.push_back     ( recoConversion.tracks()[1]->numberOfValidHits()                         );
      this->Reco_Conversion_Track2_NumOfLostHits.push_back    ( recoConversion.tracks()[1]->numberOfLostHits()                          );
      this->Reco_Conversion_Track2_NumOfValPixHits.push_back  ( recoConversion.tracks()[1]->hitPattern().numberOfValidPixelHits()       );
      this->Reco_Conversion_Track2_TrkLayersWithMea.push_back ( recoConversion.tracks()[1]->hitPattern().trackerLayersWithMeasurement() );
      this->Reco_Conversion_Track2_PixLayersWithMea.push_back ( recoConversion.tracks()[1]->hitPattern().pixelLayersWithMeasurement()   );
      this->Reco_Conversion_Track2_dXY.push_back              ( recoConversion.tracks()[1]->dxy(privtx_.position())                     );
      this->Reco_Conversion_Track2_dXYErr.push_back           ( recoConversion.tracks()[1]->dxyError()                                  );
      this->Reco_Conversion_Track2_dZ.push_back               ( recoConversion.tracks()[1]->dz(privtx_.position())                      );
      this->Reco_Conversion_Track2_dZErr.push_back            ( recoConversion.tracks()[1]->dzError()                                   );
      this->Reco_Conversion_Track2_ValFrac.push_back          ( recoConversion.tracks()[1]->validFraction()                             );
      this->Reco_Conversion_Track2_NormChi2.push_back         ( recoConversion.tracks()[1]->normalizedChi2()                            );
    }
    else {
      new ((*this->Reco_Conversion_Track2_P3)[iconv]) TVector3();
      this->Reco_Conversion_Track2_PtError.push_back          ( -1. );
      this->Reco_Conversion_Track2_NHitsBeforeVtx.push_back   ( -1. );
      this->Reco_Conversion_Track2_isHighPurity.push_back     (  0  );
      this->Reco_Conversion_Track2_NumOfValHits.push_back     ( -1  );
      this->Reco_Conversion_Track2_NumOfLostHits.push_back    ( -1  );
      this->Reco_Conversion_Track2_NumOfValPixHits.push_back  ( -1  );
      this->Reco_Conversion_Track2_TrkLayersWithMea.push_back ( -1  );
      this->Reco_Conversion_Track2_PixLayersWithMea.push_back ( -1  );
      this->Reco_Conversion_Track2_dXY.push_back              ( -99.);
      this->Reco_Conversion_Track2_dXYErr.push_back           ( -1. );
      this->Reco_Conversion_Track2_dZ.push_back               ( -99.);
      this->Reco_Conversion_Track2_dZErr.push_back            ( -1. );
      this->Reco_Conversion_Track2_ValFrac.push_back          ( -1. );
      this->Reco_Conversion_Track2_NormChi2.push_back         ( -1. );
    }
    // Reco DiConversion
    for (ushort iconv2 = iconv+1; iconv2 < recoConversionCollection.size(); iconv2++) {
      const reco::Conversion& recoConversion2 = recoConversionCollection.at(iconv2);
      // Reco DiConversion Information
      TLorentzVector conv2P4 = TLorentzVector();
      conv2P4.SetPtEtaPhiM( recoConversion2.refittedPair4Momentum().pt() , recoConversion2.refittedPair4Momentum().eta() , recoConversion2.refittedPair4Momentum().phi() , 0.0 );
      int charge2 = 0; for (auto const& track: recoConversion2.tracks()) { charge2 += track->charge(); }
      TLorentzVector diConversionP4 = TLorentzVector();
      diConversionP4 = convP4 + conv2P4;
      if (this->Reco_Conversion_isDuplicate.at(iconv) == false) { this->Reco_Conversion_isDuplicate.at(iconv) = ConversionEqualByTrack( recoConversion, recoConversion ); }
      Char_t charge = charge1 + charge2;
      TVector3 diConversionVtx = TVector3();
      Float_t vProb    = -1.;
      Float_t massWErr = -1.;
      if ( ( recoConversion.tracks()[0].isNonnull()  && recoConversion.tracks()[0].isAvailable()   ) && 
           ( recoConversion.tracks()[1].isNonnull()  && recoConversion.tracks()[1].isAvailable()   ) && 
           ( recoConversion2.tracks()[0].isNonnull()  && recoConversion2.tracks()[0].isAvailable() ) && 
           ( recoConversion2.tracks()[1].isNonnull()  && recoConversion2.tracks()[1].isAvailable() ) ) {
        std::vector<reco::TransientTrack> t_tks;
        t_tks.push_back( _theTTBuilder->build( *recoConversion.tracks()[0])  );
        t_tks.push_back( _theTTBuilder->build( *recoConversion.tracks()[1])  );
        t_tks.push_back( _theTTBuilder->build( *recoConversion2.tracks()[0]) );
        t_tks.push_back( _theTTBuilder->build( *recoConversion2.tracks()[1]) );
        KalmanVertexFitter  vtxFitter = KalmanVertexFitter(true);
        TransientVertex diConversionVertex = vtxFitter.vertex( t_tks );
        if (diConversionVertex.isValid()) {
          Float_t vChi2 = diConversionVertex.totalChiSquared();
          Float_t vNDF  = diConversionVertex.degreesOfFreedom();
          vProb = TMath::Prob(vChi2,(int)vNDF);
          diConversionVtx.SetXYZ( diConversionVertex.position().x() , diConversionVertex.position().y() , diConversionVertex.position().z() );
        }
        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );
        InvariantMassFromVertex massCalculator;
        massWErr = massCalculator.invariantMass( VtxForInvMass, _eleMasses ).error();
      }
      ushort idiconv = this->Reco_DiConversion_Charge.size();
      this->Reco_DiConversion_N = idiconv + 1;
      new ((*this->Reco_DiConversion_P4)[idiconv]) TLorentzVector( diConversionP4 );
      this->Reco_DiConversion_Charge.push_back      ( charge    );
      new ((*this->Reco_DiConversion_Vertex)[idiconv]) TVector3( diConversionVtx );
      this->Reco_DiConversion_VertexProb.push_back  ( vProb     );
      this->Reco_DiConversion_MassError.push_back   ( massWErr  );
      this->Reco_DiConversion_Conversion1_Index.push_back ( iconv  );
      this->Reco_DiConversion_Conversion2_Index.push_back ( iconv2 );
    }
    if (
        (this->Reco_Conversion_isDuplicate.at(iconv)==false) &&
        (this->Reco_Conversion_hasCompatibleInnerHits.at(iconv)==true)
        ) {
      for (ushort imuon = 0, idimuon = -1; imuon < recoMuonCollection.size(); imuon++) {
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
          idimuon += 1;
          if ( 
              ( charge == 0 && vProb > 0.01 ) && // Opposite sign dimuons and Vertex probability larger than 1%
              ( (fabs(muP4.Eta())<2.4 && muP4.Pt()>1.0) && (fabs(mu2P4.Eta())<2.4 && mu2P4.Pt()>1.0) ) && // Muons within CMS acceptance
              ( muon::isSoftMuon(recoMuon , privtx_) && muon::isSoftMuon(recoMuon2 , privtx_) ) &&           // Muons passing ID
              ( fabs(diMuonP4.M()-3.096916)<0.3 || fabs(diMuonP4.M()-9.46030)<0.3 ) // Jpsi + Upsilon
               ) {
            if ( fabs(diMuonP4.M()-9.46030)<0.3 && !( (fabs(muP4.Eta())<2.4 && muP4.Pt()>4.0) && (fabs(mu2P4.Eta())<2.4 && mu2P4.Pt()>4.0) ) ) continue;
            if ( fabs(diMuonP4.M()-3.096916)<0.3 && !( (fabs(muP4.Eta())<2.4 && muP4.Pt()>1.3) && (fabs(mu2P4.Eta())<2.4 && mu2P4.Pt()>1.3) ) ) continue;
            int idimuonconv = this->Reco_DiMuonConv_DiMuon_Index.size();
            TLorentzVector diMuonConvP4 = TLorentzVector(); diMuonConvP4 = diMuonP4 + convP4;
            this->Reco_DiMuonConv_Conversion_Index.push_back(iconv);
            this->Reco_DiMuonConv_DiMuon_Index.push_back(idimuon);
            new ((*this->Reco_DiMuonConv_P4)[idimuonconv]) TLorentzVector( diMuonConvP4 );
            this->Reco_Chi_Mass.push_back(-1.);
            this->Reco_Chi_Type.push_back(0);
            if (fabs(diMuonP4.M()-3.096916)<0.3) { 
              this->Reco_Chi_Mass.at(idimuonconv) = ( diMuonConvP4.M() - diMuonP4.M() + 3.096916);
              this->Reco_Chi_Type.at(idimuonconv) = 1;
            }
            if (fabs(diMuonP4.M()-9.46030)<0.3 ) {
              this->Reco_Chi_Mass.at(idimuonconv) = ( diMuonConvP4.M() - diMuonP4.M() + 9.46030 );
              this->Reco_Chi_Type.at(idimuonconv) = 2;
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

