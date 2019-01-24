// -*- C++ -*-
//
// Package:    TimeAnalyzer
// Class:      TimeAnalyzer
// 
/**\class TimeAnalyzer TimeAnalyzer.cc
Description: [one line class summary]
Implementation:
[Notes on implementation]
 */
//
// Original Author:  Andre, Stahl 
//         Created:  Jan  04 2019
// 
//
//

//Headers for the core items
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//Headers for the data items
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

//Headers for services and tools
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//system include files
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <map>


//
// class declaration
//

class TimeEvent
{
 public:
  // ===== Class Methods =====
  TimeEvent    ();
  ~TimeEvent   ();

  void   Clear ();

  void   SetTree      ( TTree* p ) { tree_ = p; }
  void   SetBranches  ( void );

  void   Fill         ( const reco::GenParticleRefVector& );
  void   Fill         ( const edm::View<reco::Track>&, const reco::VertexCollection&, 
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>& );
  void   Fill         ( const edm::Event&,             const reco::VertexCollection& );

 private:

  const double MASS_PION = 0.13957018;
  const double c_cm_ns = 2.99792458e1; //[cm/ns]

  TTree*                  tree_;

  // Event Info
  UInt_t                  Event_nRun;
  UShort_t                Event_nLumi;
  UInt_t                  Event_nBX;
  ULong64_t               Event_nOrbit;
  ULong64_t               Event_nEvent;
  // Primary Vertex Information
  TVector3                Event_PriVtx_Position;
  TVector3                Event_PriVtx_Error;
  UChar_t                 Event_nPV;
  // Track Info
  UInt_t                  Reco_Track_N;
  std::vector< Float_t >  Reco_Track_beta_PI;
  std::vector< Float_t >  Reco_Track_t0_PI;
  std::vector< Float_t >  Reco_Track_t0Err_PI;
  std::vector< Float_t >  Reco_Track_dt_PI;
  std::vector< Float_t >  Reco_Track_beta_PV;
  std::vector< Float_t >  Reco_Track_t0_PV;
  std::vector< Float_t >  Reco_Track_t0Err_PV;
  std::vector< Float_t >  Reco_Track_mass_PV;
  std::vector< Float_t >  Reco_Track_dt_PV;
  std::vector< Float_t >  Reco_Track_pathLength;
  std::vector< Float_t >  Reco_Track_p;
  std::vector< Float_t >  Reco_Track_tMTD;
  std::vector< Float_t >  Reco_Track_tMTDErr;
  std::vector< Float_t >  Reco_Track_eta;
  std::vector< Float_t >  Reco_Track_charge;
  std::vector< Int_t   >  Reco_Track_quality;
};

class TimeAnalyzer : public edm::EDAnalyzer
{
 public:
  explicit TimeAnalyzer(const edm::ParameterSet&);
  virtual ~TimeAnalyzer();
    
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  // TFileService
  edm::Service<TFileService> fs_;

  // Input Info
  const edm::EDGetTokenT< edm::View< reco::Track>> _recoTracksToken;
  const edm::EDGetTokenT< reco::VertexCollection > _primaryVertexToken;
  const edm::EDGetTokenT< reco::BeamSpot         > _beamSpotToken;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackBetaToken;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackT0Token;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackSigmaT0Token;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackTMTDToken;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackSigmaTMTDToken;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackMomToken;
  const edm::EDGetTokenT< edm::ValueMap< float > > _trackPathLengthToken;

  // Time Info Containers
  TTree*     timeTree_;
  TimeEvent  timeEvt_;

  // Flags
  bool firstEvent_;
};

//
// constructors and destructor
//
TimeAnalyzer::TimeAnalyzer(const edm::ParameterSet& iConfig):
  _recoTracksToken      ( consumes< edm::View< reco::Track >>( iConfig.getParameter< edm::InputTag >( "recoTracksTag"      ) )),
  _primaryVertexToken   ( consumes< reco::VertexCollection  >( iConfig.getParameter< edm::InputTag >( "primaryVertexTag"   ) )),
  _beamSpotToken        ( consumes< reco::BeamSpot          >( iConfig.getParameter< edm::InputTag >( "beamSpotTag"        ) )),
  _trackBetaToken       ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackBetaTag"       ) )),
  _trackT0Token         ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackT0Tag"         ) )),
  _trackSigmaT0Token    ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackSigmaT0Tag"    ) )),
  _trackTMTDToken       ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackTMTDTag"       ) )),
  _trackSigmaTMTDToken  ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackSigmaTMTDTag"  ) )),
  _trackMomToken        ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackMomTag"        ) )),
  _trackPathLengthToken ( consumes< edm::ValueMap< float  > >( iConfig.getParameter< edm::InputTag >( "trackPathLengthTag" ) ))
{
}

TimeAnalyzer::~TimeAnalyzer()
{
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
TimeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Initialize tree information
  timeEvt_.Clear();
  if (firstEvent_) {
    timeTree_ = fs_->make<TTree>("HiMTDTree", "");
    timeEvt_.SetTree( timeTree_ );
    timeEvt_.SetBranches();
  }

  // Extract the BeamSpot
  edm::Handle< reco::BeamSpot > beamSpotHandle;
  iEvent.getByToken(_beamSpotToken, beamSpotHandle);
  reco::Vertex beamSpot = reco::Vertex();
  if (beamSpotHandle.isValid()) beamSpot = reco::Vertex( beamSpotHandle->position(),  beamSpotHandle->covariance3D() );

  // Extract the Offline Primary Vertices
  edm::Handle< reco::VertexCollection > primaryVertexHandle;
  iEvent.getByToken(_primaryVertexToken, primaryVertexHandle);
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
  if (primaryVertexCollection.size() == 0) { primaryVertexCollection.push_back( beamSpot ); }
  timeEvt_.Fill(iEvent, primaryVertexCollection);

  // Extract the track information
  edm::Handle< edm::View< reco::Track > > recoTracksHandle;
  iEvent.getByToken(_recoTracksToken, recoTracksHandle);
  if (!recoTracksHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackBetaHandle;
  iEvent.getByToken(_trackBetaToken, trackBetaHandle);
  if (!trackBetaHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track beta valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackT0Handle;
  iEvent.getByToken(_trackT0Token, trackT0Handle);
  if (!trackT0Handle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track t0 valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackSigmaT0Handle;
  iEvent.getByToken(_trackSigmaT0Token, trackSigmaT0Handle);
  if (!trackSigmaT0Handle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track sigma t0 valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackTMTDHandle;
  iEvent.getByToken(_trackTMTDToken, trackTMTDHandle);
  if (!trackTMTDHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track tMTD valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackSigmaTMTDHandle;
  iEvent.getByToken(_trackSigmaTMTDToken, trackSigmaTMTDHandle);
  if (!trackSigmaTMTDHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track sigma tMTD valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackMomHandle;
  iEvent.getByToken(_trackMomToken, trackMomHandle);
  if (!trackMomHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track momentum valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > trackPathLengthHandle;
  iEvent.getByToken(_trackPathLengthToken, trackPathLengthHandle);
  if (!trackPathLengthHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track path length valueMap is not valid!" << std::endl; return; }
  //
  timeEvt_.Fill(*recoTracksHandle, primaryVertexCollection, *trackBetaHandle, *trackT0Handle, *trackSigmaT0Handle, *trackTMTDHandle, *trackSigmaTMTDHandle, *trackMomHandle, *trackPathLengthHandle);

  if (firstEvent_) { firstEvent_ = false; }
  timeTree_->Fill();
}

//--------------------------------------------------------------------------------------------------
void 
TimeAnalyzer::beginJob()
{
  firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ---------------------
void
TimeAnalyzer::endJob() 
{
}

//
// constructors and destructor
//
TimeEvent::TimeEvent()
{
}

TimeEvent::~TimeEvent()
{
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Fill(const edm::View<reco::Track>& recoTracks,   const reco::VertexCollection& vtxCol,
                const edm::ValueMap<float>& trackBeta,      const edm::ValueMap<float>& trackT0,
                const edm::ValueMap<float>& trackSigmaT0,   const edm::ValueMap<float>& trackTMTD,
                const edm::ValueMap<float>& trackSigmaTMTD, const edm::ValueMap<float>& trackMom,
                const edm::ValueMap<float>& trackPathLength
                )
{
  this->Reco_Track_N = 0;
  for (uint i = 0; i < recoTracks.size(); i++) {
    auto obj = recoTracks.ptrAt(i);
    const double beta_PION  = trackBeta[obj];
    const double t0_PION    = trackT0[obj];
    const double t0Err_PION = trackSigmaT0[obj];
    const double tMTD       = trackTMTD[obj];
    const double tMTDErr    = trackSigmaTMTD[obj];
    const double magp       = trackMom[obj];
    const double pathLength = trackPathLength[obj];
    const double t0_PV      = vtxCol.at(0).t();
    const double t0Err_PV   = vtxCol.at(0).tError();
    const double eta        = obj->eta();
    const double charge     = obj->charge();
    const int    quality    = obj->qualityMask();
    if (true) {//(pathLength>-1) {
      // Undo backpropagation
      const double dt_PION = (beta_PION!=0. ? (pathLength/beta_PION)*(1./c_cm_ns) : 1.0E+15);
      const double gammasq_PION = (std::abs(beta_PION)!=1. ? 1./(1.0-std::pow(beta_PION, 2.0)) : 1.0E+15);
      const double magp_T = std::sqrt(std::abs(gammasq_PION - 1.))*MASS_PION;
      const double thit_T = t0_PION + dt_PION;
      //if (pathLength>-1) {
        //if (std::abs(thit_T-tMTD)>1E-6) { std::cout << "[ERROR] thit: " << thit_T << " ,  tMTD: " << tMTD << std::endl; }
        //if (std::abs(magp_T-magp)>1E-6) { std::cout << "[ERROR] magp_T: " << magp_T << " ,  magp: " << magp << std::endl; }
      //}
      // Compute beta and mass using primary vertex info
      const double dt_PV = tMTD - t0_PV;
      const double beta_PV = (dt_PV!=0. ? (pathLength/dt_PV)*(1./c_cm_ns) : 1.0E+15);
      const double gammasq_PV = (std::abs(beta_PV)!=1. ? 1./(1.0-std::pow(beta_PV, 2.0)) : 1.0E+15);
      const double mass_PV = magp/std::sqrt(std::abs(gammasq_PV - 1.));
      // Store information
      this->Reco_Track_N++;
      this->Reco_Track_beta_PI.push_back(beta_PION);
      this->Reco_Track_t0_PI.push_back(t0_PION);
      this->Reco_Track_t0Err_PI.push_back(t0Err_PION);
      this->Reco_Track_dt_PI.push_back(dt_PION);
      this->Reco_Track_beta_PV.push_back(beta_PV);
      this->Reco_Track_t0_PV.push_back(t0_PV);
      this->Reco_Track_t0Err_PV.push_back(t0Err_PV);
      this->Reco_Track_mass_PV.push_back(mass_PV);
      this->Reco_Track_dt_PV.push_back(dt_PV);
      this->Reco_Track_pathLength.push_back(pathLength);
      this->Reco_Track_p.push_back(magp);
      this->Reco_Track_tMTD.push_back(tMTD);
      this->Reco_Track_eta.push_back(eta);
      this->Reco_Track_charge.push_back(charge);
      this->Reco_Track_quality.push_back(quality);
      // Print information
      //std::cout << "Track " << i << ">> pathLength: " << pathLength << " , thit: " << thit << " , t0Err: " << t0Err << " , beta_PION: " << beta_PION << " , t0_PION: " << t0_PION << std::endl;
      //std::cout << "Track " << i << ">> mass_PV: " << mass_PV << " , t0Err_PV: " << t0Err_PV << " , beta_PV: " << beta_PV << " , t0_PV: " << t0_PV << std::endl;
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Fill(const edm::Event& iEvent, const reco::VertexCollection& vtxCol)
{
  this->Event_nRun   = iEvent.id().run();
  this->Event_nLumi  = iEvent.luminosityBlock();
  this->Event_nBX    = std::max(iEvent.bunchCrossing(),0);
  this->Event_nOrbit = std::max(iEvent.orbitNumber(),0);
  this->Event_nEvent = iEvent.id().event();
  // Primary Vertex Information
  const reco::Vertex& privtx  = vtxCol.at(0);
  this->Event_PriVtx_Position = TVector3( privtx.x()      , privtx.y()      , privtx.z()      );
  this->Event_PriVtx_Error    = TVector3( privtx.xError() , privtx.yError() , privtx.zError() );
  this->Event_nPV = vtxCol.size();
}

//--------------------------------------------------------------------------------------------------
void 
TimeEvent::SetBranches(void)
{
  // Event Info
  tree_->Branch("Event_Run",             &(this->Event_nRun),   "Event_Run/i"    );
  tree_->Branch("Event_Lumi",            &(this->Event_nLumi),  "Event_Lumi/s"   );
  tree_->Branch("Event_Bx",              &(this->Event_nBX),    "Event_Bx/i"     );
  tree_->Branch("Event_Orbit",           &(this->Event_nOrbit), "Event_Orbit/l"  );
  tree_->Branch("Event_Number",          &(this->Event_nEvent), "Event_Number/l" );
  tree_->Branch("Event_nPV",             &(this->Event_nPV),    "Event_nPV/b"    );
  tree_->Branch("Event_PriVtx_Pos",   "TVector3", &(this->Event_PriVtx_Position) );
  tree_->Branch("Event_PriVtx_Err",   "TVector3", &(this->Event_PriVtx_Error )   );
  // Track Info
  tree_->Branch("Reco_Track_N",          &(this->Reco_Track_N), "Reco_Track_N/i" );
  tree_->Branch("Reco_Track_beta_PI",    &(this->Reco_Track_beta_PI)             );
  tree_->Branch("Reco_Track_t0_PI",      &(this->Reco_Track_t0_PI)               );
  tree_->Branch("Reco_Track_t0Err_PI",   &(this->Reco_Track_t0Err_PI)            );
  tree_->Branch("Reco_Track_dt_PI",      &(this->Reco_Track_dt_PI)               );
  tree_->Branch("Reco_Track_beta_PV",    &(this->Reco_Track_beta_PV)             );
  tree_->Branch("Reco_Track_t0_PV",      &(this->Reco_Track_t0_PV)               );
  tree_->Branch("Reco_Track_t0Err_PV",   &(this->Reco_Track_t0Err_PV)            );
  tree_->Branch("Reco_Track_mass_PV",    &(this->Reco_Track_mass_PV)             );
  tree_->Branch("Reco_Track_dt_PV",      &(this->Reco_Track_dt_PV)               );
  tree_->Branch("Reco_Track_pathLength", &(this->Reco_Track_pathLength)          );
  tree_->Branch("Reco_Track_p",          &(this->Reco_Track_p)                   );
  tree_->Branch("Reco_Track_tMTD",       &(this->Reco_Track_tMTD)                );
  tree_->Branch("Reco_Track_eta",        &(this->Reco_Track_eta)                 );
  tree_->Branch("Reco_Track_charge",     &(this->Reco_Track_charge)              );
  tree_->Branch("Reco_Track_quality",    &(this->Reco_Track_quality)             );
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Clear(void)
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
  // Track Info
  this->Reco_Track_N = 0;
  this->Reco_Track_beta_PI.clear();
  this->Reco_Track_t0_PI.clear();
  this->Reco_Track_t0Err_PI.clear();
  this->Reco_Track_dt_PI.clear();
  this->Reco_Track_beta_PV.clear();
  this->Reco_Track_t0_PV.clear();
  this->Reco_Track_t0Err_PV.clear();
  this->Reco_Track_mass_PV.clear();
  this->Reco_Track_dt_PV.clear();
  this->Reco_Track_pathLength.clear();
  this->Reco_Track_p.clear();
  this->Reco_Track_tMTD.clear();
  this->Reco_Track_eta.clear();
  this->Reco_Track_charge.clear();
  this->Reco_Track_quality.clear();
}

DEFINE_FWK_MODULE(TimeAnalyzer);
