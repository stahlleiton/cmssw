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
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

//Headers for the data items
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

//Headers for services and tools
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

//system include files
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
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
  void   incDiTracks  ( const bool& b ) { doDiTracks_ = b; }

  void   Fill         ( const edm::View<reco::Track>&, const reco::VertexCollection&, 
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<reco::DeDxData>&,
                        const edm::Handle<reco::GenParticleCollection>&,
                        const edm::Handle<float>&  ,   const TransientTrackBuilder& );
  void   Fill         ( const edm::Event&,             const reco::VertexCollection&, const float&, const int& );

 private:

  reco::GenParticleRef findMother(const reco::GenParticleRef&);

  const double MASS_PION = 0.13957018;
  const double MASS_KAON = 0.493677;
  const double MASS_PROT = 0.93827208;
  const double MASS_LC   = 2.28646;
  const double MASS_D0   = 1.86484;
  const double c_cm_ns   = 2.99792458e1; //[cm/ns]

  TTree*                  tree_;

  // Flags
  bool                    doDiTracks_ = false;

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
  // Centrality Information
  Float_t                 Event_hiHF;
  Int_t                   Event_hiBin;
  // Track Info
  UInt_t                  Reco_Track_N;
  std::vector< Float_t >  Reco_Track_beta_PI;
  std::vector< Float_t >  Reco_Track_t0_PI;
  std::vector< Float_t >  Reco_Track_t0Err_PI;
  std::vector< Float_t >  Reco_Track_beta_PV;
  std::vector< Float_t >  Reco_Track_t0_PV;
  std::vector< Float_t >  Reco_Track_t0Err_PV;
  std::vector< Float_t >  Reco_Track_massSq_PV;
  std::vector< Float_t >  Reco_Track_beta_PID;
  std::vector< Float_t >  Reco_Track_t0_PID;
  std::vector< Float_t >  Reco_Track_t0Err_PID;
  std::vector< Float_t >  Reco_Track_massSq_PID;
  std::vector< Float_t >  Reco_Track_pathLength;
  std::vector< Float_t >  Reco_Track_p;
  std::vector< Float_t >  Reco_Track_tMTD;
  std::vector< Float_t >  Reco_Track_tMTDErr;
  std::vector< Float_t >  Reco_Track_eta;
  std::vector< Float_t >  Reco_Track_phi;
  std::vector< Float_t >  Reco_Track_charge;
  std::vector< Int_t   >  Reco_Track_quality;
  std::vector< Float_t >  Reco_Track_pt;
  std::vector< Float_t >  Reco_Track_ptErr;
  std::vector< Float_t >  Reco_Track_dcaDz;
  std::vector< Float_t >  Reco_Track_dcaDxy;
  std::vector< Float_t >  Reco_Track_normChi2;
  std::vector< Int_t   >  Reco_Track_nHits;
  // Track dEdx
  std::vector< Float_t >  Reco_Track_dEdx;
  std::vector< Float_t >  Reco_Track_dEdxErr;
  std::vector< Int_t   >  Reco_Track_nSatMea;
  std::vector< Int_t   >  Reco_Track_nMea;
  // Gen Particle Info
  std::vector< Float_t >  Gen_Track_pdgId;
  std::vector< Float_t >  Gen_Track_pt;
  std::vector< Float_t >  Gen_Track_eta;
  std::vector< Float_t >  Gen_Track_phi;
  std::vector< Float_t >  Gen_Track_t0;
  std::vector< Float_t >  Gen_Track_charge;
  std::vector< Float_t >  Gen_Track_mass;
  std::vector< Float_t >  Gen_Track_momPdgId;
  // Pion+Kaon Info
  UInt_t                  Reco_DiTrack_N;
  std::vector< Int_t   >  Reco_DiTrack_Pion_idx;
  std::vector< Int_t   >  Reco_DiTrack_Kaon_idx;
  std::vector< Float_t >  Reco_DiTrack_pt;
  std::vector< Float_t >  Reco_DiTrack_rap;
  std::vector< Float_t >  Reco_DiTrack_mass;
  std::vector< Float_t >  Reco_DiTrack_vProb;
  std::vector< Float_t >  Reco_DiTrack_alpha;
  std::vector< Float_t >  Reco_DiTrack_d0Sig;
  // Gen Pion+Kaon Info
  std::vector< Float_t >  Gen_DiTrack_pdgId;
  std::vector< Float_t >  Gen_DiTrack_pt;
  std::vector< Float_t >  Gen_DiTrack_eta;
  std::vector< Float_t >  Gen_DiTrack_phi;
  std::vector< Float_t >  Gen_DiTrack_charge;
  std::vector< Float_t >  Gen_DiTrack_mass;
  std::vector< Int_t   >  Gen_DiTrack_isSwap;
  // Proton+Pion+Kaon Info
  UInt_t                  Reco_TriTrack_N;
  std::vector< Int_t   >  Reco_TriTrack_Pion_idx;
  std::vector< Int_t   >  Reco_TriTrack_Kaon_idx;
  std::vector< Int_t   >  Reco_TriTrack_Proton_idx;
  std::vector< Float_t >  Reco_TriTrack_pt;
  std::vector< Float_t >  Reco_TriTrack_rap;
  std::vector< Float_t >  Reco_TriTrack_mass;
  std::vector< Float_t >  Reco_TriTrack_vProb;
  std::vector< Float_t >  Reco_TriTrack_alpha;
  std::vector< Float_t >  Reco_TriTrack_d0Sig;
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
  const edm::EDGetTokenT< edm::ValueMap< float > > _tofPIDT0Token;
  const edm::EDGetTokenT< edm::ValueMap< float > > _tofPIDSigmaT0Token;
  const edm::EDGetTokenT< reco::GenParticleCollection> _genParticlesToken;
  const edm::EDGetTokenT< float                  > _genT0Token;
  const edm::EDGetTokenT<reco::Centrality        > _centralitySrcToken;
  const edm::EDGetTokenT< int                    > _centralityBinToken;
  const edm::EDGetTokenT< edm::ValueMap< reco::DeDxData > > _dEdxToken;

  // Flags
  const bool _doDiTracks;

  // Time Info Containers
  TTree*    timeTree_;
  TimeEvent timeEvt_;

  // Flags
  bool firstEvent_;
};

//
// constructors and destructor
//
TimeAnalyzer::TimeAnalyzer(const edm::ParameterSet& iConfig):
  _recoTracksToken      ( consumes< edm::View< reco::Track   > >( iConfig.getParameter< edm::InputTag >( "recoTracksTag"      ) )),
  _primaryVertexToken   ( consumes< reco::VertexCollection     >( iConfig.getParameter< edm::InputTag >( "primaryVertexTag"   ) )),
  _beamSpotToken        ( consumes< reco::BeamSpot             >( iConfig.getParameter< edm::InputTag >( "beamSpotTag"        ) )),
  _trackBetaToken       ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackBetaTag"       ) )),
  _trackT0Token         ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackT0Tag"         ) )),
  _trackSigmaT0Token    ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackSigmaT0Tag"    ) )),
  _trackTMTDToken       ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackTMTDTag"       ) )),
  _trackSigmaTMTDToken  ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackSigmaTMTDTag"  ) )),
  _trackMomToken        ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackMomTag"        ) )),
  _trackPathLengthToken ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "trackPathLengthTag" ) )),
  _tofPIDT0Token        ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "tofPIDT0Tag"        ) )),
  _tofPIDSigmaT0Token   ( consumes< edm::ValueMap< float     > >( iConfig.getParameter< edm::InputTag >( "tofPIDSigmaT0Tag"   ) )),
  _genParticlesToken    ( consumes< reco::GenParticleCollection>( iConfig.getParameter< edm::InputTag >( "genParticlesTag"    ) )),
  _genT0Token           ( consumes< float                      >( iConfig.getParameter< edm::InputTag >( "genT0Tag"           ) )),
  _centralitySrcToken   ( consumes< reco::Centrality           >( iConfig.getParameter< edm::InputTag >( "centralitySrcTag"   ) )),
  _centralityBinToken   ( consumes< int                        >( iConfig.getParameter< edm::InputTag >( "centralityBinTag"   ) )),
  _dEdxToken            ( consumes< edm::ValueMap< reco::DeDxData > >( iConfig.getParameter< edm::InputTag >( "dEdxTag"       ) )),
  _doDiTracks           ( iConfig.getParameter<bool>("doDiTracks") )
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
    timeEvt_.incDiTracks(_doDiTracks);
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
  if (primaryVertexHandle.isValid() && primaryVertexHandle->size()>0) {
    for (uint i = 0; i < primaryVertexHandle->size(); i++ ) {
      const reco::Vertex& primaryVertex = primaryVertexHandle->at(i);
      if ( (primaryVertex.isFake() == false) && (primaryVertex.tracksSize() >= 2) ) {
        primaryVertexCollection.push_back( primaryVertex );
      }
    }
  }
  if (primaryVertexCollection.size() == 0) { primaryVertexCollection.push_back( beamSpot ); }

  // Extract the Centrality Information
  edm::Handle<reco::Centrality> centralitySrcHandle;
  iEvent.getByToken(_centralitySrcToken, centralitySrcHandle);
  float hiHF = -1.;
  if (centralitySrcHandle.isValid()) { hiHF = centralitySrcHandle->EtHFtowerSum(); }

  // Extract the Centrality Bin
  edm::Handle<int> centralityBinHandle;
  iEvent.getByToken(_centralityBinToken, centralityBinHandle);
  int hiBin = -1;
  if (centralityBinHandle.isValid()) { hiBin = *centralityBinHandle; }

  timeEvt_.Fill(iEvent, primaryVertexCollection, hiHF, hiBin);

  // Extract the track information
  edm::Handle< edm::View< reco::Track > > recoTracksHandle;
  iEvent.getByToken(_recoTracksToken, recoTracksHandle);
  if (!recoTracksHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track collection is not valid!" << std::endl; return; }
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
  edm::Handle< edm::ValueMap< float > > tofPIDT0Handle;
  iEvent.getByToken(_tofPIDT0Token, tofPIDT0Handle);
  if (!tofPIDT0Handle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "TOF PID t0 valueMap is not valid!" << std::endl; return; }
  //
  edm::Handle< edm::ValueMap< float > > tofPIDSigmaT0Handle;
  iEvent.getByToken(_tofPIDSigmaT0Token, tofPIDSigmaT0Handle);
  if (!tofPIDSigmaT0Handle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "TOF PID sigma t0 valueMap is not valid!" << std::endl; return; }
  //
  // Extract the track dEdx information
  edm::Handle< edm::ValueMap< reco::DeDxData > > dEdxHandle;
  iEvent.getByToken(_dEdxToken, dEdxHandle);
  if (!dEdxHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "dEdx valueMap is not valid!" << std::endl; return; }
  //
  // Extract the gen information
  edm::Handle< reco::GenParticleCollection > genParticlesHandle;
  iEvent.getByToken(_genParticlesToken, genParticlesHandle);
  //
  edm::Handle< float > genT0Handle;
  iEvent.getByToken(_genT0Token, genT0Handle);
  //
  // Extract the track builder
  edm::ESHandle< TransientTrackBuilder > theTTBuilderHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTBuilderHandle);
  if (!theTTBuilderHandle.isValid()) { edm::LogError("TimeAnalyzer_RecordNotValid") << "Transient track builder record is not valid!" << std::endl; return; }
  //
  timeEvt_.Fill(*recoTracksHandle, primaryVertexCollection, *trackBetaHandle, *trackT0Handle,
                *trackSigmaT0Handle, *trackTMTDHandle, *trackSigmaTMTDHandle, *trackMomHandle,
                *trackPathLengthHandle, *tofPIDT0Handle, *tofPIDSigmaT0Handle, *dEdxHandle,
                genParticlesHandle, genT0Handle, *theTTBuilderHandle);

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
reco::GenParticleRef
TimeEvent::findMother(const reco::GenParticleRef& genParRef)
{
  reco::GenParticleRef genMomRef = genParRef;
  int pdg = genParRef->pdgId(); const int pdg_OLD = pdg;
  while(pdg==pdg_OLD && genMomRef->numberOfMothers()>0) {
    genMomRef = genMomRef->motherRef(0);
    pdg = genMomRef->pdgId();
  }
  return ( (pdg_OLD==pdg) ? genParRef : genMomRef );
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Fill(const edm::View<reco::Track>& recoTracks,    const reco::VertexCollection& vtxCol,
                const edm::ValueMap<float>& trackBeta,       const edm::ValueMap<float>& trackT0,
                const edm::ValueMap<float>& trackSigmaT0,    const edm::ValueMap<float>& trackTMTD,
                const edm::ValueMap<float>& trackSigmaTMTD,  const edm::ValueMap<float>& trackMom,
                const edm::ValueMap<float>& trackPathLength, const edm::ValueMap<float>& tofPIDT0,
                const edm::ValueMap<float>& tofPIDSigmaT0,   const edm::ValueMap<reco::DeDxData>& dEdxData,
                const edm::Handle<reco::GenParticleCollection>& genParticlesHandle,
                const edm::Handle<float>& genT0Handle,       const TransientTrackBuilder& theTTBuilder)
{
  // Select the vertex
  unsigned int vtxIdx=0;
  //for (unsigned int iVtx = 0; iVtx < vtxCol.size(); ++iVtx) { if (vtxCol.at(iVtx).tracksSize() > vtxCol.at(vtxIdx).tracksSize()) { vtxIdx = iVtx; } }
  const reco::Vertex& priVtx = vtxCol.at(vtxIdx);
  // Fill the track information
  this->Reco_Track_N = 0;
  this->Reco_DiTrack_N = 0;
  this->Reco_TriTrack_N = 0;
  std::vector< int > recoGenIdx;
  reco::TrackCollection recoTracksColl;
  for (uint i = 0; i < recoTracks.size(); i++) {
    auto obj = recoTracks.ptrAt(i);
    // Select tracks with p > 700 MeV
    //if (obj->p() <= 0.7) continue;
    // Select tracks with pTError < 10% of pt
    if ((obj->ptError()/obj->pt()) >= 0.1) continue;
    // Select tracks passing high purity
    if (obj->quality(reco::TrackBase::TrackQuality::highPurity) == false) continue;
    // Select tracks with more than or equal 11 hits
    if (obj->numberOfValidHits() < 11) continue;
    // Select tracks fitted with normalized chi2 lower than 7
    if (obj->normalizedChi2() >= 7.) continue;
    // Compute the distance to primary vertex
    const double vx = priVtx.x();
    const double vy = priVtx.y();
    const double vz = priVtx.z();
    const double vxError = priVtx.xError();
    const double vyError = priVtx.yError();
    const double vzError = priVtx.zError();
    const math::XYZPoint vtx(vx, vy, vz);
    const double dz = obj->dz(vtx);
    const double dzsigma = std::sqrt(std::pow(obj->dzError(), 2.0) + std::pow(vzError, 2.0));
    const double dcaDz = dz/dzsigma;
    const double dxy = obj->dxy(vtx);
    const double dxysigma = std::sqrt(std::pow(obj->dxyError(), 2.0) + vxError*vyError);
    const double dcaDxy = dxy/dxysigma;
    // Select tracks with DCA less than 3
    //if (std::abs(dcaDz) >= 3.0 || std::abs(dcaDxy) >= 3.0) continue;
    // Extract the reco information
    const double beta_PION  = trackBeta[obj];
    const double t0_PION    = trackT0[obj];
    const double t0Err_PION = trackSigmaT0[obj];
    const double tMTD       = trackTMTD[obj];
    const double tMTDErr    = trackSigmaTMTD[obj];
    const double magp       = trackMom[obj];
    const double pathLength = trackPathLength[obj];
    const double t0_PID     = tofPIDT0[obj];
    const double t0Err_PID  = tofPIDSigmaT0[obj];
    const double t0_PV      = priVtx.t();
    const double t0Err_PV   = priVtx.tError();
    const double pt         = obj->pt();
    const double ptErr      = obj->ptError();
    const double eta        = obj->eta();
    const double phi        = obj->phi();
    const double charge     = obj->charge();
    const double normChi2   = obj->normalizedChi2();
    const int    quality    = obj->qualityMask();
    const int    nHits      = obj->numberOfValidHits();
    // Compute beta and mass using primary vertex info
    const double dt_PV = tMTD - t0_PV;
    const double beta_PV = (tMTDErr>=0. ? (pathLength/dt_PV)*(1./c_cm_ns) : -99.);
    const double gammaSq_PV = (tMTDErr>=0. ? 1./(1.0-std::pow(beta_PV, 2.0)) : -99.);
    const double massSq_PV  = (tMTDErr>=0. ? std::pow(magp, 2.0)/(gammaSq_PV - 1.) : -99.);
    // Compute beta and mass using PID info
    const double dt_PID = tMTD - t0_PID;
    const double beta_PID = (tMTDErr>=0. ? (pathLength/dt_PID)*(1./c_cm_ns) : -99.);
    const double gammaSq_PID = (tMTDErr>=0. ? 1./(1.0-std::pow(beta_PID, 2.0)) : -99.);
    const double massSq_PID  = (tMTDErr>=0. ? std::pow(magp, 2.0)/(gammaSq_PID - 1.) : -99.);
    // Store reco information
    this->Reco_Track_N++;
    this->Reco_Track_beta_PI.push_back(beta_PION);
    this->Reco_Track_t0_PI.push_back(t0_PION);
    this->Reco_Track_t0Err_PI.push_back(t0Err_PION);
    this->Reco_Track_beta_PV.push_back(beta_PV);
    this->Reco_Track_t0_PV.push_back(t0_PV);
    this->Reco_Track_t0Err_PV.push_back(t0Err_PV);
    this->Reco_Track_massSq_PV.push_back(massSq_PV);
    this->Reco_Track_beta_PID.push_back(beta_PID);
    this->Reco_Track_t0_PID.push_back(t0_PID);
    this->Reco_Track_t0Err_PID.push_back(t0Err_PID);
    this->Reco_Track_massSq_PID.push_back(massSq_PID);
    this->Reco_Track_pathLength.push_back(pathLength);
    this->Reco_Track_p.push_back(magp);
    this->Reco_Track_tMTD.push_back(tMTD);
    this->Reco_Track_tMTDErr.push_back(tMTDErr);
    this->Reco_Track_eta.push_back(eta);
    this->Reco_Track_phi.push_back(phi);
    this->Reco_Track_charge.push_back(charge);
    this->Reco_Track_quality.push_back(quality);
    this->Reco_Track_pt.push_back(pt);
    this->Reco_Track_ptErr.push_back(ptErr);
    this->Reco_Track_dcaDz.push_back(dcaDz);
    this->Reco_Track_dcaDxy.push_back(dcaDxy);
    this->Reco_Track_nHits.push_back(nHits);
    this->Reco_Track_normChi2.push_back(normChi2);
    // Store reco dEdx information
    this->Reco_Track_dEdx.push_back(dEdxData[obj].dEdx());
    this->Reco_Track_dEdxErr.push_back(dEdxData[obj].dEdxError());
    this->Reco_Track_nSatMea.push_back(dEdxData[obj].numberOfSaturatedMeasurements());
    this->Reco_Track_nMea.push_back(dEdxData[obj].numberOfMeasurements());
    // Keep the track
    recoTracksColl.push_back(*obj);
    //
    // Do the gen matching
    // 
    if (genParticlesHandle.isValid()) {
      int genIdx = -1;
      double matchDeltaR = 999999999.;
      TLorentzVector recMom; recMom.SetPtEtaPhiM(obj->pt(), obj->eta(), obj->phi(), 0.0);
      const auto& genParticles = *genParticlesHandle;
      for (uint iGen = 0; iGen < genParticles.size(); iGen++) {
        const auto& gen = genParticles.at(iGen);
        if (gen.status()==1 && gen.charge()!=0 && std::abs(gen.eta())<4.3 && gen.p()>0.5) {
          // Lorentz vector
          TLorentzVector genMom; genMom.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
          const double deltaR = genMom.DeltaR(recMom);
          // Check if reco and gen match
          if (
              (deltaR < 0.03) &&
              ((std::abs(pt - gen.pt())/gen.pt()) < 0.5) &&
              (deltaR < matchDeltaR)
              )
            {
              matchDeltaR = deltaR;
              genIdx = iGen;
            }
        }
      }
      // Extract the gen information
      double gen_pdgId=-999. , gen_pt=-1. , gen_eta=-9. , gen_phi=-9. , gen_t0=-9. , gen_charge=-9. , gen_mass=-1. , gen_momPdgId=-999.;
      if (genIdx > -1) {
        const auto& gen = reco::GenParticleRef(genParticlesHandle, genIdx);
        gen_pdgId  = gen->pdgId();
        gen_pt     = gen->pt();
        gen_eta    = gen->eta();
        gen_phi    = gen->phi();
        gen_t0     = (genT0Handle.isValid() ? *genT0Handle : -9.);
        gen_charge = gen->charge();
        gen_mass   = gen->mass();
        gen_momPdgId = findMother(gen)->pdgId();
      }
      // Store gen information
      this->Gen_Track_pdgId.push_back(gen_pdgId);
      this->Gen_Track_pt.push_back(gen_pt);
      this->Gen_Track_eta.push_back(gen_eta);
      this->Gen_Track_phi.push_back(gen_phi);
      this->Gen_Track_t0.push_back(gen_t0);
      this->Gen_Track_charge.push_back(gen_charge);
      this->Gen_Track_mass.push_back(gen_mass);
      this->Gen_Track_momPdgId.push_back(gen_momPdgId);
      recoGenIdx.push_back(genIdx);
    }
  }
  if (doDiTracks_) {
    // Loop over selected tracks
    for (uint i = 0; i < recoTracksColl.size(); i++) {
      const auto& piTrk = recoTracksColl.at(i);
      TLorentzVector piMom; piMom.SetPtEtaPhiM(piTrk.pt(), piTrk.eta(), piTrk.phi(), MASS_PION);
      //
      // Reconstruct D0 mesons
      //
      for (uint j = 0; j < recoTracksColl.size(); j++) {
        if (j==i) continue;
        const auto& kaTrk = recoTracksColl.at(j);
        // Select opposite sign tracks between pion and kaon
        if (piTrk.charge()!=kaTrk.charge()) {
          TLorentzVector kaMom; kaMom.SetPtEtaPhiM(kaTrk.pt(), kaTrk.eta(), kaTrk.phi(), MASS_KAON);
          // Compute lorentz vector
          const TLorentzVector diTrackMom = piMom + kaMom;
          // Apply the D0 mass cut
          if (std::abs(diTrackMom.M() - MASS_D0) < 0.15) {
            // Initialize topological variables
            float vProb_PiKa=-1., alpha_PiKa=99., d0Sig_PiKa=-99.;
            // compute pion+kaon vertex fit
            std::vector<reco::TransientTrack> t_tks;
            t_tks.push_back( theTTBuilder.build( piTrk ) );
            t_tks.push_back( theTTBuilder.build( kaTrk ) );
            KalmanVertexFitter vtxFitter = KalmanVertexFitter(true);
            const TransientVertex& diTrackVertex = vtxFitter.vertex( t_tks );
            if (diTrackVertex.isValid()) {
              const float vChi2 = diTrackVertex.totalChiSquared();
              const float vNDF  = diTrackVertex.degreesOfFreedom();
              const TVector3 diTrackVtx( diTrackVertex.position().x() , diTrackVertex.position().y() , diTrackVertex.position().z() );
              const TVector3 vDiff = diTrackVtx - TVector3(priVtx.position().x(), priVtx.position().y(), priVtx.position().z());
              const AlgebraicSymMatrix33 vErr3D = GlobalError(priVtx.error()).matrix() + GlobalError((reco::Vertex(diTrackVertex)).error()).matrix();
              const AlgebraicVector3 vDiff3D( vDiff.x() , vDiff.y() , vDiff.z() );
              //
              vProb_PiKa = TMath::Prob( vChi2 , (int)vNDF );
              alpha_PiKa = vDiff.Angle(diTrackMom.Vect());
              d0Sig_PiKa = vDiff.Mag2() / std::sqrt(ROOT::Math::Similarity(vErr3D, vDiff3D));
            }
            else continue;
            // Store pion+kaon information
            // Official: d0Sig > 3.5 , vProb > 0.08 , alpha < 0.15
            // V0Analyzer: d0Sig > 0. , vProb > 0.0001 , alpha < 999.0
            if ((d0Sig_PiKa > 0.) && (vProb_PiKa > 0.0001) && (alpha_PiKa < 999.)) {
              // Store pion+kaon reconstruction variables
              this->Reco_DiTrack_N++;
              this->Reco_DiTrack_Pion_idx.push_back(i);
              this->Reco_DiTrack_Kaon_idx.push_back(j);
              this->Reco_DiTrack_pt.push_back(diTrackMom.Pt());
              this->Reco_DiTrack_rap.push_back(diTrackMom.Rapidity());
              this->Reco_DiTrack_mass.push_back(diTrackMom.M());
              this->Reco_DiTrack_vProb.push_back(vProb_PiKa);
              this->Reco_DiTrack_alpha.push_back(alpha_PiKa);
              this->Reco_DiTrack_d0Sig.push_back(d0Sig_PiKa);
              // Store pion+kaon generator variables
              if (genParticlesHandle.isValid()) {
                // Initialize the gen information
                int isSwap = -1;
                double genMom_pdgId=-999. , genMom_pt=-1. , genMom_eta=-9. , genMom_phi=-9. , genMom_charge=-9. , genMom_mass=-1.;
                const auto& genPionIdx = recoGenIdx[i];
                const auto& genKaonIdx = recoGenIdx[j];
                if (genPionIdx>=0 && genKaonIdx>=0) {
                  const auto& genPion    = reco::GenParticleRef(genParticlesHandle, genPionIdx);
                  const auto& genKaon    = reco::GenParticleRef(genParticlesHandle, genKaonIdx);
                  // Check if it is swapped
                  if      (std::abs(genPion->pdgId())==211 && std::abs(genKaon->pdgId())==321) { isSwap = 0; }
                  else if (std::abs(genPion->pdgId())==321 && std::abs(genKaon->pdgId())==211) { isSwap = 1; }
                  else    { isSwap = -1; }
                  // Find mother of both tracks
                  const auto& genPionMom = findMother(genPion);
                  const auto& genKaonMom = findMother(genKaon);
                  if (genPionMom == genKaonMom) {
                    genMom_pdgId  = genPionMom->pdgId();
                    genMom_pt     = genPionMom->pt();
                    genMom_eta    = genPionMom->eta();
                    genMom_phi    = genPionMom->phi();
                    genMom_charge = genPionMom->charge();
                    genMom_mass   = genPionMom->mass();
                  }
                }
                // Store gen information
                this->Gen_DiTrack_pdgId.push_back(genMom_pdgId);
                this->Gen_DiTrack_pt.push_back(genMom_pt);
                this->Gen_DiTrack_eta.push_back(genMom_eta);
                this->Gen_DiTrack_phi.push_back(genMom_phi);
                this->Gen_DiTrack_charge.push_back(genMom_charge);
                this->Gen_DiTrack_mass.push_back(genMom_mass);
                this->Gen_DiTrack_isSwap.push_back(isSwap);
              }
            }
          }
          /*
          //
          // Reconstruct LambdaC mesons
          //
          for (uint k = 0; k < recoTracksColl.size(); k++) {
          if (k==i || k==j) continue;
          const auto& prTrk = recoTracksColl.at(k);
          // Select opposite sign tracks between proton and kaon
          if (prTrk.charge()!=kaTrk.charge() && prTrk.charge()>0) {
          TLorentzVector prMom; prMom.SetPtEtaPhiM(prTrk.pt(), prTrk.eta(), prTrk.phi(), MASS_PROT);
          // Compute lorentz vector
          const TLorentzVector triTrackMom = prMom + diTrackMom;
          // Apply invariant mass cuts
          if (std::abs(triTrackMom.M() - MASS_LC) > 0.15) continue;
          // Apply cuts on daughter pt
          if ((piTrk.pt()/triTrackMom.Pt()) <= 0.12) continue;
          if ((kaTrk.pt()/triTrackMom.Pt()) <= 0.14) continue;
          //if ((prTrk.pt()/triTrackMom.Pt()) <= 0.28) continue;
          // Initialize topological variables
          float vProb_PrPiKa=-1., alpha_PrPiKa=99., d0Sig_PrPiKa=-99.;
          // compute proton+pion+kaon vertex fit
          std::vector<reco::TransientTrack> t_3tks;
          t_3tks.push_back( theTTBuilder.build( piTrk ) );
          t_3tks.push_back( theTTBuilder.build( kaTrk ) );
          t_3tks.push_back( theTTBuilder.build( prTrk ) );
          KalmanVertexFitter vtxFitter3Trk = KalmanVertexFitter(true);
          const TransientVertex& triTrackVertex = vtxFitter3Trk.vertex( t_3tks );
          if (triTrackVertex.isValid()) {
          const float vChi2 = triTrackVertex.totalChiSquared();
          const float vNDF  = triTrackVertex.degreesOfFreedom();
          const TVector3 triTrackVtx( triTrackVertex.position().x() , triTrackVertex.position().y() , triTrackVertex.position().z() );
          const TVector3 vDiff = triTrackVtx - TVector3(priVtx.position().x(), priVtx.position().y(), priVtx.position().z());
          const AlgebraicSymMatrix33 vErr3D = GlobalError(priVtx.error()).matrix() + GlobalError((reco::Vertex(triTrackVertex)).error()).matrix();
          const AlgebraicVector3 vDiff3D( vDiff.x() , vDiff.y() , vDiff.z() );
          //
          vProb_PrPiKa = TMath::Prob( vChi2 , (int)vNDF );
          alpha_PrPiKa = vDiff.Angle(triTrackMom.Vect());
          d0Sig_PrPiKa = vDiff.Mag2() / std::sqrt(ROOT::Math::Similarity(vErr3D, vDiff3D));
          }
          else continue;
          // Store pion+kaon information
          // Official: d0Sig > 3.75 , vProb > 0.2 , alpha < 0.1
          if ((d0Sig_PrPiKa > 3.5) && (vProb_PrPiKa > 0.08) && (alpha_PrPiKa < 0.15)) {
          this->Reco_TriTrack_N++;
          this->Reco_TriTrack_Pion_idx.push_back(i);
          this->Reco_TriTrack_Kaon_idx.push_back(j);
          this->Reco_TriTrack_Proton_idx.push_back(k);
          this->Reco_TriTrack_pt.push_back(triTrackMom.Pt());
          this->Reco_TriTrack_rap.push_back(triTrackMom.Rapidity());
          this->Reco_TriTrack_mass.push_back(triTrackMom.M());
          this->Reco_TriTrack_vProb.push_back(vProb_PrPiKa);
          this->Reco_TriTrack_alpha.push_back(alpha_PrPiKa);
          this->Reco_TriTrack_d0Sig.push_back(d0Sig_PrPiKa);
          }
          }
          }
          */
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Fill(const edm::Event& iEvent, const reco::VertexCollection& vtxCol, const float& hiHF, const int& hiBin)
{
  this->Event_nRun   = iEvent.id().run();
  this->Event_nLumi  = iEvent.luminosityBlock();
  this->Event_nBX    = std::max(iEvent.bunchCrossing(),0);
  this->Event_nOrbit = std::max(iEvent.orbitNumber(),0);
  this->Event_nEvent = iEvent.id().event();
  // Primary Vertex Information
  unsigned int vtxIdx=0;
  //for (unsigned int iVtx = 0; iVtx < vtxCol.size(); ++iVtx) { if (vtxCol.at(iVtx).tracksSize() > vtxCol.at(vtxIdx).tracksSize()) { vtxIdx = iVtx; } }
  const reco::Vertex& priVtx  = vtxCol.at(vtxIdx);
  this->Event_PriVtx_Position = TVector3( priVtx.x()      , priVtx.y()      , priVtx.z()      );
  this->Event_PriVtx_Error    = TVector3( priVtx.xError() , priVtx.yError() , priVtx.zError() );
  this->Event_nPV = vtxCol.size();
  // Centrality Information
  this->Event_hiHF   = hiHF;
  this->Event_hiBin  = hiBin;
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
  tree_->Branch("Event_hiHF",            &(this->Event_hiHF),   "Event_hiHF/F"   );
  tree_->Branch("Event_hiBin",           &(this->Event_hiBin),  "Event_hiBin/I"  );
  // Track Info
  tree_->Branch("Reco_Track_N",          &(this->Reco_Track_N), "Reco_Track_N/i" );
  tree_->Branch("Reco_Track_beta_PI",    &(this->Reco_Track_beta_PI)             );
  tree_->Branch("Reco_Track_t0_PI",      &(this->Reco_Track_t0_PI)               );
  tree_->Branch("Reco_Track_t0Err_PI",   &(this->Reco_Track_t0Err_PI)            );
  tree_->Branch("Reco_Track_beta_PV",    &(this->Reco_Track_beta_PV)             );
  tree_->Branch("Reco_Track_t0_PV",      &(this->Reco_Track_t0_PV)               );
  tree_->Branch("Reco_Track_t0Err_PV",   &(this->Reco_Track_t0Err_PV)            );
  tree_->Branch("Reco_Track_massSq_PV",  &(this->Reco_Track_massSq_PV)           );
  tree_->Branch("Reco_Track_beta_PID",   &(this->Reco_Track_beta_PID)            );
  tree_->Branch("Reco_Track_t0_PID",     &(this->Reco_Track_t0_PID)              );
  tree_->Branch("Reco_Track_t0Err_PID",  &(this->Reco_Track_t0Err_PID)           );
  tree_->Branch("Reco_Track_massSq_PID", &(this->Reco_Track_massSq_PID)          );
  tree_->Branch("Reco_Track_pathLength", &(this->Reco_Track_pathLength)          );
  tree_->Branch("Reco_Track_p",          &(this->Reco_Track_p)                   );
  tree_->Branch("Reco_Track_tMTD",       &(this->Reco_Track_tMTD)                );
  tree_->Branch("Reco_Track_tMTDErr",    &(this->Reco_Track_tMTDErr)             );
  tree_->Branch("Reco_Track_eta",        &(this->Reco_Track_eta)                 );
  tree_->Branch("Reco_Track_phi",        &(this->Reco_Track_phi)                 );
  tree_->Branch("Reco_Track_charge",     &(this->Reco_Track_charge)              );
  tree_->Branch("Reco_Track_quality",    &(this->Reco_Track_quality)             );
  tree_->Branch("Reco_Track_pt",         &(this->Reco_Track_pt)                  );
  tree_->Branch("Reco_Track_ptErr",      &(this->Reco_Track_ptErr)               );
  tree_->Branch("Reco_Track_dcaDz",      &(this->Reco_Track_dcaDz)               );
  tree_->Branch("Reco_Track_dcaDxy",     &(this->Reco_Track_dcaDxy)              );
  tree_->Branch("Reco_Track_normChi2",   &(this->Reco_Track_normChi2)            );
  tree_->Branch("Reco_Track_nHits",      &(this->Reco_Track_nHits)               );
  // Track dEdx Info
  tree_->Branch("Reco_Track_dEdx",       &(this->Reco_Track_dEdx)                );
  tree_->Branch("Reco_Track_dEdxErr",    &(this->Reco_Track_dEdxErr)             );
  tree_->Branch("Reco_Track_nSatMea",    &(this->Reco_Track_nSatMea)             );
  tree_->Branch("Reco_Track_nMea",       &(this->Reco_Track_nMea)                );
  // Gen Particle Info
  tree_->Branch("Gen_Track_pdgId",       &(this->Gen_Track_pdgId)                );
  tree_->Branch("Gen_Track_pt",          &(this->Gen_Track_pt)                   );
  tree_->Branch("Gen_Track_eta",         &(this->Gen_Track_eta)                  );
  tree_->Branch("Gen_Track_phi",         &(this->Gen_Track_phi)                  );
  tree_->Branch("Gen_Track_t0",          &(this->Gen_Track_t0)                   );
  tree_->Branch("Gen_Track_charge",      &(this->Gen_Track_charge)               );
  tree_->Branch("Gen_Track_mass",        &(this->Gen_Track_mass)                 );
  tree_->Branch("Gen_Track_momPdgId",    &(this->Gen_Track_momPdgId)             );
  // Pion+Kaon Info
  tree_->Branch("Reco_DiTrack_N",        &(this->Reco_DiTrack_N), "Reco_DiTrack_N/i");
  tree_->Branch("Reco_DiTrack_Pion_idx", &(this->Reco_DiTrack_Pion_idx)          );
  tree_->Branch("Reco_DiTrack_Kaon_idx", &(this->Reco_DiTrack_Kaon_idx)          );
  tree_->Branch("Reco_DiTrack_pt",       &(this->Reco_DiTrack_pt)                );
  tree_->Branch("Reco_DiTrack_rap",      &(this->Reco_DiTrack_rap)               );
  tree_->Branch("Reco_DiTrack_mass",     &(this->Reco_DiTrack_mass)              );
  tree_->Branch("Reco_DiTrack_vProb",    &(this->Reco_DiTrack_vProb)             );
  tree_->Branch("Reco_DiTrack_alpha",    &(this->Reco_DiTrack_alpha)             );
  tree_->Branch("Reco_DiTrack_d0Sig",    &(this->Reco_DiTrack_d0Sig)             );
  // Gen Pion+Kaon Info
  tree_->Branch("Gen_DiTrack_pdgId",     &(this->Gen_DiTrack_pdgId)              );
  tree_->Branch("Gen_DiTrack_pt",        &(this->Gen_DiTrack_pt)                 );
  tree_->Branch("Gen_DiTrack_eta",       &(this->Gen_DiTrack_eta)                );
  tree_->Branch("Gen_DiTrack_phi",       &(this->Gen_DiTrack_phi)                );
  tree_->Branch("Gen_DiTrack_charge",    &(this->Gen_DiTrack_charge)             );
  tree_->Branch("Gen_DiTrack_mass",      &(this->Gen_DiTrack_mass)               );
  tree_->Branch("Gen_DiTrack_isSwap",    &(this->Gen_DiTrack_isSwap)             );
  // Proton+Pion+Kaon Info
  tree_->Branch("Reco_TriTrack_N",       &(this->Reco_TriTrack_N), "Reco_TriTrack_N/i");
  tree_->Branch("Reco_TriTrack_Pion_idx",&(this->Reco_TriTrack_Pion_idx)         );
  tree_->Branch("Reco_TriTrack_Kaon_idx",&(this->Reco_TriTrack_Kaon_idx)         );
  tree_->Branch("Reco_TriTrack_Proton_idx",&(this->Reco_TriTrack_Proton_idx)     );
  tree_->Branch("Reco_TriTrack_pt",      &(this->Reco_TriTrack_pt)               );
  tree_->Branch("Reco_TriTrack_rap",     &(this->Reco_TriTrack_rap)              );
  tree_->Branch("Reco_TriTrack_mass",    &(this->Reco_TriTrack_mass)             );
  tree_->Branch("Reco_TriTrack_vProb",   &(this->Reco_TriTrack_vProb)            );
  tree_->Branch("Reco_TriTrack_alpha",   &(this->Reco_TriTrack_alpha)            );
  tree_->Branch("Reco_TriTrack_d0Sig",   &(this->Reco_TriTrack_d0Sig)            );
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
  this->Event_hiHF   = -1.;
  this->Event_hiBin  = -1;
  // Track Info
  this->Reco_Track_N = 0;
  this->Reco_Track_beta_PI.clear();
  this->Reco_Track_t0_PI.clear();
  this->Reco_Track_t0Err_PI.clear();
  this->Reco_Track_beta_PV.clear();
  this->Reco_Track_t0_PV.clear();
  this->Reco_Track_t0Err_PV.clear();
  this->Reco_Track_massSq_PV.clear();
  this->Reco_Track_beta_PID.clear();
  this->Reco_Track_t0_PID.clear();
  this->Reco_Track_t0Err_PID.clear();
  this->Reco_Track_massSq_PID.clear();
  this->Reco_Track_pathLength.clear();
  this->Reco_Track_p.clear();
  this->Reco_Track_tMTD.clear();
  this->Reco_Track_tMTDErr.clear();
  this->Reco_Track_eta.clear();
  this->Reco_Track_phi.clear();
  this->Reco_Track_charge.clear();
  this->Reco_Track_quality.clear();
  this->Reco_Track_pt.clear();
  this->Reco_Track_ptErr.clear();
  this->Reco_Track_dcaDz.clear();
  this->Reco_Track_dcaDxy.clear();
  this->Reco_Track_normChi2.clear();
  this->Reco_Track_nHits.clear();
  // Track dEdx Info
  this->Reco_Track_dEdx.clear();
  this->Reco_Track_dEdxErr.clear();
  this->Reco_Track_nSatMea.clear();
  this->Reco_Track_nMea.clear();
  // Gen Particle Info
  this->Gen_Track_pdgId.clear();
  this->Gen_Track_pt.clear();
  this->Gen_Track_eta.clear();
  this->Gen_Track_phi.clear();
  this->Gen_Track_t0.clear();
  this->Gen_Track_charge.clear();
  this->Gen_Track_mass.clear();
  this->Gen_Track_momPdgId.clear();
  // Pion+Kaon Info
  this->Reco_DiTrack_N = 0;
  this->Reco_DiTrack_Pion_idx.clear();
  this->Reco_DiTrack_Kaon_idx.clear();
  this->Reco_DiTrack_pt.clear();
  this->Reco_DiTrack_rap.clear();
  this->Reco_DiTrack_mass.clear();
  this->Reco_DiTrack_vProb.clear();
  this->Reco_DiTrack_alpha.clear();
  this->Reco_DiTrack_d0Sig.clear();
  // Gen Pion+Kaon Info
  this->Gen_DiTrack_pdgId.clear();
  this->Gen_DiTrack_pt.clear();
  this->Gen_DiTrack_eta.clear();
  this->Gen_DiTrack_phi.clear();
  this->Gen_DiTrack_charge.clear();
  this->Gen_DiTrack_mass.clear();
  this->Gen_DiTrack_isSwap.clear();
  // Proton+Pion+Kaon Info
  this->Reco_TriTrack_N = 0;
  this->Reco_TriTrack_Pion_idx.clear();
  this->Reco_TriTrack_Kaon_idx.clear();
  this->Reco_TriTrack_Proton_idx.clear();
  this->Reco_TriTrack_pt.clear();
  this->Reco_TriTrack_rap.clear();
  this->Reco_TriTrack_mass.clear();
  this->Reco_TriTrack_vProb.clear();
  this->Reco_TriTrack_alpha.clear();
  this->Reco_TriTrack_d0Sig.clear();
}

DEFINE_FWK_MODULE(TimeAnalyzer);
