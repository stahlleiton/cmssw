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
#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

//Headers for services and tools
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/TrackAssociation/interface/trackAssociationChi2.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//system include files
#include <TTree.h>
#include <TLorentzVector.h>
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

  void   Clear ( void );
  void   Init  ( const size_t& );

  void   SetTree      ( TTree* p ) { tree_ = p; }
  void   SetBranches  ( void );

  void   Fill         ( const edm::Handle<reco::TrackCollection>&, const std::vector<reco::VertexRef>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const edm::ValueMap<float>&,
                        const edm::ValueMap<float>&,   const std::vector<edm::ValueMap<reco::DeDxData> >&,
                        const edm::Handle<reco::GenParticleCollection>& genParticlesHandle,
                        const edm::Handle<reco::RecoToSimCollection>&,
                        const reco::BeamSpot&, const MagneticField& );
  void   Fill         ( const edm::Event&, const std::vector<reco::VertexRef>&, const std::map<std::string, std::map<reco::VertexRef, size_t> >&, const float&, const int&, const float& );

 private:

  reco::GenParticleRef findMother(const reco::GenParticleRef&);

  const float c_cm_ns   = 2.99792458e1; //[cm/ns]
  const std::vector<float> MASS = {0.13957018, 0.493677, 0.9382720813, 1.87561, 2.80925, 2.80923/2., 3.72742/2.};

  TTree*                  tree_;

  // Event Info
  UInt_t                  Event_nRun;
  UShort_t                Event_nLumi;
  UInt_t                  Event_nBX;
  ULong64_t               Event_nOrbit;
  ULong64_t               Event_nEvent;
  // Centrality Information
  Float_t                 Event_hiHF;
  Short_t                 Event_hiBin;
  // Generator Information
  Float_t                 Event_genT0;
  // Vertex Info
  std::vector< Float_t >  Reco_Vertex_x;
  std::vector< Float_t >  Reco_Vertex_y;
  std::vector< Float_t >  Reco_Vertex_z;
  std::vector< Float_t >  Reco_Vertex_t0;
  std::vector< Float_t >  Reco_Vertex_t0Err;
 public:
  std::map<std::string, std::vector<Int_t> >  Reco_Vertex_nTrk;
 private:
  // Track Info
  std::vector< Int_t   >  Reco_Track_vtxIdx;
  std::vector< Float_t >  Reco_Track_beta_PID;
  std::vector< Float_t >  Reco_Track_t0_PID;
  std::vector< Float_t >  Reco_Track_t0Err_PID;
  std::vector< Float_t >  Reco_Track_beta_PV;
  std::vector< Float_t >  Reco_Track_t0_TOF;
  std::vector< Float_t >  Reco_Track_t0Err_TOF;
  std::vector< Float_t >  Reco_Track_pathLength;
  std::vector< Float_t >  Reco_Track_p;
  std::vector< Float_t >  Reco_Track_tMTD;
  std::vector< Float_t >  Reco_Track_tMTDErr;
  std::vector< Float_t >  Reco_Track_eta;
  std::vector< Float_t >  Reco_Track_phi;
  std::vector< Short_t >  Reco_Track_charge;
  std::vector< Int_t   >  Reco_Track_quality;
  std::vector< Float_t >  Reco_Track_pt;
  std::vector< Float_t >  Reco_Track_ptErr;
  std::vector< Float_t >  Reco_Track_dcaDz;
  std::vector< Float_t >  Reco_Track_dcaDxy;
  std::vector< Float_t >  Reco_Track_normChi2;
  std::vector< Int_t   >  Reco_Track_nHits;
  // Track dEdx
  std::vector< std::vector< Float_t > > Reco_Track_dEdx;
  std::vector< Int_t   >  Reco_Track_genIdx;
  std::vector< Float_t >  Reco_Track_genChi2;
  // Gen Particle Info
  std::vector< Int_t   >  Gen_Track_recoIdx;
  std::vector< Int_t   >  Gen_Track_pdgId;
  std::vector< Float_t >  Gen_Track_pt;
  std::vector< Float_t >  Gen_Track_eta;
  std::vector< Float_t >  Gen_Track_phi;
  std::vector< Short_t >  Gen_Track_charge;
  std::vector< Float_t >  Gen_Track_mass;
  std::vector< Int_t   >  Gen_Track_momPdgId;
  std::vector< Bool_t  >  Gen_Track_isGen;
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
  const edm::EDGetTokenT< reco::TrackCollection  > _recoTracksToken;
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
  const edm::EDGetTokenT< reco::GenParticleCollection > _genParticlesToken;
  const edm::EDGetTokenT< float                  > _genT0Token;
  const edm::EDGetTokenT< reco::Centrality       > _centralitySrcToken;
  const edm::EDGetTokenT< int                    > _centralityBinToken;
  const edm::EDGetTokenT< reco::RecoToSimCollection > _associatorMapRS;
  std::vector< edm::EDGetTokenT< edm::ValueMap< reco::DeDxData > > > _dEdxTokens;

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
  _recoTracksToken      ( consumes< reco::TrackCollection      >( iConfig.getParameter< edm::InputTag >( "recoTracksTag"      ) )),
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
  _associatorMapRS      ( consumes< reco::RecoToSimCollection  >( iConfig.getParameter< edm::InputTag >( "associatorMap"      ) ))
{
  for (const auto& p : iConfig.getParameter< std::vector<edm::InputTag> >("dEdxTags"))
    _dEdxTokens.push_back(consumes< edm::ValueMap< reco::DeDxData > >(p));
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
    timeEvt_.Init(_dEdxTokens.size());
    timeEvt_.SetTree( timeTree_ );
    timeEvt_.SetBranches();
  }

  edm::ESHandle<MagneticField> magField;
  iSetup.get<IdealMagneticFieldRecord>().get(magField);

  // Extract the BeamSpot
  edm::Handle< reco::BeamSpot > beamSpotHandle;
  iEvent.getByToken(_beamSpotToken, beamSpotHandle);

  // Extract the Offline Primary Vertices
  edm::Handle< reco::VertexCollection > primaryVertexHandle;
  iEvent.getByToken(_primaryVertexToken, primaryVertexHandle);
  std::vector<reco::VertexRef> primaryVertexCollection;
  if (primaryVertexHandle.isValid() && primaryVertexHandle->size()>0) {
    for (size_t i = 0; i < primaryVertexHandle->size(); i++) {
      const reco::Vertex& primaryVertex = primaryVertexHandle->at(i);
      if (!primaryVertex.isFake() && primaryVertex.tracksSize()>=2) {
        primaryVertexCollection.push_back(reco::VertexRef(primaryVertexHandle, i));
      }
    }
  }

  // Extract the Offline Tracks
  edm::Handle< reco::TrackCollection > recoTracksHandle;
  iEvent.getByToken(_recoTracksToken, recoTracksHandle);
  if (!recoTracksHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "General track collection is not valid!" << std::endl; return; }

  // Derive track multiplicity vertex map
  std::map<std::string, std::map<reco::VertexRef, size_t> > vtxNTrk;
  for (const auto& p : timeEvt_.Reco_Vertex_nTrk) {
    for (const auto& pv : primaryVertexCollection) vtxNTrk[p.first][pv] = 0;
  }
  if (recoTracksHandle.isValid() && !primaryVertexCollection.empty()) {
    for (const auto& trk : *recoTracksHandle) {
      if (!trk.quality(reco::TrackBase::highPurity)) continue;
      const auto absEta = std::abs(trk.eta());
      if (trk.charge() == 0 || absEta >= 3.0) continue;
      if (trk.ptError()/trk.pt() >= 0.1) continue;
      const auto dzTrkErr2 = trk.dzError()*trk.dzError();
      const auto dxyTrkErr2 = trk.dxyError()*trk.dxyError();
      // loop over primary vertices
      for (const auto& pv : primaryVertexCollection) {
        const auto dz = trk.dz(pv->position());
        const auto dzErr2 = dzTrkErr2 + pv->zError()*pv->zError();
        if (dz*dz >= 9.0*dzErr2) continue;
        const auto dxy = trk.dxy(pv->position());
        const auto dxyErr2 = dxyTrkErr2 + pv->xError()*pv->yError();
        if (dxy*dxy >= 9.0*dxyErr2) continue;
        if (trk.pt() > 0.4) vtxNTrk["AEta3"][pv] += 1;
        if (absEta < 0.5) vtxNTrk["AEta0p5"][pv] += 1;
        if (absEta < 0.8) vtxNTrk["AEta0p8"][pv] += 1;
      }
    }
  }

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

  // Extract gen information
  edm::Handle< float > genT0Handle;
  iEvent.getByToken(_genT0Token, genT0Handle);
  const auto genT0 = (genT0Handle.isValid() ? *genT0Handle : -9.);
  edm::Handle< reco::GenParticleCollection > genParticlesHandle;
  iEvent.getByToken(_genParticlesToken, genParticlesHandle);

  timeEvt_.Fill(iEvent, primaryVertexCollection, vtxNTrk, hiHF, hiBin, genT0);

  // Extract the track information
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
  std::vector< edm::ValueMap<reco::DeDxData> > dEdxMaps;
  for (size_t i=0; i<_dEdxTokens.size(); i++) {
    edm::Handle< edm::ValueMap< reco::DeDxData > > dEdxHandle;
    iEvent.getByToken(_dEdxTokens[i], dEdxHandle);
    if (!dEdxHandle.isValid()) { edm::LogError("TimeAnalyzer_ProductNotValid") << "dEdx valueMap is not valid!" << std::endl; return; }
    dEdxMaps.push_back(*dEdxHandle);
  }
  //
  // Extract sim information
  edm::Handle< reco::RecoToSimCollection > recToSimMapHandle;
  iEvent.getByToken(_associatorMapRS, recToSimMapHandle);
  //
  timeEvt_.Fill(recoTracksHandle, primaryVertexCollection, *trackBetaHandle, *trackT0Handle,
                *trackSigmaT0Handle, *trackTMTDHandle, *trackSigmaTMTDHandle, *trackMomHandle,
                *trackPathLengthHandle, *tofPIDT0Handle, *tofPIDSigmaT0Handle, dEdxMaps,
                genParticlesHandle, recToSimMapHandle, *beamSpotHandle, *magField);
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
TimeEvent::Fill(const edm::Handle<reco::TrackCollection>& recoTracks, const std::vector<reco::VertexRef>& vtxCol,
                const edm::ValueMap<float>& trackBeta,       const edm::ValueMap<float>& trackT0,
                const edm::ValueMap<float>& trackSigmaT0,    const edm::ValueMap<float>& trackTMTD,
                const edm::ValueMap<float>& trackSigmaTMTD,  const edm::ValueMap<float>& trackMom,
                const edm::ValueMap<float>& trackPathLength, const edm::ValueMap<float>& tofPIDT0,
                const edm::ValueMap<float>& tofPIDSigmaT0,   const std::vector<edm::ValueMap<reco::DeDxData> >& dEdxMaps,
                const edm::Handle<reco::GenParticleCollection>& genParticlesHandle,
                const edm::Handle<reco::RecoToSimCollection>& recToSimMapHandle,
                const reco::BeamSpot& bs, const MagneticField& mf)
{
  //
  // Do gen-reco matching
  //
  std::set<reco::GenParticleRef> genPars;
  if (genParticlesHandle.isValid()) {
    for (size_t iGen = 0; iGen < genParticlesHandle->size(); iGen++) {
      const auto& gen = reco::GenParticleRef(genParticlesHandle, iGen);
      if (gen->isLastCopy() && gen->status()==1 && gen->charge()!=0 && std::abs(gen->eta())<3.2) {
        genPars.insert(gen);
      }
    }
  }
  std::map<reco::GenParticleRef, bool> genMatch;
  std::map<reco::TrackRef, const TrackingParticle*> simMap;
  std::map<reco::GenParticleRef, std::vector<reco::TrackRef> > genTrkM;
  for (size_t i = 0; i < recoTracks->size(); i++) {
    const auto& trk = reco::TrackRef(recoTracks, i);
    // First try with sim info
    if (recToSimMapHandle.isValid()) {
      const auto& matchedSim = recToSimMapHandle->find(edm::RefToBase<reco::Track>(trk));
      if (matchedSim != recToSimMapHandle->end()) {
        const auto& sim = matchedSim->val.front().first.get();
        const auto& genV = sim->genParticles();
        simMap[trk] = sim;
        if (!genV.empty()) {
          for (const auto& gen : genV) {
            if (gen->status()==1) {
              genTrkM[gen].push_back(trk);
              genMatch[gen] = true;
              break;
            }
          }
        }
      }
    }
    // Then try using gen info
    if (!simMap[trk] && !genPars.empty()) {
      std::pair<double, reco::GenParticleRef> minDR(0.03, {});
      for (const auto& gen : genPars) {
        const auto deltaEta = trk->eta() - gen->eta();
        const auto deltaPhi = TVector2::Phi_mpi_pi(trk->phi() - gen->phi());
        const auto deltaR = std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
        const auto deltaPt = std::abs(trk->pt()*std::abs(gen->charge()) - gen->pt())/gen->pt();
        // Check if reco and gen match
        if (deltaR < minDR.first && deltaPt < 0.15 && (gen->charge()>0)==(trk->charge()>0))
          minDR = std::make_pair(deltaR, gen);
      }
      const auto& gen = minDR.second;
      if (gen.isNonnull() && !genMatch[gen]) genTrkM[gen].push_back(trk);
    }
  }
  std::map<reco::TrackRef, reco::GenParticleRef> genMap;
  for (const auto& g : genTrkM) {
    const auto& gen = g.first;
    if (gen.isNull()) continue;
    auto trk = g.second[0];
    if (g.second.size()>1) {
      for (const auto& obj : g.second) {
        const auto minDPt  = (std::abs(trk->pt()*std::abs(gen->charge()) - gen->pt())/gen->pt());
        const auto deltaPt = (std::abs(obj->pt()*std::abs(gen->charge()) - gen->pt())/gen->pt());
        if (deltaPt < minDPt) trk = obj;
      }
    }
    genMap[trk] = gen;
  }
  std::vector<reco::TrackRef> objV;
  std::set<const TrackingParticle*> simPars;
  // Fill the track information
  for (size_t iTrk = 0; iTrk < recoTracks->size(); iTrk++) {
    auto obj = reco::TrackRef(recoTracks, iTrk);
    // Select tracks with pTError < 10% of pt
    if ((obj->ptError()/obj->pt()) >= 0.1) continue;
    // Select tracks passing high purity
    if (obj->quality(reco::TrackBase::TrackQuality::highPurity) == false) continue;
    // Select tracks with more than 11 hits
    if (obj->numberOfValidHits() <= 11) continue;
    // Select tracks fitted with normalized chi2 lower than 7
    if (obj->normalizedChi2() >= 7.) continue;    
    // Extract the reco information
    const float beta_PID   = trackBeta[obj];
    const float t0_PID     = trackT0[obj];
    const float t0Err_PID  = trackSigmaT0[obj];
    const float tMTD       = trackTMTD[obj];
    const float tMTDErr    = trackSigmaTMTD[obj];
    const float magp       = trackMom[obj];
    const float pathLength = trackPathLength[obj];
    const float t0_TOF     = tofPIDT0[obj];
    const float t0Err_TOF  = tofPIDSigmaT0[obj];
    const float pt         = obj->pt();
    const float ptErr      = obj->ptError();
    const float eta        = obj->eta();
    const float phi        = obj->phi();
    const short charge     = obj->charge();
    const float normChi2   = obj->normalizedChi2();
    const int   quality    = obj->qualityMask();
    const int   nHits      = obj->numberOfValidHits();
    //find associated vertices
    int vtxIdx = -1;
    if (vtxCol.empty()) break;
    else if (vtxCol.size()==1) vtxIdx = 0;
    else {
      //first try based on association weights, but keep track of closest in z and z-t as well
      const auto rsigmazsq = 1./obj->dzError()/obj->dzError();
      int vtxIdxmindz = -1, vtxIdxminchisq = -1;
      double mindz = 0.1, minchisq = std::numeric_limits<double>::max(), maxW = 0.5;
      for (size_t ivtx = 0; ivtx<vtxCol.size(); ++ivtx) {
        const auto& vtx = *vtxCol[ivtx];
        const auto w = vtx.trackWeight(obj);
        if (w>maxW) {
          maxW = w;
          vtxIdx = ivtx;
        }
        const auto dz = std::abs(obj->dz(vtx.position()));
        if (dz<mindz) {
          mindz = dz;
          vtxIdxmindz = ivtx;
        }
        if (vtx.tError()>0. && vtx.tError()<0.025) {
          for (const auto& m : MASS) {
            const auto invBeta = std::sqrt(m*m + obj->p()*obj->p())/obj->p();
            const auto t0 = tMTD - (pathLength*invBeta/c_cm_ns);
            const auto dt = std::abs(t0-vtx.t());
            const auto dtsig = dt/tMTDErr;
            const auto chisq = dz*dz*rsigmazsq + dtsig*dtsig;
            if (chisq<minchisq) {
              minchisq = chisq;
              vtxIdxminchisq = ivtx;
            }
          }
        }
      }
      //if no vertex found based on association weights, fall back to closest in z or z-t
      if (vtxIdx<0) {
        //if closest vertex in z does not have valid time information, just use it, 
        //otherwise use the closest vertex in z-t plane with timing info, with a fallback to the closest in z
        const bool isBadT = (vtxIdxmindz>=0 && !(vtxCol[vtxIdxmindz]->tError()>0. && vtxCol[vtxIdxmindz]->tError()<0.025));
        if (isBadT) vtxIdx = vtxIdxmindz;
        else if (vtxIdxminchisq>=0) vtxIdx = vtxIdxminchisq;
        else if (vtxIdxmindz>=0) vtxIdx = vtxIdxmindz;
      }
      if (vtxIdx<0) continue; 
    }
    const auto& priVtx = *vtxCol.at(vtxIdx);
    const float t0_PV  = priVtx.t();
    // Compute the distance to primary vertex
    const float vx = priVtx.x();
    const float vy = priVtx.y();
    const float vz = priVtx.z();
    const float vxError = priVtx.xError();
    const float vyError = priVtx.yError();
    const float vzError = priVtx.zError();
    const math::XYZPoint vtx(vx, vy, vz);
    const float dz = obj->dz(vtx);
    const float dzsigma = std::sqrt(obj->dzError()*obj->dzError() + vzError*vzError);
    const float dcaDz = dz/dzsigma;
    const float dxy = obj->dxy(vtx);
    const float dxysigma = std::sqrt(obj->dxyError()*obj->dxyError() + vxError*vyError);
    const float dcaDxy = dxy/dxysigma;
    // Compute beta and mass using primary vertex info
    const float dt_PV = tMTD - t0_PV;
    const float beta_PV = (tMTDErr>=0. ? (pathLength/dt_PV)*(1./c_cm_ns) : -99.);
    // Store reco information
    this->Reco_Track_vtxIdx.push_back(vtxIdx);
    this->Reco_Track_beta_PID.push_back(beta_PID);
    this->Reco_Track_t0_PID.push_back(t0_PID);
    this->Reco_Track_t0Err_PID.push_back(t0Err_PID);
    this->Reco_Track_beta_PV.push_back(beta_PV);
    this->Reco_Track_t0_TOF.push_back(t0_TOF);
    this->Reco_Track_t0Err_TOF.push_back(t0Err_TOF);
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
    for (size_t i=0; i<dEdxMaps.size(); i++) {
      this->Reco_Track_dEdx[i].push_back(dEdxMaps[i][obj].dEdx());
    }
    if (genParticlesHandle.isValid()) {
      objV.push_back(obj);
      if (simMap[obj]) {
        simPars.insert(simMap[obj]);
      }
    } 
  }
  // Fill the reco-gen information
  std::map<reco::GenParticleRef, std::pair<size_t, size_t> > genTrkIdx;
  std::map<const TrackingParticle*, std::pair<size_t, size_t> > simTrkIdx;
  for (size_t iTrk=0; iTrk<objV.size(); iTrk++) {
    const auto& obj = objV[iTrk];
    if (!genParticlesHandle.isValid()) break;
    // Extract the gen information
    int gen_idx(-1);
    float gen_chi2(-1.);
    if (genMap[obj].isNonnull()) {
      const auto& gen = genMap[obj];
      if (genTrkIdx.find(gen)==genTrkIdx.end()) {
        genTrkIdx[gen] = std::make_pair(iTrk, std::distance(genPars.begin(), std::find(genPars.begin(), genPars.end(), gen)));
      }
      gen_idx = genTrkIdx[gen].second;
      // Compute chi2
      auto cov = obj->covariance(); cov.Invert();
      gen_chi2 = track_associator::trackAssociationChi2(obj->parameters(), cov, Basic3DVector<double>(gen->momentum()), Basic3DVector<double>(gen->vertex()), gen->charge(), mf, bs);
    }
    else if (simMap[obj]) {
      const auto& gen = simMap[obj];
      if (simTrkIdx.find(gen)==simTrkIdx.end()) {
        simTrkIdx[gen] = std::make_pair(iTrk, std::distance(simPars.begin(), std::find(simPars.begin(), simPars.end(), gen)));
      }
      gen_idx = genPars.size() + simTrkIdx[gen].second;
      gen_chi2 = 0;
    }
    // Store gen information
    this->Reco_Track_genIdx.push_back(gen_idx);
    this->Reco_Track_genChi2.push_back(gen_chi2);
  }
  // Fill the gen information
  for (const auto& gen : genPars) {
    this->Gen_Track_recoIdx.push_back(genTrkIdx.find(gen)!=genTrkIdx.end() ? genTrkIdx[gen].first : -1);
    this->Gen_Track_pdgId.push_back(gen->pdgId());
    this->Gen_Track_pt.push_back(gen->pt());
    this->Gen_Track_eta.push_back(gen->eta());
    this->Gen_Track_phi.push_back(gen->phi());
    this->Gen_Track_charge.push_back(gen->charge());
    this->Gen_Track_mass.push_back(gen->mass());
    this->Gen_Track_momPdgId.push_back(findMother(gen)->pdgId());
    this->Gen_Track_isGen.push_back(true);
  }
  // Fill the sim information
  for (const auto& gen : simPars) {
    this->Gen_Track_recoIdx.push_back(simTrkIdx.find(gen)!=simTrkIdx.end() ? simTrkIdx[gen].first : -1);
    this->Gen_Track_pdgId.push_back(gen->pdgId());
    this->Gen_Track_pt.push_back(gen->pt());
    this->Gen_Track_eta.push_back(gen->eta());
    this->Gen_Track_phi.push_back(gen->phi());
    this->Gen_Track_charge.push_back(gen->charge());
    this->Gen_Track_mass.push_back(gen->mass());
    auto mom = gen->parentVertex();
    this->Gen_Track_momPdgId.push_back(mom.isNonnull() && mom.isAvailable() ? mom->sourceTracks()[0]->pdgId() : -999);
    this->Gen_Track_isGen.push_back(false);
  }
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Fill(const edm::Event& iEvent, const std::vector<reco::VertexRef>& vtxCol, const std::map<std::string, std::map<reco::VertexRef, size_t> >& vtxNTrk, const float& hiHF, const int& hiBin, const float& genT0)
{
  this->Event_nRun   = iEvent.id().run();
  this->Event_nLumi  = iEvent.luminosityBlock();
  this->Event_nBX    = std::max(iEvent.bunchCrossing(),0);
  this->Event_nOrbit = std::max(iEvent.orbitNumber(),0);
  this->Event_nEvent = iEvent.id().event();
  // Primary Vertex Information
  for (const auto& pv : vtxCol) {
    this->Reco_Vertex_x.push_back(pv->x());
    this->Reco_Vertex_y.push_back(pv->y());
    this->Reco_Vertex_z.push_back(pv->z());
    this->Reco_Vertex_t0.push_back(pv->t());
    this->Reco_Vertex_t0Err.push_back(pv->tError());
    for (auto& p : this->Reco_Vertex_nTrk) p.second.push_back(vtxNTrk.at(p.first).at(pv));
  }
  // Centrality Information
  this->Event_hiHF   = hiHF;
  this->Event_hiBin  = hiBin;
  // Generator information
  this->Event_genT0  = genT0;
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
  tree_->Branch("Event_hiHF",            &(this->Event_hiHF),   "Event_hiHF/F"   );
  tree_->Branch("Event_hiBin",           &(this->Event_hiBin),  "Event_hiBin/S"  );
  tree_->Branch("Event_genT0",           &(this->Event_genT0),  "Event_genT0/F"  );
  // Vertex Info
  tree_->Branch("Reco_Vertex_x",         &(this->Reco_Vertex_x)                  );
  tree_->Branch("Reco_Vertex_y",         &(this->Reco_Vertex_y)                  );
  tree_->Branch("Reco_Vertex_z",         &(this->Reco_Vertex_z)                  );
  tree_->Branch("Reco_Vertex_t0",        &(this->Reco_Vertex_t0)                 );
  tree_->Branch("Reco_Vertex_t0Err",     &(this->Reco_Vertex_t0Err)              );
  for (auto& p : this->Reco_Vertex_nTrk) tree_->Branch(("Reco_Vertex_nTrk_"+p.first).c_str(), &(p.second));
  // Track Info
  tree_->Branch("Reco_Track_vtxIdx",     &(this->Reco_Track_vtxIdx)              );
  tree_->Branch("Reco_Track_beta_PID",   &(this->Reco_Track_beta_PID)            );
  tree_->Branch("Reco_Track_t0_PID",     &(this->Reco_Track_t0_PID)              );
  tree_->Branch("Reco_Track_t0Err_PID",  &(this->Reco_Track_t0Err_PID)           );
  tree_->Branch("Reco_Track_beta_PV",    &(this->Reco_Track_beta_PV)             );
  tree_->Branch("Reco_Track_t0_TOF",     &(this->Reco_Track_t0_TOF)              );
  tree_->Branch("Reco_Track_t0Err_TOF",  &(this->Reco_Track_t0Err_TOF)           );
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
  for (size_t i=0; i<this->Reco_Track_dEdx.size(); i++) {
    tree_->Branch(Form("Reco_Track_dEdx%lu",i), &(this->Reco_Track_dEdx[i])      );
  }
  tree_->Branch("Reco_Track_genIdx",     &(this->Reco_Track_genIdx)              );
  tree_->Branch("Reco_Track_genChi2",    &(this->Reco_Track_genChi2)             );
  // Gen Particle Info
  tree_->Branch("Gen_Track_recoIdx",     &(this->Gen_Track_recoIdx)              );
  tree_->Branch("Gen_Track_pdgId",       &(this->Gen_Track_pdgId)                );
  tree_->Branch("Gen_Track_pt",          &(this->Gen_Track_pt)                   );
  tree_->Branch("Gen_Track_eta",         &(this->Gen_Track_eta)                  );
  tree_->Branch("Gen_Track_phi",         &(this->Gen_Track_phi)                  );
  tree_->Branch("Gen_Track_charge",      &(this->Gen_Track_charge)               );
  tree_->Branch("Gen_Track_mass",        &(this->Gen_Track_mass)                 );
  tree_->Branch("Gen_Track_momPdgId",    &(this->Gen_Track_momPdgId)             );
  tree_->Branch("Gen_Track_isGen",       &(this->Gen_Track_isGen)                );
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
  this->Event_hiHF   = -1.;
  this->Event_hiBin  = -1;
  this->Event_genT0  = -1.;
  // Vertex Info
  this->Reco_Vertex_x.clear();
  this->Reco_Vertex_y.clear();
  this->Reco_Vertex_z.clear();
  this->Reco_Vertex_t0.clear();
  this->Reco_Vertex_t0Err.clear();
  for (auto& p : this->Reco_Vertex_nTrk) p.second.clear();
  // Track Info
  this->Reco_Track_vtxIdx.clear();
  this->Reco_Track_beta_PID.clear();
  this->Reco_Track_t0_PID.clear();
  this->Reco_Track_t0Err_PID.clear();
  this->Reco_Track_beta_PV.clear();
  this->Reco_Track_t0_TOF.clear();
  this->Reco_Track_t0Err_TOF.clear();
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
  for (size_t i=0; i<this->Reco_Track_dEdx.size(); i++) {
    this->Reco_Track_dEdx[i].clear();
  }
  this->Reco_Track_genIdx.clear();
  this->Reco_Track_genChi2.clear();
  // Gen Particle Info
  this->Gen_Track_recoIdx.clear();
  this->Gen_Track_pdgId.clear();
  this->Gen_Track_pt.clear();
  this->Gen_Track_eta.clear();
  this->Gen_Track_phi.clear();
  this->Gen_Track_charge.clear();
  this->Gen_Track_mass.clear();
  this->Gen_Track_momPdgId.clear();
  this->Gen_Track_isGen.clear();
}

//--------------------------------------------------------------------------------------------------
void
TimeEvent::Init(const size_t& dEdxSize)
{
  // Track dEdx Info
  if (this->Reco_Track_dEdx.empty()) {
    this->Reco_Track_dEdx.resize(dEdxSize);
  }
  if (this->Reco_Vertex_nTrk.empty()) {
    this->Reco_Vertex_nTrk = {{"AEta3", {}}, {"AEta0p5", {}}, {"AEta0p8", {}}};
  }
}

DEFINE_FWK_MODULE(TimeAnalyzer);
