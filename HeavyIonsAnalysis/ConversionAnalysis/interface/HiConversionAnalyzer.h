#ifndef HeavyIonsAnalysis_ConversionAnalysis_HiConversionAnalizer_h
#define HeavyIonsAnalysis_ConversionAnalysis_HiConversionAnalizer_h

//Headers for the core items
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Headers for the data items
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//Headers for services and tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//system include files
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <type_traits>
#include <vector>
#include <string>
#include <map>


typedef std::map< std::string, std::vector< Char_t > >   IndexMap;
typedef std::vector< std::string >                       StringVector;
typedef std::map< std::string , bool >                   StringBoolMap;
typedef edm::ESHandle<TransientTrackBuilder>             ESTransientTrackBuilder;


class HiConversionEvent
{
 public:
  // ===== Class Methods =====
  HiConversionEvent    ( void );
  ~HiConversionEvent   ( void );

  void   Clear   ( void );


  void   SetVertexCollection    ( const reco::VertexCollection& p ) { vtxCol_  = p;   }
  void   SetPrimaryVertex       ( const reco::Vertex&           p ) { privtx_  = p;   }
  void   SetTree                (       TTree*                  p ) { tree_    = p;   }
  void   SetBranches            ( std::string , const StringBoolMap&);
  void   IniTrackBuilder        ( const edm::EventSetup& );
  void   IniArrays              ( void );

  void   Fill                   ( const reco::ConversionCollection&,  const IndexMap& , const reco::PFCandidateCollection&, const reco::MuonCollection&  );
  void   Fill                   ( const reco::GenParticleRefVector&,  const IndexMap&        );
  void   Fill                   ( const reco::PFCandidateCollection&, const IndexMap&        );
  void   Fill                   ( const edm::Event&,                  const edm::EventSetup& );

 private:


  bool   CheckTkVtxCompatible   ( const reco::Conversion& );
  bool   FoundCompInnerHits     ( const reco::HitPattern& ,           const reco::HitPattern& );
  bool   ConversionEqualByTrack ( const reco::Conversion& ,           const reco::Conversion& );

  reco::GenParticleRef  findMotherRef ( const reco::GenParticle& );

  ESTransientTrackBuilder         _theTTBuilder;
  const std::vector < double >    _phoMasses = { 0.0 , 0.0 };
  const std::vector < double >    _eleMasses = { 0.000511 , 0.000511 , 0.000511 , 0.000511 };

  TTree*                          tree_;
  reco::Vertex                    privtx_;
  reco::VertexCollection          vtxCol_;

  // Event Info
  UInt_t                          Event_nRun;
  UShort_t                        Event_nLumi;
  UInt_t                          Event_nBX;
  ULong64_t                       Event_nOrbit;
  ULong64_t                       Event_nEvent;
  // Primary Vertex Information
  TVector3                        Event_PriVtx_Position;
  TVector3                        Event_PriVtx_Error;
  UChar_t                         Event_nPV;

  // Reco Conversion Kinematic
  UShort_t                   Reco_Conversion_N;
  TClonesArray*              Reco_Conversion_P4;
  TClonesArray*              Reco_Conversion_Vertex;
  std::vector < Char_t    >  Reco_Conversion_Charge;
  // Reco Convertion Matched Index
  std::vector < Short_t   >  Reco_Conversion_Gen_Index;
  std::vector < Short_t   >  Reco_Conversion_PF_Index;
  // Reco Conversion ID Flags
  std::vector < Bool_t    >  Reco_Conversion_isDuplicate;
  std::vector < Bool_t    >  Reco_Conversion_isPionLeg;
  std::vector < Bool_t    >  Reco_Conversion_isPriVtxCompatible;
  std::vector < Bool_t    >  Reco_Conversion_hasCompatibleInnerHits;
  std::vector < UInt_t    >  Reco_Conversion_Quality;
  // Reco Conversion ID vars
  std::vector < UShort_t  >  Reco_Conversion_NSharedHits;
  std::vector < UChar_t   >  Reco_Conversion_Algo;
  std::vector < Float_t   >  Reco_Conversion_MVA;
  std::vector < Float_t   >  Reco_Conversion_DCA;
  std::vector < Float_t   >  Reco_Conversion_VertexProb;
  std::vector < Float_t   >  Reco_Conversion_dXY;
  std::vector < Float_t   >  Reco_Conversion_dZ;
  std::vector < Float_t   >  Reco_Conversion_lXY;
  std::vector < Float_t   >  Reco_Conversion_lZ;
  std::vector < Float_t   >  Reco_Conversion_dPhiTracksAtVtx;
  std::vector < Float_t   >  Reco_Conversion_pairCotThetaSep;
  std::vector < Float_t   >  Reco_Conversion_EoverP;
  std::vector < UChar_t   >  Reco_Conversion_NTracks;
  // Reco Conversion Track 1
  TClonesArray*              Reco_Conversion_Track1_P3;
  std::vector < Float_t   >  Reco_Conversion_Track1_PtError;
  std::vector < Float_t   >  Reco_Conversion_Track1_NHitsBeforeVtx;
  std::vector < Float_t   >  Reco_Conversion_Track1_isHighPurity;
  std::vector < Float_t   >  Reco_Conversion_Track1_NumOfValHits;
  std::vector < Float_t   >  Reco_Conversion_Track1_NumOfLostHits;
  std::vector < Float_t   >  Reco_Conversion_Track1_NumOfValPixHits;
  std::vector < Float_t   >  Reco_Conversion_Track1_TrkLayersWithMea;
  std::vector < Float_t   >  Reco_Conversion_Track1_PixLayersWithMea;
  std::vector < Float_t   >  Reco_Conversion_Track1_dXY;
  std::vector < Float_t   >  Reco_Conversion_Track1_dXYErr;
  std::vector < Float_t   >  Reco_Conversion_Track1_dZ;
  std::vector < Float_t   >  Reco_Conversion_Track1_dZErr;
  std::vector < Float_t   >  Reco_Conversion_Track1_ValFrac;
  std::vector < Float_t   >  Reco_Conversion_Track1_NormChi2;
  // Reco Conversion Track 2
  TClonesArray*              Reco_Conversion_Track2_P3;
  std::vector < Float_t   >  Reco_Conversion_Track2_PtError;
  std::vector < Float_t   >  Reco_Conversion_Track2_NHitsBeforeVtx;
  std::vector < Float_t   >  Reco_Conversion_Track2_isHighPurity;
  std::vector < Float_t   >  Reco_Conversion_Track2_NumOfValHits;
  std::vector < Float_t   >  Reco_Conversion_Track2_NumOfLostHits;
  std::vector < Float_t   >  Reco_Conversion_Track2_NumOfValPixHits;
  std::vector < Float_t   >  Reco_Conversion_Track2_TrkLayersWithMea;
  std::vector < Float_t   >  Reco_Conversion_Track2_PixLayersWithMea;
  std::vector < Float_t   >  Reco_Conversion_Track2_dXY;
  std::vector < Float_t   >  Reco_Conversion_Track2_dXYErr;
  std::vector < Float_t   >  Reco_Conversion_Track2_dZ;
  std::vector < Float_t   >  Reco_Conversion_Track2_dZErr;
  std::vector < Float_t   >  Reco_Conversion_Track2_ValFrac;
  std::vector < Float_t   >  Reco_Conversion_Track2_NormChi2;
  // Reco DiConversion
  UInt_t                     Reco_DiConversion_N;
  TClonesArray*              Reco_DiConversion_P4;
  std::vector < Char_t    >  Reco_DiConversion_Charge;
  std::vector < UShort_t  >  Reco_DiConversion_Conversion1_Index;
  std::vector < UShort_t  >  Reco_DiConversion_Conversion2_Index;
  TClonesArray*              Reco_DiConversion_Vertex;
  std::vector < Float_t   >  Reco_DiConversion_VertexProb;
  std::vector < Float_t   >  Reco_DiConversion_MassError;
  // JUST FOR FUN, DIMUON+CONV
  TClonesArray*              Reco_DiMuonConv_P4;
  std::vector < UShort_t  >  Reco_DiMuonConv_Conversion_Index;
  std::vector < UShort_t  >  Reco_DiMuonConv_DiMuon_Index;
  std::vector < Float_t   >  Reco_Chi_Mass;
  std::vector < UChar_t   >  Reco_Chi_Type;

  // PF Photon
  UShort_t                   PF_Photon_N;
  TClonesArray*              PF_Photon_P4;
  std::vector < Short_t   >  PF_Photon_Gen_Index;
  std::vector < Short_t   >  PF_Photon_Reco_Index;
  // PF DiPhoton
  UInt_t                     PF_DiPhoton_N;
  TClonesArray*              PF_DiPhoton_P4;
  std::vector < UShort_t  >  PF_DiPhoton_Photon1_Index;
  std::vector < UShort_t  >  PF_DiPhoton_Photon2_Index;
 
  // Gen Photon
  UChar_t                    Gen_Photon_N;
  TClonesArray*              Gen_Photon_P4;
  std::vector < UShort_t  >  Gen_Photon_Particle_Index;
  std::vector < Short_t   >  Gen_Photon_Reco_Index;
  std::vector < Short_t   >  Gen_Photon_PF_Index;
  // Gen DiPhoton
  UChar_t                    Gen_DiPhoton_N;
  TClonesArray*              Gen_DiPhoton_P4;
  std::vector < Int_t     >  Gen_DiPhoton_PdgId;
  std::vector < UShort_t  >  Gen_DiPhoton_Photon1_Index;
  std::vector < UShort_t  >  Gen_DiPhoton_Photon2_Index;
};


class HiConversionAnalyzer : public edm::EDAnalyzer 
{
 public:
  explicit HiConversionAnalyzer(const edm::ParameterSet&);
  virtual ~HiConversionAnalyzer();

 private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // TFileService
  edm::Service<TFileService> fs_;

  // Input Info

  const edm::EDGetTokenT< edm::View < reco::Muon        > >   _recoMuonsToken;
  const edm::EDGetTokenT< edm::View < reco::Conversion  > >   _recoConversionsToken;
  const edm::EDGetTokenT< reco::GenParticleCollection     >   _genParticlesToken;
  const edm::EDGetTokenT< edm::View < reco::PFCandidate > >   _pfCandidatesToken;
  const edm::EDGetTokenT< edm::View < reco::Vertex      > >   _primaryVertexToken;
  const edm::EDGetTokenT< reco::BeamSpot                  >   _beamSpotToken;

  const bool            _doAll;

  // Conversion Containers
  std::map< std::string , TTree* >             convTree_;
  std::map< std::string , HiConversionEvent >  convEvt_;

  // Flags
  bool firstEvent_;

  // Comparators
  GreaterByPt<reco::PFCandidate>  pfCandidatePTComparator_;
  GreaterByPt<reco::GenParticle>  genParticlePTComparator_;
  static bool recoConversionGreaterByChi2_ (const reco::Conversion& c1, const reco::Conversion& c2) 
  {
    return TMath::Prob(c1.conversionVertex().chi2(),c1.conversionVertex().ndof()) < TMath::Prob(c2.conversionVertex().chi2(),c2.conversionVertex().ndof());
  };
  GreaterByPt<reco::Muon>         recoMuonPTComparator_;

};

// Global Functions  
template <class T>
static inline
void 
getCollection(const edm::Event & event, const edm::EDGetTokenT<T> token, edm::Handle<T> & handle) 
{
  event.getByToken(token, handle);
  if (!handle.isValid()) { handle.clear(); }
}

// Comparators
bool lowerByFirstElement(std::pair<double,short> a, std::pair<double,short> b) { return a.first < b.first; }

// define operator== for conversions, those with at least one track in common
namespace reco {
  bool operator==(const reco::Conversion& c1, const reco::Conversion& c2) {
    return c1.tracks()[0] == c2.tracks()[0] ||
           c1.tracks()[1] == c2.tracks()[1] ||
	   c1.tracks()[1] == c2.tracks()[0] ||
           c1.tracks()[0] == c2.tracks()[1] ;
  }
}

template<class T>  reco::LeafCandidate getRecoCandidate(const T& a) { return a; }
template<>         reco::LeafCandidate getRecoCandidate(const reco::Conversion& a) {
  reco::Candidate::LorentzVector p4;
  if (a.conversionVertex().isValid()) { p4.SetPxPyPzE( a.refittedPair4Momentum().Px() , a.refittedPair4Momentum().Py() , a.refittedPair4Momentum().Pz() , a.refittedPair4Momentum().E() ); }
  else { p4.SetPxPyPzE( a.pairMomentum().X() , a.pairMomentum().Y() , a.pairMomentum().Z() , a.pairMomentum().Rho() ); }
  reco::LeafCandidate o; o.setP4(p4); o.setCharge(0.0); return o;
}

template < class inT , class mT >
static inline
std::vector< std::vector< Char_t > >
doMatching(const std::vector<inT>& inC, const std::vector<mT>& mC, double maxDeltaR, double maxDPtRel)
{
  std::vector< Char_t > inV;
  std::map< Char_t, Char_t > mMap;
  for (ushort icand = 0; icand < inC.size(); icand++) {
    const reco::Candidate& inCand = getRecoCandidate(inC.at(icand));
    std::vector< std::pair< float , Char_t > > indexPair;
    for (ushort icand2 = 0; icand2 < mC.size(); icand2++) {
      const reco::Candidate& mCand = getRecoCandidate(mC.at(icand2));
      double deltaR = reco::deltaR( inCand.eta(), inCand.phi(), mCand.eta(), mCand.phi());
      double dPtRel = abs( inCand.pt() - mCand.pt() ) / (mCand.pt()+1E-9);
      if ( (deltaR < maxDeltaR) && (dPtRel < maxDPtRel) && (inCand.charge() == mCand.charge()) ) { indexPair.push_back( std::make_pair( deltaR , icand2 ) ); }
    }
    std::sort(indexPair.begin(), indexPair.end());
    inV.push_back( (indexPair.size() > 0) ? indexPair[0].second : -1 );
    if (indexPair.size() > 0) mMap[indexPair[0].second] = icand;
  }
  std::vector< Char_t > mV;
  for (ushort icand = 0; icand < mC.size(); icand++) {
    mV.push_back( (mMap.count(icand) > 0) ? mMap[icand] : -1 );
  }
  std::vector< std::vector< Char_t > > output;
  output.push_back( inV );
  output.push_back( mV );
  return output;
}


#endif
