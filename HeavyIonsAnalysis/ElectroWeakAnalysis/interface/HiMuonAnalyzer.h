#ifndef HeavyIonsAnalysis_ElectroWeakAnalysis_HiMuonAnalizer_h
#define HeavyIonsAnalysis_ElectroWeakAnalysis_HiMuonAnalizer_h

//Headers for the core items
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Headers for the data items
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

//Headers for services and tools
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//system include files
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <utility>
#include <vector>
#include <string>
#include <map>


typedef std::map< std::string, std::vector< Char_t > >   IndexMap;
typedef std::vector< std::string >                       StringVector;
typedef std::map< std::string , bool >                   StringBoolMap;
typedef edm::Handle< edm::TriggerResults >               TriggerResultHandle;
typedef edm::ESHandle<TransientTrackBuilder>             ESTransientTrackBuilder;

// Histograms
TH1F*                 hStats_;

class HiMuonEvent
{
 public:
  // ===== Class Methods =====
  HiMuonEvent    ( void );
  ~HiMuonEvent   ( void );

  void   Clear   ( void );

  void   SetPrimaryVertex       ( const reco::Vertex&        p ) { privtx_               = p;   }
  void   SetHLTTriggerNames     ( const StringVector&        p ) { triggerNames_         = p;   }
  void   SetHLTFilterNames      ( const StringVector&        p ) { filterNames_          = p;   }
  void   SetHLTTriggerResults   ( const TriggerResultHandle& p ) { triggerResultsHandle_ = p;   }
  void   SetHLTPrescaleProvider (       HLTPrescaleProvider* p ) { hltPrescaleProvider_  = p;   }
  void   SetTree                (       TTree*               p ) { tree_                 = p;   }
  void   SetBranches            ( std::string , const StringBoolMap&);
  void   IniTrackBuilder        ( const edm::EventSetup& );
  void   IniArrays              ( void );

  void   Fill                   ( const pat::MuonCollection&,         const IndexMap& );
  void   Fill                   ( const reco::GenParticleRefVector&,  const IndexMap& );
  void   Fill                   ( const reco::MuonCollection&,        const IndexMap&,        const short muonIndex = -1    );
  void   Fill                   ( const reco::PFCandidateCollection&, const IndexMap&,        const reco::VertexCollection& );
  void   Fill                   ( const edm::Event&,                  const edm::EventSetup&, const reco::VertexCollection& );

 private:

  bool   isMatched              ( const reco::Candidate&,             const reco::Candidate&,             double, double    );
  double pfIsolation            ( const reco::Candidate&,             const reco::PFCandidateCollection&, double, double    );
  short  findGenIndex           ( const reco::GenParticleRef&,        const reco::GenParticleRefVector&                     );

  ESTransientTrackBuilder         _theTTBuilder;
  const std::vector < double >    _muMasses = { 0.1056583715 , 0.1056583715 };

  TTree*                          tree_;
  reco::Vertex                    privtx_;
  HLTPrescaleProvider*            hltPrescaleProvider_;
  TriggerResultHandle             triggerResultsHandle_;
  StringVector                    triggerNames_;
  StringVector                    filterNames_;

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
  // Trigger Info
  std::vector < Bool_t         >  Event_Trigger_isFired;
  std::vector < Int_t          >  Event_Trigger_Prescale;

  // PAT Muon Info
  std::vector < std::vector < UChar_t > >  Pat_Muon_TriggerMatched;
  std::vector < Float_t        >  Pat_Muon_dB;
  std::vector < Float_t        >  Pat_Muon_dBErr;
  
  // Reco Muon Kinematic
  UChar_t                         Reco_Muon_N;
  TClonesArray*                   Reco_Muon_P4;
  std::vector < Char_t         >  Reco_Muon_Charge;
  // Reco Muon Matched Index
  std::vector < Char_t         >  Reco_Muon_Gen_Index;
  std::vector < Char_t         >  Reco_Muon_PF_Index;
  // Reco Muon ID Flags
  std::vector < Bool_t         >  Reco_Muon_isPFMuon;
  std::vector < Bool_t         >  Reco_Muon_isGlobalMuon;
  std::vector < Bool_t         >  Reco_Muon_isTrackerMuon;
  std::vector < Bool_t         >  Reco_Muon_isStandAloneMuon;
  std::vector < Bool_t         >  Reco_Muon_isLooseMuon;
  std::vector < Bool_t         >  Reco_Muon_isMediumMuon;
  std::vector < Bool_t         >  Reco_Muon_isHighPtMuon;
  std::vector < Bool_t         >  Reco_Muon_isSoftMuon;
  std::vector < Bool_t         >  Reco_Muon_isTightMuon;
  std::vector < Bool_t         >  Reco_Muon_isArbitrated;
  std::vector < Bool_t         >  Reco_Muon_TrackerMuonArbitrated;
  std::vector < Bool_t         >  Reco_Muon_GlobalMuonPromptTight;
  std::vector < Bool_t         >  Reco_Muon_TMLastStationLoose;
  std::vector < Bool_t         >  Reco_Muon_TMLastStationTight;
  std::vector < Bool_t         >  Reco_Muon_TM2DCompatibilityLoose;
  std::vector < Bool_t         >  Reco_Muon_TM2DCompatibilityTight;
  std::vector < Bool_t         >  Reco_Muon_TMOneStationLoose;
  std::vector < Bool_t         >  Reco_Muon_TMOneStationTight;
  std::vector < Bool_t         >  Reco_Muon_GMTkChiCompatibility;
  std::vector < Bool_t         >  Reco_Muon_GMStaChiCompatibility;
  std::vector < Bool_t         >  Reco_Muon_GMTkKinkTight;
  std::vector < Bool_t         >  Reco_Muon_TMLastStationAngLoose;
  std::vector < Bool_t         >  Reco_Muon_TMLastStationAngTight;
  std::vector < Bool_t         >  Reco_Muon_TMOneStationAngLoose;
  std::vector < Bool_t         >  Reco_Muon_TMOneStationAngTight;
  // Reco Muon ID vars
  std::vector < Short_t        >  Reco_Muon_NumberOfMatchedStations;
  std::vector < Short_t        >  Reco_Muon_NumberOfMatches;
  std::vector < Float_t        >  Reco_Muon_SegmentCompatibility;
  std::vector < Float_t        >  Reco_Muon_Chi2LocalPosition;
  std::vector < Float_t        >  Reco_Muon_TrkKink;
  // Reco Muon Inner Track
  TClonesArray*                   Reco_Muon_InnerTrack_P3;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_PtError;
  std::vector < Bool_t         >  Reco_Muon_InnerTrack_isHighPurity;
  std::vector < Short_t        >  Reco_Muon_InnerTrack_NumOfValHits;
  std::vector < Short_t        >  Reco_Muon_InnerTrack_NumOfLostHits;
  std::vector < Short_t        >  Reco_Muon_InnerTrack_NumOfValPixHits;
  std::vector < Short_t        >  Reco_Muon_InnerTrack_TrkLayersWithMea;
  std::vector < Short_t        >  Reco_Muon_InnerTrack_PixLayersWithMea;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_dXY;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_dXYErr;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_dZ;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_dZErr;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_ValFrac;
  std::vector < Float_t        >  Reco_Muon_InnerTrack_NormChi2;
  // Reco Muon Global Track
  TClonesArray*                   Reco_Muon_GlobalTrack_P3;
  std::vector < Float_t        >  Reco_Muon_GlobalTrack_PtError;
  std::vector < Short_t        >  Reco_Muon_GlobalTrack_NumOfValMuonHits;
  std::vector < Float_t        >  Reco_Muon_GlobalTrack_NormChi2;
  // Reco Muon Best Track
  TClonesArray*                   Reco_Muon_BestTrack_P3;
  TClonesArray*                   Reco_Muon_BestTrack_Vertex;
  std::vector < Char_t         >  Reco_Muon_BestTrack_Type;
  std::vector < Float_t        >  Reco_Muon_BestTrack_PtError;
  std::vector < Float_t        >  Reco_Muon_BestTrack_dXY;
  std::vector < Float_t        >  Reco_Muon_BestTrack_dXYErr;
  std::vector < Float_t        >  Reco_Muon_BestTrack_dZ;
  std::vector < Float_t        >  Reco_Muon_BestTrack_dZErr;
  // Reco Muon Isolation
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_ChargedEM_SumPt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_NeutralEM_SumEt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_ChargedHad_SumPt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_NeutralHad_SumEt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_ChargedHadPU_SumPt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_IsoPUCorr;
  std::vector < Float_t        >  Reco_Muon_PFIsoR03_IsoNoPUCorr;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_ChargedEM_SumPt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_ChargedHad_SumPt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_NeutralHad_SumEt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_NeutralEM_SumEt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_ChargedHadPU_SumPt;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_IsoPUCorr;
  std::vector < Float_t        >  Reco_Muon_PFIsoR04_IsoNoPUCorr;
  std::vector < Float_t        >  Reco_Muon_IsoR03_SumPt;
  std::vector < Float_t        >  Reco_Muon_IsoR03_Iso;
  std::vector < Float_t        >  Reco_Muon_IsoR05_SumPt;
  std::vector < Float_t        >  Reco_Muon_IsoR05_Iso;
  // Reco DiMuon
  UShort_t                        Reco_DiMuon_N;
  TClonesArray*                   Reco_DiMuon_P4;
  std::vector < Char_t         >  Reco_DiMuon_Charge;
  std::vector < UChar_t        >  Reco_DiMuon_Muon1_Index;
  std::vector < UChar_t        >  Reco_DiMuon_Muon2_Index;
  std::vector < Bool_t         >  Reco_DiMuon_isCowBoy;
  TClonesArray*                   Reco_DiMuon_Vertex;
  std::vector < Float_t        >  Reco_DiMuon_VertexProb;
  std::vector < Float_t        >  Reco_DiMuon_DCA;
  std::vector < Float_t        >  Reco_DiMuon_CTau;
  std::vector < Float_t        >  Reco_DiMuon_CTauErr;
  std::vector < Float_t        >  Reco_DiMuon_CosAlpha;
  std::vector < Float_t        >  Reco_DiMuon_MassError;

  // PF Candidate
  std::vector < Bool_t         >  PF_Candidate_isPU;
  std::vector < UChar_t        >  PF_Candidate_Id;
  std::vector < Float_t        >  PF_Candidate_Eta;
  std::vector < Float_t        >  PF_Candidate_Phi;
  std::vector < Float_t        >  PF_Candidate_Pt;
  // PF Muon
  UChar_t                         PF_Muon_N;
  TClonesArray*                   PF_Muon_P4;
  std::vector < Char_t         >  PF_Muon_Charge;
  std::vector < Char_t         >  PF_Muon_Gen_Index;
  std::vector < Char_t         >  PF_Muon_Reco_Index;
  // Reco Muon Isolation
  std::vector < Float_t        >  PF_Muon_PFIsoR03_ChargedEM_SumPt;
  std::vector < Float_t        >  PF_Muon_PFIsoR03_NeutralEM_SumEt;
  std::vector < Float_t        >  PF_Muon_PFIsoR03_ChargedHad_SumPt;
  std::vector < Float_t        >  PF_Muon_PFIsoR03_NeutralHad_SumEt;
  std::vector < Float_t        >  PF_Muon_PFIsoR03_ChargedHadPU_SumPt;
  std::vector < Float_t        >  PF_Muon_PFIsoR03_IsoPUCorr;
  std::vector < Float_t        >  PF_Muon_PFIsoR03_IsoNoPUCorr;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_ChargedEM_SumPt;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_ChargedHad_SumPt;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_NeutralHad_SumEt;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_NeutralEM_SumEt;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_ChargedHadPU_SumPt;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_IsoPUCorr;
  std::vector < Float_t        >  PF_Muon_PFIsoR04_IsoNoPUCorr;
  // PF DiMuon
  UShort_t                        PF_DiMuon_N;
  TClonesArray*                   PF_DiMuon_P4;
  std::vector < Char_t         >  PF_DiMuon_Charge;
  std::vector < UChar_t        >  PF_DiMuon_Muon1_Index;
  std::vector < UChar_t        >  PF_DiMuon_Muon2_Index;
  TClonesArray*                   PF_DiMuon_Vertex;
  std::vector < Float_t        >  PF_DiMuon_VertexProb;
  std::vector < Float_t        >  PF_DiMuon_DCA;
  std::vector < Float_t        >  PF_DiMuon_MassError;
  // PF MET
  TClonesArray*                   PF_MET_P2;
  // PF Muon-MET
  TClonesArray*                   PF_MuonMET_P4T;
  // Gen Particle
  TClonesArray*                   Gen_Particle_P4;
  std::vector < Int_t          >  Gen_Particle_PdgId;
  std::vector < UChar_t        >  Gen_Particle_Status;
  std::vector < std::vector < UShort_t > >  Gen_Particle_Mother_Index;
  std::vector < std::vector < UShort_t > >  Gen_Particle_Daughter_Index;
  // Gen Muon
  UChar_t                         Gen_Muon_N;
  TClonesArray*                   Gen_Muon_P4;
  std::vector < Char_t         >  Gen_Muon_Charge;
  std::vector < UShort_t       >  Gen_Muon_Particle_Index;
  std::vector < Char_t         >  Gen_Muon_Reco_Index;
  std::vector < Char_t         >  Gen_Muon_PF_Index;
};


class HiMuonAnalyzer : public edm::EDAnalyzer 
{
 public:
  explicit HiMuonAnalyzer(const edm::ParameterSet&);
  virtual ~HiMuonAnalyzer();

 private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // TFileService
  edm::Service<TFileService> fs_;

  // Input Info
  const edm::EDGetTokenT< edm::View < pat::Muon         > >   _patMuonsToken;
  const edm::EDGetTokenT< edm::View < reco::Muon        > >   _recoMuonsToken;
  const edm::EDGetTokenT< reco::GenParticleCollection     >   _genParticlesToken;
  const edm::EDGetTokenT< edm::View < reco::PFCandidate > >   _pfCandidatesToken;
  const edm::EDGetTokenT< edm::View < reco::Vertex      > >   _primaryVertexToken;
  const edm::EDGetTokenT< reco::BeamSpot                  >   _beamSpotToken;
  const edm::EDGetTokenT< edm::TriggerResults             >   _triggerResultsToken;

  const StringVector    _triggerPathNames;
  const StringVector    _triggerFilterNames;
  HLTPrescaleProvider   _hltPrescaleProvider;
  const bool            _doAll;

  // Muon Containers
  std::map< std::string , TTree* >        muonTree_;
  std::map< std::string , HiMuonEvent >   muonEvt_;

  // Flags
  bool firstEvent_;

  // PT Comparators
  GreaterByPt<pat::Muon>          patMuonPTComparator_;
  GreaterByPt<reco::Muon>         recoMuonPTComparator_;
  GreaterByPt<reco::PFCandidate>  pfCandidatePTComparator_;
  GreaterByPt<reco::GenParticle>  genParticlePTComparator_;
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

bool checkObject(const reco::Muon& inC        , const reco::PFCandidate& mC ) { return inC.isPFMuon(); }
bool checkObject(const reco::PFCandidate& inC , const reco::Muon& mC        ) { return mC.isPFMuon();  }
bool checkObject(const reco::Candidate& inC   , const reco::Candidate& mC   ) { return true;           }

template < class inT , class mT >
static inline
std::vector< std::vector< Char_t > >
doMatching(const std::vector<inT>& inC, const std::vector<mT>& mC, double maxDeltaR, double maxDPtRel)
{
  std::vector< Char_t > inV;
  std::map< Char_t, Char_t > mMap;
  for (ushort icand = 0; icand < inC.size(); icand++) {
    const reco::Candidate& inCand = inC.at(icand);
    std::vector< std::pair< float , Char_t > > indexPair;
    for (ushort icand2 = 0; icand2 < mC.size(); icand2++) {
      const reco::Candidate& mCand = mC.at(icand2);
      if ( !checkObject(inC.at(icand), mC.at(icand2)) ) continue;
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
