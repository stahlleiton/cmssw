#ifndef HeavyIonsAnalysis_ElectroWeakAnalysis_HiMETAnalizer_h
#define HeavyIonsAnalysis_ElectroWeakAnalysis_HiMETAnalizer_h

//Headers for the core items
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

//Headers for the data items
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//Headers for services and tools
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
#include <TVector2.h>
#include <TMatrixD.h>

//system include files
#include <string>
#include <map>


typedef std::map< std::string , bool >  StringBoolMap;


class HiMETEvent
{
 public:
  // ===== Class Methods =====
  HiMETEvent    ( void );
  ~HiMETEvent   ( void );

  void  Clear       ( void );
  void  SetTree     ( TTree* tree ) {  this->tree = tree; }
  void  SetBranches ( std::string , const StringBoolMap&);

  void  Fill        ( const edm::Event&    );
  void  Fill        ( const reco::PFMET&   );
  void  Fill        ( const reco::CaloMET& );
  void  Fill        ( const reco::GenMET&  );
  void  Fill        ( const pat::MET& , std::string );
  void  Fill        ( const std::map< std::string, Bool_t >& );
  void  FillReco    ( const reco::MET&     );
  
 private:
  TTree*         tree;

  UInt_t         Event_nRun;
  UShort_t       Event_nLumi;
  UInt_t         Event_nBX;
  ULong64_t      Event_nEvent;

  TVector2       Reco_P2;
  TMatrixD       Reco_SigMatrix;
  Float_t        Reco_SumEt;
  Float_t        Reco_Significance;
  Float_t        Reco_mEtSig;

  TVector2       PF_P2;
  Float_t        PF_MuonEt;
  Float_t        PF_MuonEtFraction;
  Float_t        PF_ChargedEMEt;
  Float_t        PF_ChargedEMEtFraction;
  Float_t        PF_NeutralEMEt;
  Float_t        PF_NeutralEMEtFraction;
  Float_t        PF_ChargedHadEt;
  Float_t        PF_ChargedHadEtFraction;
  Float_t        PF_NeutralHadEt;
  Float_t        PF_NeutralHadEtFraction;
  Float_t        PF_HFEMEt;
  Float_t        PF_HFEMEtFraction;
  Float_t        PF_HFHadronEt;
  Float_t        PF_HFHadronEtFraction;

  TVector2       Calo_P2;
  Float_t        Calo_METInmHF;
  Float_t        Calo_METPhiInmHF;
  Float_t        Calo_SETInmHF;
  Float_t        Calo_METInpHF;
  Float_t        Calo_METPhiInpHF;
  Float_t        Calo_SETInpHF;
  Float_t        Calo_MaxEtInEmTowers;
  Float_t        Calo_MaxEtInHadTowers;
  Float_t        Calo_EtFractionHadronic;
  Float_t        Calo_EMEtFraction;
  Float_t        Calo_HadEtInHB;
  Float_t        Calo_HadEtInHE;
  Float_t        Calo_HadEtInHF;
  Float_t        Calo_HadEtInHO;
  Float_t        Calo_EMEtInEB;
  Float_t        Calo_EMEtInEE;
  Float_t        Calo_EMEtInHF;
  Float_t        Calo_MetSignificance;

  TVector2       Gen_P2;
  Float_t        Gen_InvisibleEt;
  Float_t        Gen_InvisibleEtFraction;
  Float_t        Gen_MuonEt;
  Float_t        Gen_MuonEtFraction;
  Float_t        Gen_ChargedEMEt;
  Float_t        Gen_ChargedEMEtFraction;
  Float_t        Gen_NeutralEMEt;
  Float_t        Gen_NeutralEMEtFraction;
  Float_t        Gen_ChargedHadEt;
  Float_t        Gen_ChargedHadEtFraction;
  Float_t        Gen_NeutralHadEt;
  Float_t        Gen_NeutralHadEtFraction;

  std::map< std::string , std::map< std::string , TVector2 > >  shiftedP2;
  std::map< std::string , std::map< std::string , Float_t > >   shiftedSumEt;
  std::map< std::string , Bool_t >                              Filter;
};


class HiMETAnalyzer : public edm::EDAnalyzer 
{
 public:
  explicit HiMETAnalyzer(const edm::ParameterSet&);
  virtual ~HiMETAnalyzer();

 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // TFileService
  edm::Service<TFileService> fs_;

  // Input Info
  const edm::EDGetTokenT< edm::View<pat::MET> >        _patMETToken;
  const edm::EDGetTokenT< edm::View<reco::PFMET> >     _pfMETToken;
  const edm::EDGetTokenT< edm::View<reco::CaloMET> >   _caloMETToken;
  const edm::EDGetTokenT< edm::View<reco::GenMET> >    _genMETToken;
  const edm::EDGetTokenT< edm::TriggerResults >        _filterResultToken;

  const std::vector< std::string >    _metFilters;
  const std::vector< std::string >    _corrNames;
  const bool                          _doAll;
  const std::string                   _eventFilter;

  // MET Containers
  std::map< std::string , TTree* >       metTree_;
  std::map< std::string , HiMETEvent >   metEvt_;

  // Flags
  bool firstEvent_;
};


// Global Variables
std::map< std::string, pat::MET::METCorrectionLevel> METCorrectionLevelMap = 
  {
    {"PF"             , pat::MET::Raw},
    {"Calo"           , pat::MET::RawCalo},
    {"Type1"          , pat::MET::Type1},
    {"Type01"         , pat::MET::Type01},
    {"TypeXY"         , pat::MET::TypeXY},
    {"Type1XY"        , pat::MET::Type1XY},
    {"Type01XY"       , pat::MET::Type01XY},
    {"Type1Smear"     , pat::MET::Type1Smear},
    {"Type01Smear"    , pat::MET::Type01Smear},
    {"Type1SmearXY"   , pat::MET::Type1SmearXY},
    {"Type01SmearXY"  , pat::MET::Type01SmearXY}
  };

std::map< std::string , pat::MET::METUncertainty>  METUncertaintyMap = 
  {
    {"NoShift"         , pat::MET::NoShift},
    {"JetResUp"        , pat::MET::JetResUp},
    {"JetResDown"      , pat::MET::JetResDown},
    {"JetEnUp"         , pat::MET::JetEnUp},
    {"JetEnDown"       , pat::MET::JetEnDown},
    {"MuonEnUp"        , pat::MET::MuonEnUp},
    {"MuonEnDown"      , pat::MET::MuonEnDown},
    {"ElectronEnUp"    , pat::MET::ElectronEnUp},
    {"ElectronEnDown"  , pat::MET::ElectronEnDown},
    {"TauEnUp"         , pat::MET::TauEnUp},
    {"TauEnDown"       , pat::MET::TauEnDown},
    {"UnclusEnUp"      , pat::MET::UnclusteredEnUp},
    {"UnclusEnDown"    , pat::MET::UnclusteredEnDown},
    {"PhotonEnUp"      , pat::MET::PhotonEnUp},
    {"PhotonEnDown"    , pat::MET::PhotonEnDown}
  };


// Global Functions  
template <class T>
static inline
void getCollection(const edm::Event & event, const edm::EDGetTokenT<T> token, edm::Handle<T> & handle) 
{
  event.getByToken(token, handle);
  if (!handle.isValid()) { handle.clear(); }
}



#endif
