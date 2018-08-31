#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>
#include <cmath>

#include "TTree.h"

//rounds to first few nonzero sig figs
float rndSF(float value, int nSignificantDigits) 
{
  if(value==0) return 0; 

  float dSign = (value > 0.0) ? 1 : -1; 
  value *= dSign; 

  int nOffset = static_cast<int>(log10(value)); 
  if(nOffset>=0) ++nOffset;  

  float dScale = pow(10.0,nSignificantDigits-nOffset); 

  return dSign * static_cast<float>(round(value*dScale))/dScale;    
}

//keeps first n digits after the decimal place
inline float rndDP(float value, int nPlaces)
{
  return float(int(value*pow(10,nPlaces))/pow(10,nPlaces));
}

class TreePFCandEventData
{
  public:
    void SetTree(TTree * t) { tree_ = t; }
    void SetBranches(bool doJets, bool doMC);
    void Clear();

    Int_t nPFpart_;
    std::vector<Int_t> pfId_;
    std::vector<Float_t> pfPt_;
    std::vector<Float_t> pfEnergy_;
    std::vector<Float_t> pfEta_;
    std::vector<Float_t> pfPhi_;
    std::vector<Float_t> pfM_;

    Int_t nGENpart_;
    std::vector<Int_t> genPDGId_;
    std::vector<Float_t> genPt_;
    std::vector<Float_t> genEta_;
    std::vector<Float_t> genPhi_;

    Int_t njets_;
    std::vector<Float_t> jetEnergy_;
    std::vector<Float_t> jetPt_;
    std::vector<Float_t> jetEta_;
    std::vector<Float_t> jetPhi_;

  private:
    TTree* tree_;
};

class HiPFCandAnalyzer : public edm::EDAnalyzer {
  public:
    explicit HiPFCandAnalyzer(const edm::ParameterSet&);
    ~HiPFCandAnalyzer();

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::Service<TFileService> fs;

    // Event Info
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatePF_;
    edm::EDGetTokenT<reco::GenParticleCollection> genLabel_;
    edm::EDGetTokenT<pat::JetCollection> jetLabel_;

    TreePFCandEventData pfEvt_;
    TTree *pfTree_;

    // cuts
    Double_t pfPtMin_;
    Double_t pfAbsEtaMax_;
    Double_t jetPtMin_;
    Double_t genPtMin_;

    bool doJets_;
    bool doMC_;
    bool skipCharged_;
};
