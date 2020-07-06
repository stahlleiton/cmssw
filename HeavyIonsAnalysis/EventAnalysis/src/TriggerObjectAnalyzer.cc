// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"

#include "TNtuple.h"
#include "TRegexp.h"

using namespace std;
using namespace edm;

//
// class declaration
//

class TriggerObjectAnalyzer : public edm::EDAnalyzer {
public:
  explicit TriggerObjectAnalyzer(const edm::ParameterSet&);
  ~TriggerObjectAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------

  const std::string processName_;
  const std::vector<std::string> triggerNames_;
  const edm::InputTag triggerResultsTag_;
  const edm::InputTag triggerEventTag_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;

  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;

  HLTConfigProvider hltConfig_;

  vector<string> moduleLabels_;

  edm::Service<TFileService> fs_;
  vector<TTree*> nt_;
  int verbose_;

  std::vector<std::string> triggerNamesInMenu_;
  std::map<std::string, bool> triggerInMenu_;

  struct OBJ {
    vector<short> id;
    vector<float> pt;
    vector<float> eta;
    vector<float> phi;
    vector<float> mass;
  };

  vector<OBJ> trgInfo_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerObjectAnalyzer::TriggerObjectAnalyzer(const edm::ParameterSet& ps):
  processName_(ps.getParameter<std::string>("processName")),
  triggerNames_(ps.getParameter<std::vector<std::string> >("triggerNames")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent")),
  triggerResultsToken_(consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerEventToken_(consumes<trigger::TriggerEvent>(triggerEventTag_))
{
  //now do what ever initialization is needed
  nt_.resize(triggerNames_.size());
  trgInfo_.resize(triggerNames_.size());
  triggerNamesInMenu_ = triggerNames_;
  for(size_t itrig=0; itrig<triggerNames_.size(); itrig++){
    nt_[itrig] = fs_->make<TTree>(triggerNames_.at(itrig).c_str(),Form("trigger %lud",itrig));

    auto& trg = trgInfo_[itrig];
    nt_[itrig]->Branch("TriggerObjID",&trg.id);
    nt_[itrig]->Branch("pt",&trg.pt);
    nt_[itrig]->Branch("eta",&trg.eta);
    nt_[itrig]->Branch("phi",&trg.phi);
    nt_[itrig]->Branch("mass",&trg.mass);
  }


  verbose_ = 0;
}


TriggerObjectAnalyzer::~TriggerObjectAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(hltConfig_.size() > 0){

    using namespace edm;
    iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
    iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);

    for(size_t itrig=0; itrig<triggerNames_.size(); itrig++) {
      auto& trg = trgInfo_[itrig];
      size_t filterIndex = 999999;
      const auto& triggerName = triggerNamesInMenu_[itrig];
      if (triggerName.rfind("HLT_", 0)==0) {
        if (triggerInMenu_.find(triggerName)==triggerInMenu_.end()) continue;
	    const auto& triggerIndex_ = hltConfig_.triggerIndex(triggerName);
        if (!triggerResultsHandle_->accept(triggerIndex_)) continue;
        const auto& mIndex = hltConfig_.moduleLabels(triggerIndex_).size()-1;
        // start from last modulein the path, iterate back until finding the filter that was run last
        for (int j = mIndex; j >= 0; --j) {
          const auto& filter = hltConfig_.moduleLabels(triggerIndex_).at(j); //this is simple to put into a loop to get all triggers...
          const auto& filterIdx = triggerEventHandle_->filterIndex(InputTag(filter,"",processName_));
          if (filterIdx<triggerEventHandle_->sizeFilters()) { filterIndex = filterIdx; break; }
        }
        assert(filterIndex<999999);
      }
      else if (triggerName.rfind("hlt", 0)==0) {
        const auto& filterIdx = triggerEventHandle_->filterIndex(InputTag(triggerName,"",processName_));
        if (filterIdx<triggerEventHandle_->sizeFilters()) { filterIndex = filterIdx; }
      }
      if (filterIndex<999999) {
	    const auto& VIDS = triggerEventHandle_->filterIds(filterIndex);
	    const auto& KEYS = triggerEventHandle_->filterKeys(filterIndex);
	    const auto& nI = VIDS.size();
	    const auto& nK = KEYS.size();
	    assert(nI==nK);
	    const auto& TOC = triggerEventHandle_->getObjects();
	    for (size_t i=0; i<nI; i++) {
	      const auto& TO = TOC[KEYS[i]];
	      trg.id.push_back(TO.id());
	      trg.pt.push_back(TO.pt());
	      trg.eta.push_back(TO.eta());
	      trg.phi.push_back(TO.phi());
	      trg.mass.push_back(TO.mass());
	    }
	  }
	}
  }

  for(size_t itrig=0; itrig<nt_.size(); itrig++){
	nt_[itrig]->Fill();
    trgInfo_[itrig] = OBJ();
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
TriggerObjectAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggerObjectAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
TriggerObjectAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const auto& activeHLTPathsInThisEvent = hltConfig_.triggerNames();
      triggerInMenu_.clear();
      for(size_t itrig=0; itrig<triggerNames_.size(); itrig++){
        auto& triggerName = triggerNamesInMenu_[itrig];
	    for (const auto& iHLT : activeHLTPathsInThisEvent) {
          //matching with regexp filter name. More than 1 matching filter is allowed so trig versioning is transparent to analyzer
          if (TString(iHLT).Contains(TRegexp(TString(triggerNames_[itrig])))){
            triggerInMenu_[iHLT] = true;
            triggerName = iHLT;
          }
        }
        if (triggerName.rfind("HLT_",0)==0) {
	      if (triggerInMenu_.find(triggerName)==triggerInMenu_.end()) {
            cout << "<HLT Object Analyzer> Warning! Trigger " << triggerName << " not found in HLTMenu. Skipping..." << endl;
          }
        }
        else if (triggerName.rfind("hlt",0)!=0) {
          cout << "<HLT Object Analyzer> Warning! Trigger name " << triggerName << " is not valid. Skipping..." << endl;
        }
      }
      if(verbose_){
	    hltConfig_.dump("ProcessName");
	    hltConfig_.dump("GlobalTag");
	    hltConfig_.dump("TableName");
	    hltConfig_.dump("Streams");
	    hltConfig_.dump("Datasets");
	    hltConfig_.dump("PrescaleTable");
	    hltConfig_.dump("ProcessPSet");
      }
    }
  } else {
    cout << "HLTObjectAnalyzer::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }
}

// ------------ method called when ending the processing of a run  ------------
void
TriggerObjectAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TriggerObjectAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TriggerObjectAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerObjectAnalyzer);
