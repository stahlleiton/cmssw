#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "DataFormats/Common/interface/Handle.h"
//#include "FWCore/Common/interface/Provenance.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "TTree.h"


class PixelClusterAnalyzer : public edm::one::EDAnalyzer<>
{
 public:
  PixelClusterAnalyzer(edm::ParameterSet const& conf);
  ~PixelClusterAnalyzer() override {};

  void analyze(const edm::Event& e, const edm::EventSetup& iSetup) override;

 private:
  TTree* t_;

  std::map<std::string, int> eventInfo;
  std::map<std::string, std::vector<float> > pixelInfo;

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geom_esToken_;
  const edm::ESGetToken<PixelClusterParameterEstimator, TkPixelCPERecord> pixelCPE_esToken_;
  const edm::EDGetTokenT<SiPixelClusterCollectionNew> tok_pixelClusters_;
};

PixelClusterAnalyzer::PixelClusterAnalyzer(edm::ParameterSet const& conf) :
  geom_esToken_(esConsumes()),
  pixelCPE_esToken_(esConsumes<PixelClusterParameterEstimator, TkPixelCPERecord>(edm::ESInputTag("", "PixelCPEGeneric"))),
  tok_pixelClusters_(consumes<SiPixelClusterCollectionNew>(edm::InputTag("siPixelClusters")))
{
  // open the tree file and initialize the tree
  edm::Service<TFileService> fs;
  t_ = fs->make<TTree>("pixelCluster", "");
  for (const auto& n : {"lumiBlock", "run", "bx", "orbit"})
    t_->Branch(Form("%s", n), &(eventInfo[n]), Form("%s/I", n));
  for (const auto& n : {"eta", "phi", "charge"})
    t_->Branch(Form("pixel_%s", n), &(pixelInfo[n]));
}

void PixelClusterAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{
  // fill event information
  eventInfo["lumiBlock"] = iEvent.luminosityBlock();
  eventInfo["run"] = iEvent.id().run();
  eventInfo["bx"] = iEvent.bunchCrossing();
  eventInfo["orbit"] = iEvent.orbitNumber();

  // fill pixel information
  const auto& trackerGeometry = iSetup.getData(geom_esToken_);
  const auto& pixelParam = iSetup.getData(pixelCPE_esToken_);
  const auto& pixelClusters = iEvent.get(tok_pixelClusters_);

  for (auto& p : pixelInfo)
    p.second.clear();

  int nPixelClusters(0);
  for (auto it = pixelClusters.begin(); it != pixelClusters.end(); it++) {
    const auto& detset = *it;
    const auto& id = it->detId();
    const auto& surface = trackerGeometry.idToDet(id)->surface();
    const auto& detUnit = *trackerGeometry.idToDetUnit(id);
    for (const auto& pixelCluster : detset) {
      const auto& pos = surface.toGlobal(pixelParam.localParametersV(pixelCluster, detUnit)[0].first);
      pixelInfo["eta"].emplace_back(pos.eta());
      pixelInfo["phi"].emplace_back(pos.phi());
      pixelInfo["charge"].emplace_back(pixelCluster.charge());
      nPixelClusters += 1;
    }
  }
  eventInfo["nPixelClusters"] = nPixelClusters;

  // fill the variables tree
  t_->Fill();
}

// declare this class as a framework plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PixelClusterAnalyzer);
