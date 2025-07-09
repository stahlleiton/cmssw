#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackDeDxHits.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "CondFormats/SiPixelObjects/interface/PixelIndices.h"

#include "DataFormats/SiStripCommon/interface/ConstantsForHardwareSystems.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"

// My classes
#include "../interface/TBunchCrossing.h"
#include "../interface/TVertex.h"
#include "../interface/TTrack.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <fstream>
using namespace std;
using namespace reco;

/*****************************************************************************/
class MicroDstProducer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>
{
 public:
  explicit MicroDstProducer(const edm::ParameterSet& pset);
  ~MicroDstProducer();

  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(edm::Run const&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, const edm::EventSetup&) {};
  virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);

 private:
  void processHits(TTrack & r, const reco::DeDxHitInfo & info, const std::vector<float> & mom);
  void processEventRelated  (const edm::Event& ev);
  bool processVerticesTracks(const edm::EventSetup& es, const edm::Event& ev);

  // Root
  TTree * tree;
  TBunchCrossing * bunchCrossing;

  //
  const double MeVPerElectron_;
  const int VCaltoElectronGain_, VCaltoElectronGain_L1_, VCaltoElectronOffset_, VCaltoElectronOffset_L1_;
  const int pixelSaturationThr_;
  
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<std::vector<edm::Ptr<pat::PackedCandidate> > > track2pcSrc_;
  const edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxHitInfoToken_;
  const edm::EDGetTokenT<edm::ValueMap<std::vector<float>>> dedxHitMomToken_;

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tkTopoToken_;

  SiPixelGainCalibrationOfflineService pixelCalib_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  edm::ESHandle<TrackerTopology> tkTopo_;

  //
  //string outFile;
};

/*****************************************************************************/
MicroDstProducer::MicroDstProducer(const edm::ParameterSet& iConfig)
  : MeVPerElectron_(iConfig.getParameter<double>("MeVPerElectron")),
    VCaltoElectronGain_(iConfig.getParameter<int>("VCaltoElectronGain")),
    VCaltoElectronGain_L1_(iConfig.getParameter<int>("VCaltoElectronGain_L1")),
    VCaltoElectronOffset_(iConfig.getParameter<int>("VCaltoElectronOffset")),
    VCaltoElectronOffset_L1_(iConfig.getParameter<int>("VCaltoElectronOffset_L1")),
    pixelSaturationThr_(iConfig.getParameter<int>("pixelSaturationThr")),
    tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackProducer"))),
    track2pcSrc_(consumes<std::vector<edm::Ptr<pat::PackedCandidate> > >(iConfig.getParameter<edm::InputTag>("trackProducer"))),
    dedxHitInfoToken_(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedxHitInfo"))),
    dedxHitMomToken_(consumes<edm::ValueMap<std::vector<float>>>(iConfig.getParameter<edm::InputTag>("dedxMomentum"))),
    tkGeomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
    tkTopoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
    pixelCalib_(iConfig, consumesCollector())
{
  usesResource(TFileService::kSharedResource);
}

/*****************************************************************************/
MicroDstProducer::~MicroDstProducer()
{
}

/*****************************************************************************/
void MicroDstProducer::beginJob()
{
  // Root
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("hadronTree","HadronTree");
  bunchCrossing = new TBunchCrossing();
  tree->Branch("bunchCrossing", "TBunchCrossing", &bunchCrossing, 16000, 2);
}

/*****************************************************************************/
void MicroDstProducer::endJob()
{
}

/*****************************************************************************/
void MicroDstProducer::beginRun(edm::Run const&, const edm::EventSetup& iSetup) {
  tkGeom_ = iSetup.getHandle(tkGeomToken_);
  tkTopo_ = iSetup.getHandle(tkTopoToken_);
}

/*****************************************************************************/
void MicroDstProducer::processHits(TTrack & r, const reco::DeDxHitInfo & info, const std::vector<float> & mom)
{
  for(size_t i = 0; i < info.size(); i++)
  {
    const auto& type = info.type(i);

    const DetId& detId = info.detId(i);

    // Strip
    if (const auto& stripCluster = info.stripCluster(i))
    { 
      TStripHit hit;

      hit.thickness = dynamic_cast<const StripGeomDetUnit*>(tkGeom_->idToDet(detId))->surface().bounds().thickness();

      hit.chip = stripCluster->barycenter() / sistrip::STRIPS_PER_APV;

      // Collect adc
      for (const auto& adc : stripCluster->amplitudes())
        hit.adc.emplace_back(adc);

      //
      hit.detId = detId;

      hit.forCalib = (type & (1 << reco::DeDxHitInfo::Calibration));
      hit.forEloss = (type & (1 << reco::DeDxHitInfo::Complete)) &&
                     (type & (1 << reco::DeDxHitInfo::Compatible));

      hit.meas = 0; // dummy, measured cluster size

      hit.x = info.pathlength(i);

      hit.r = 0; // dummy, global r coordinate of the hit
      hit.p = mom[i];

      r.stripHits.emplace_back(hit);
    }

    // Pixel
    if (const auto& pixelCluster = info.pixelCluster(i))
    {
      TPixelHit hit;

      hit.thickness = dynamic_cast<const PixelGeomDetUnit*>(tkGeom_->idToDet(detId))->surface().bounds().thickness(); 

      hit.chip = (int(pixelCluster->x() / ROCSizeInX) << 3)
                + int(pixelCluster->y() / ROCSizeInY);

      double delta = 0;
      hit.isSaturated = false;

      for(size_t j = 0; j < pixelCluster->pixelADC().size(); j++)
      {
        const auto & elec = pixelCluster->pixelADC()[j];
        delta += elec * MeVPerElectron_;

        if(hit.isSaturated) continue; 

        const auto& row = pixelCluster->minPixelRow()
                        + pixelCluster->pixelOffset()[2 * j];
        const auto& col = pixelCluster->minPixelCol()
                        + pixelCluster->pixelOffset()[2 * j + 1];

        // Go back to adc
        const auto& DBgain = pixelCalib_.getGain(detId, col, row);
        const auto& DBpedestal = pixelCalib_.getPedestal(detId, col, row);

        if (elec == std::numeric_limits<uint16_t>::max())
          hit.isSaturated = true;
        else if (DBgain > 0.) {
          double vcal;
          const auto& theLayer = (detId.subdetId() == 1) ? tkTopo_->pxbLayer(detId) : 0;
          if (theLayer == 1)
            vcal = (elec - VCaltoElectronOffset_L1_) / VCaltoElectronGain_L1_;
          else
            vcal = (elec - VCaltoElectronOffset_) / VCaltoElectronGain_;

          const auto adc = std::round(DBpedestal + vcal / DBgain);

          if (adc > pixelSaturationThr_)
             hit.isSaturated = true;
        }
      }

      hit.Delta = delta;

      //
      hit.detId = detId;

      hit.forCalib = (type & (1 << reco::DeDxHitInfo::Calibration));
      hit.forEloss = (type & (1 << reco::DeDxHitInfo::Complete)) &&
                     (type & (1 << reco::DeDxHitInfo::Compatible));

      hit.meas = pair<short int,short int>(0,0); // dummy, measured cluster size

      hit.x = info.pathlength(i);

      hit.r = 0; // dummy, global r coordinate of the hit
      hit.p = mom[i];

      hit.nChannels = pixelCluster->pixelADC().size();

      r.pixelHits.emplace_back(hit);
    }
  }
}

/*****************************************************************************/
bool MicroDstProducer::processVerticesTracks
  (const edm::EventSetup& es, const edm::Event& ev)
{
  const auto & tracks = ev.getHandle(tracksToken_);
  const auto & track2pc = ev.getHandle(track2pcSrc_);
  const auto & dedxHitInfo = ev.get(dedxHitInfoToken_);
  const auto & dedxHitMom = ev.get(dedxHitMomToken_);

  pixelCalib_.setESObjects(es);

  LogDebug("microDstProducer")
         << " [MicroDstProducer] rectracks = " << tracks->size();

  // Process recvertices, rectracks 
  if (!tracks->empty() && tracks->size() < 100) // not too many
  //if(tracks->size() < 200) // not too many
  {
    TVertex recVertex;

    for (size_t i = 0; i < tracks->size(); i++)
    {
      const auto & track = reco::TrackRef(tracks, i);
      const auto & cand = track2pc.isValid() ? track2pc->at(i) : edm::Ptr<pat::PackedCandidate>();
      if (cand.isNonnull() && !dedxHitInfo.contains(cand.id()))
        continue;
      const auto & dedxHits = cand.isNonnull() ? dedxHitInfo[cand] : dedxHitInfo[track];
      if (dedxHits.isNull())
        continue;
      const auto & dedxMom = dedxHitMom[dedxHits];

      // track
      TTrack r;

      // rec
      r.isHighPurity = track->quality(reco::TrackBase::highPurity); 

      r.charge = track->charge();

      r.eta    = track->eta();
      r.pt     = track->pt();
      r.phi    = track->phi();

      r.chi2   = track->chi2();
      r.ndf    = track->ndof();
  
      processHits(r, *dedxHits, dedxMom);

      // Store track
      recVertex.tracks.push_back(r);
    }

    // Store vertex 
    bunchCrossing->recVertices.push_back(recVertex);
  }

  if(bunchCrossing->recVertices.empty())
    return false;

  LogDebug("microDstProducer")
    << " [MicroDstProducer] only "
    << bunchCrossing->recVertices.size() << " vertex with "
    << bunchCrossing->recVertices.front().tracks.size() << " tracks";

  return true;
}

/*****************************************************************************/
void MicroDstProducer::processEventRelated(const edm::Event& ev)
{
  bunchCrossing->runNumber   = ev.run();
  bunchCrossing->lumiSection = ev.luminosityBlock();
  bunchCrossing->bxNumber    = ev.bunchCrossing();
}

/*****************************************************************************/
void MicroDstProducer::analyze
  (const edm::Event& ev, const edm::EventSetup& es)
{
  LogDebug("microDstProducer") << "[MicroDstProducer]";

  LogDebug("microDstProducer")
    << " [MicroDstProducer] bunchCrossing number = " << ev.bunchCrossing();

  // Analyze
  processEventRelated(ev);

  if(processVerticesTracks(es,ev))
  {
    // Fill tree
    tree->Fill();
  }

  bunchCrossing->Clear();
}

DEFINE_FWK_MODULE(MicroDstProducer);

