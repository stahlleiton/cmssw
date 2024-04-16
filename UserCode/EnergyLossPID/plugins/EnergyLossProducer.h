#ifndef EnergyLossProducer_H
#define EnergyLossProducer_H

#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"

//
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

//
#include "RecoTracker/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"

#include <fstream>

#include "../interface/gzstream.h"

namespace edm { class Event; class EventSetup; class Run; }

class Trajectory;
class DataHandler;
class HadronClusters;
class TTrack;

class EnergyLossProducer : public edm::one::EDProducer<>
{
public:
  explicit EnergyLossProducer(const edm::ParameterSet& ps);
  ~EnergyLossProducer();
  virtual void produce(edm::Event& ev, const edm::EventSetup& es);

private:
  void beginJob();
  void endJob();

  void readFitResults();
//  void beginRun(edm::Run const & run, edm::EventSetup const & es);
  std::pair<std::vector<float>, std::vector<float> > getProbabilities
      (const TTrack & track, const std::pair<double,double> & epsilon);

  //
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeom_;
  SiPixelGainCalibrationOfflineService theCalib;
  edm::ESGetToken<ClusterShapeHitFilter, CkfComponentsRecord> clusterShape_;

  //
  edm::EDGetTokenT<reco::TrackCollection>         trackProducer_;
  edm::EDGetTokenT<reco::TrackCollection>         trackProducer2_;
  edm::EDGetTokenT<std::vector<Trajectory> > trajectoryProducer_;
  edm::EDGetTokenT<SiPixelClusterShapeCache>  clusterShapeCache_;
  edm::EDGetTokenT<reco::DeDxDataValueMap>         dedxProducer_;

  //
  std::string tag;
  std::string outFile;

  DataHandler     * theDataHandler;
  HadronClusters  * theClusters;

  // charge, ieta, ipt, part, (amp,mean)
  std::map<std::pair<std::pair<int,int>,std::pair<int,int> >, double> amp, mean;
  // charge, ieta, ipt, (scale)
  std::map<std::pair<std::pair<int,int>,int>, double> scale;

  const int K = 3;

  const int    etaBins = 20;
  const double etaMin  = -1.;
  const double etaMax  =  1.;

  const int     ptBins = 40;
  const double  ptMin  = 0.;
  const double  ptMax  = 2.;

  ogzstream file;
};
#endif

