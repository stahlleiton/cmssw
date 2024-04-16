#ifndef _HadronClusters_h_
#define _HadronClusters_h_

#include <fstream>
#include <vector>

namespace edm { class EventSetup; class ParameterSet; }

class SiPixelRecHit;
class SiStripRecHit2D;
class TrackingRecHit;


#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/GeometryVector/interface/LocalVector.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "RecoLocalTracker/SiPixelClusterizer/plugins/PixelClusterizerBase.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

class TrackerGeometry;
class ClusterShapeHitFilter;
class Trajectory;

class TTrack;

class HadronClusters
{
 public:
   HadronClusters(const TrackerGeometry * theTracker,
                  SiPixelGainCalibrationOfflineService & theCalib,
                  const ClusterShapeHitFilter * theClusterShape);

   ~HadronClusters();

  void analyzeRecTrajectory
    (const SiPixelClusterShapeCache & clusterShapeCache,
     const Trajectory & trajectory, TTrack & r);

 private:
  void processRec(
    const SiPixelClusterShapeCache & clusterShapeCache,
    const SiPixelRecHit &   recHit,
    LocalVector ldir, TTrack & r, float, float);

  void processRec(
    const SiStripRecHit2D & recHit,
    LocalVector ldir, TTrack & r, float, float);

  const TrackerGeometry * theTracker;
  SiPixelGainCalibrationOfflineService theCalib;
  const ClusterShapeHitFilter * theClusterShape;
};

#endif
