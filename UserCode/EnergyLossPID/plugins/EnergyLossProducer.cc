#include "EnergyLossProducer.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

//
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

//
#include "../interface/DataHandler.h"
#include "../interface/HadronClusters.h"
#include "../interface/TTrack.h"

using namespace std;

#undef DebugOutput

/*****************************************************************************/
EnergyLossProducer::EnergyLossProducer(const edm::ParameterSet& pset)
  : trackerGeom_( esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
    theCalib(pset, consumesCollector()),
    clusterShape_(esConsumes<ClusterShapeHitFilter, CkfComponentsRecord>())
{
  cerr << "\033[22;31m" << "Getting parameters.."
       << "\033[22;0m"  << endl;

  trackProducer_      = consumes<reco::TrackCollection>(
                        pset.getParameter<edm::InputTag>("trackProducer"));
  trackProducer2_     = consumes<reco::TrackCollection>(
                        pset.getParameter<edm::InputTag>("trackProducer2"));

  trajectoryProducer_ = consumes<vector<Trajectory> >(
                        pset.getParameter<edm::InputTag>("trajectoryProducer"));

  clusterShapeCache_  = consumes<SiPixelClusterShapeCache>(
                        pset.getParameter<edm::InputTag>("clusterShapeCache"));

  dedxProducer_       = consumes<reco::DeDxDataValueMap>(
                        pset.getParameter<edm::InputTag>("dedxProducer"));

  //
  tag           = pset.getParameter<string>("tag");
  outFile       = pset.getParameter<string>("outFile");

  cerr << "\033[22;31m" << "Setting up products.."
       << "\033[22;0m"  << endl;
  produces<reco::DeDxDataValueMap>("energyLossPixHits");
  produces<reco::DeDxDataValueMap>("energyLossStrHits");
  produces<reco::DeDxDataValueMap>("energyLossAllHits");

  cerr << "\033[22;31m" << "Starting DataHandler (" << tag << ").."
       << "\033[22;0m"  << endl;

  theDataHandler = new DataHandler(tag);
  theDataHandler->applyGain = true;
  theDataHandler->beginJob();

  cerr << "\033[22;31m" << "Starting HadronClusters.."
       << "\033[22;0m"  << endl;

  // read fit results 
  cerr << "\033[22;31m" << "Reading fit amplitudes.."
       << "\033[22;0m"  << endl;
  readFitResults();
}

/*****************************************************************************/
EnergyLossProducer::~EnergyLossProducer()
{
}

/*****************************************************************************/
void EnergyLossProducer::beginJob()
{
#ifdef DebugOutput
  file.open(outFile.c_str());
#endif
}

/*****************************************************************************/
void EnergyLossProducer::endJob()
{
#ifdef DebugOutput
  file.close();
#endif
}

/*****************************************************************************/
void EnergyLossProducer::readFitResults()
{
  char fileName[256];
  sprintf(fileName,"UserCode/EnergyLossPID/data/mostprob_%s.dat", tag.c_str());

  edm::FileInPath fileInPath(fileName);
  ifstream file(fileInPath.fullPath().c_str());

  cerr << " from " << fileName << endl;

  while(!file.eof())
  {
    int icha,ieta,ipt, d;
    float mue, p;

    file >> icha >> ieta >> ipt >> p;

    for(int k = 0; k < K ; k++)
    {
      pair<pair<int,int>,pair<int,int> > key(pair<int,int>(icha,ieta),
                                             pair<int,int>( ipt,   k));

      file >> mean[key] >> mue >> amp[key];
    }

    pair<pair<int,int>,int> key(pair<int,int>(icha,ieta), ipt);

    file >> scale[key] >> d >> d;
  }

  file.close();
}

/*****************************************************************************/
pair<vector<float>, vector<float> > EnergyLossProducer::getProbabilities
   (const TTrack & track, const pair<double,double> & epsilon)
{
  pair<vector<float>,vector<float> >
   res(vector<float>(3,0.), vector<float>(3,999.));

  if(track.eta > etaMin && track.eta < etaMax &&
     track.pt  >  ptMin && track.pt  <  ptMax)
  {
    //
    int icha = (track.charge > 0 ? pos : neg);
    int ieta = int((track.eta - etaMin)/(etaMax - etaMin) * etaBins);
    int ipt  = int((track.pt  -  ptMin)/( ptMax -  ptMin) *  ptBins);

    //
    double    logde = log(epsilon.first);
    double siglogde = sqrt(epsilon.second) / epsilon.first;

    //
    pair<pair<int,int>,int> key(pair<int,int>(icha,ieta), ipt);

    double sig = scale[key] * siglogde;

    for(int k = 0; k < K; k++)
    {
      pair<pair<int,int>,pair<int,int> > key(pair<int,int>(icha,ieta),
                                             pair<int,int>( ipt,   k));

      double q = (logde - mean[key]) / sig;

      if(fabs(q) < 10) // FIXME
      {
        res.first[k]  = amp[key] * exp(-0.5*q*q); // division by sig not needed
        res.second[k] = q;
      }
    }

    double s = 0.;
      for(int k = 0; k < K; k++) s += res.first[k];

    if(s > 0.)
      for(int k = 0; k < K; k++) res.first[k] /= s;
  }

  return res;
}

/*****************************************************************************/
void EnergyLossProducer::produce(edm::Event& ev, const edm::EventSetup& es)
{
  //
  const TrackerGeometry       * theTracker      = &es.getData(trackerGeom_);
  const ClusterShapeHitFilter * theClusterShape = &es.getData(clusterShape_);

  theCalib.setESObjects(es);
  theClusters = new HadronClusters(theTracker, theCalib, theClusterShape);

  edm::Handle<SiPixelClusterShapeCache>  clusterShapeCache;
  ev.getByToken(     clusterShapeCache_, clusterShapeCache);

  // Get track collection
  edm::Handle<reco::TrackCollection> trackHandle;
  ev.getByToken(trackProducer_,      trackHandle);
  const vector<reco::Track> & trackCollection =
                                   *(trackHandle.product());
  const auto& trackHandle2 = ev.getHandle(trackProducer2_);
  const auto& trackCollection2 = ev.get(trackProducer2_);

  // Not too many
  if(trackCollection.size() == 0 ||
     trackCollection.size() >= 100) return;

  // Get official dE/dx collection
  edm::Handle<reco::DeDxDataValueMap> dedxHandle;
  ev.getByToken(dedxProducer_,        dedxHandle);

  auto outputPix = std::make_unique<reco::DeDxDataValueMap>(); 
  auto outputStr = std::make_unique<reco::DeDxDataValueMap>(); 
  auto outputAll = std::make_unique<reco::DeDxDataValueMap>(); 

  reco::DeDxDataValueMap::Filler fillerPix(*outputPix);
  reco::DeDxDataValueMap::Filler fillerStr(*outputStr);
  reco::DeDxDataValueMap::Filler fillerAll(*outputAll);

  // Get trajectory collection
  edm::Handle<vector<Trajectory> >   trajeHandle;
  ev.getByToken(trajectoryProducer_, trajeHandle);
  const vector<Trajectory> & trajeCollection =
                           *(trajeHandle.product());

  //
  vector<reco::DeDxData> estimatePix;
  vector<reco::DeDxData> estimateStr;
  vector<reco::DeDxData> estimateAll;

  // Associate tracks
  std::map<size_t, std::pair<size_t, double>> trkMap, trkMap2;
  for (size_t i=0; i<trackCollection.size(); i++) {
    const auto& trk = trackCollection[i];
    for (size_t j=0; j<trackCollection2.size(); j++) {
      const auto& trk2 = trackCollection2[j];
      if (trk.charge() != trk2.charge() || abs(trk2.pt() - trk.pt())/(trk2.pt()<0.5 ? 1. : trk2.pt()) > 1.)
        continue;
      const auto dR = reco::deltaR(trk, trk2);
      const auto& it = trkMap.find(j);
      const auto minDR = (it!=trkMap.end() ? it->second.second : 0.5);
      if (dR < minDR)
        trkMap[j] = {i, dR};
      auto it2 = trkMap2.find(i);
      const auto minDR2 = (it2!=trkMap2.end() ? it2->second.second : 0.5);
      if (dR < minDR2)
        trkMap2[i] = {j, dR};
    }
  }
  for (const auto& t : trkMap)
    if (t.first != trkMap2[t.second.first].first)
      trkMap[t.first] = {-1, -1.};

  // Take all trajectories
  for (size_t iTrk=0; iTrk<trackCollection2.size(); iTrk++)
  {
    TTrack track;

    const auto& it = trkMap.find(iTrk);
    const int j = (it!=trkMap.end() ? it->second.first : -1);
    if (j<0) {
      const auto& trkRef = reco::TrackRef(trackHandle2, iTrk);
      const auto& dedx = (*dedxHandle)[trkRef];
      estimatePix.emplace_back(dedx);
      estimateStr.emplace_back(dedx);
      estimateAll.emplace_back(dedx);
      continue;
    }

    const auto& trajectory = &trajeCollection[j];
    track.charge = trackCollection[j].charge();
    track.eta    = trackCollection[j].eta();
    track.pt     = trackCollection[j].pt();

    theClusters->analyzeRecTrajectory(*clusterShapeCache, *trajectory,track);

    pair<double,double> epsilon = 
      theDataHandler->processTrack(track, estimatePix,
                                          estimateStr,
                                          estimateAll);

    //
    pair<vector<float>, vector<float> > res = getProbabilities(track, epsilon);

    vector<float> pro = res.first;
    vector<float> dev = res.second;

#ifdef DebugOutput
    double p = track.pt * cosh(track.eta);

    reco::TrackRef trackref = reco::TrackRef( trackHandle, j );

    const reco::DeDxData & dedxOfficial = dedxHandle->get(trackref.key()); 

    // DataFormats/TrackReco/src/DeDxData.cc always returns dEdxError() == -1
    file << " " << p
         << " " << estimateAll.back().dEdx()                 // epsilon
         << " " << sqrt(epsilon.second)                      // sigma(epsilon)
         << " " << estimateAll.back().numberOfMeasurements() // nhits
         << " " << pro[0] // prob pion
         << " " << pro[1] // prob kaon
         << " " << pro[2] // prob prot (sum = 1)
         << " " << dev[0] // deviation pion
         << " " << dev[1] // deviation kaon
         << " " << dev[2] // deviation prot (in units of sigma)
         << " " << dedxOfficial.dEdx()
         //
         << endl; 
#endif
  }

  fillerPix.insert(trackHandle2, estimatePix.begin(), estimatePix.end());
  fillerStr.insert(trackHandle2, estimateStr.begin(), estimateStr.end());
  fillerAll.insert(trackHandle2, estimateAll.begin(), estimateAll.end());

  fillerPix.fill();
  fillerStr.fill();
  fillerAll.fill();

  // Put back result to event
  ev.put(std::move(outputPix), "energyLossPixHits");
  ev.put(std::move(outputStr), "energyLossStrHits");
  ev.put(std::move(outputAll), "energyLossAllHits");
}

