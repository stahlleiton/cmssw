//Headers for the core items
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Headers for the data items
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

//Headers for services and tools
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "PIDCutParam.h"
#include "TF1.h"


class TrackPIDSelector : public edm::EDProducer {
public:
  explicit TrackPIDSelector(const edm::ParameterSet&);
  ~TrackPIDSelector();

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // Cut info
  const std::string particle_;
  const double maxInvBetaSignificance_;
  const double maxDeDxSignificance_;

  // Input info
  const edm::EDGetTokenT< reco::TrackCollection  > recoTracksToken_;
  const edm::EDGetTokenT< reco::VertexCollection > primaryVertexToken_;
  const edm::EDGetTokenT< reco::BeamSpot         > beamSpotToken_;
  const edm::EDGetTokenT< edm::ValueMap< reco::DeDxData > > dEdxToken_;
  std::map<std::string, edm::EDGetTokenT< edm::ValueMap< float > > > trackInfoToken_;

  // String selector
  const StringCutObjectSelector<reco::Track, true> trackSelection_;

  // Cut function
  std::map<std::string, std::map<std::string, std::map<std::string, TF1> > > cutFunction_;

  // Constants
  const float c_cm_ns   = 2.99792458e1; //[cm/ns]
  const std::vector<float> MASS = {0.13957018, 0.493677, 0.9382720813, 1.87561, 2.80925, 2.80923/2., 3.72742/2.};
  const std::vector<std::string> NAME = {"Pion", "Kaon", "Proton", "Deuteron", "Triton", "Helium3", "Helium4"};
};


// Constructor
TrackPIDSelector::TrackPIDSelector(const edm::ParameterSet& iConfig) :
  particle_(iConfig.getParameter<std::string>("particle")),
  maxInvBetaSignificance_(iConfig.getParameter<double>("maxInvBetaSignificance")),
  maxDeDxSignificance_(iConfig.getParameter<double>("maxDeDxSignificance")),
  recoTracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("recoTracksTag"))),
  primaryVertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
  dEdxToken_(consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("dEdxTag"))),
  trackSelection_(iConfig.getParameter<std::string>("trackSelection"))
{
  // Initialize track information
  for (const auto& s : {"TMTD", "SigmaTMTD", "PathLength"}) {
    trackInfoToken_[s] = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>(Form("track%sTag", s)));
  }

  // Initialize cut function
  for (const auto& e : {"BTL", "ETL"}) {
    for (const auto& p : {"invBeta", "dEdx"}) {
      cutFunction_[e][p].emplace("Mean", TF1("", FUNCMAP_[e][particle_][p]["Mean"].c_str(), 0.7, 200.));
      cutFunction_[e][p].emplace("Error", TF1("", FUNCMAP_[e][particle_][p]["Error"].c_str(), 0.7, 200.));
    }
  }

  // Check particle name
  if (std::find(NAME.begin(), NAME.end(), particle_) == NAME.end()) {
    std::string msg = "Particle label: "+particle_+" not valid!\nOptions are: ";
    for (const auto& p : NAME) msg += p+" , ";
    edm::LogError("TrackPIDSelector_InvalidLabel") << msg << std::endl; return;
  }

  produces<reco::TrackCollection>();
}


// Destructor
TrackPIDSelector::~TrackPIDSelector()
{
}


//
// Methods
//

// Producer method
void TrackPIDSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Extract primary vertices
  reco::VertexCollection vertices;
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(primaryVertexToken_, vertexHandle);
  for (const auto& pv : *vertexHandle) {
    if (pv.isFake() || pv.tracksSize()<2) continue;
    vertices.push_back(pv);
  }
  if (vertices.empty()) {
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    vertices.push_back(reco::Vertex(beamSpotHandle->position(), beamSpotHandle->rotatedCovariance3D()));
  }

  // Extract tracks
  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByToken(recoTracksToken_, trackHandle);
  if (!trackHandle.isValid()) { edm::LogError("TrackPIDSelector_ProductNotValid") << "General track collection is not valid!" << std::endl; return; }

  // Extract track MTD information
  std::map<std::string, edm::ValueMap<float> > trackInfoMap;
  for (const auto& s : trackInfoToken_) {
    edm::Handle<edm::ValueMap<float> > trackInfoHandle;
    iEvent.getByToken(s.second, trackInfoHandle);
    if (!trackInfoHandle.isValid()) { edm::LogError("TrackPIDSelector_ProductNotValid") << s.first << " valueMap is not valid!" << std::endl; return; }
    trackInfoMap[s.first] = *trackInfoHandle;
  }

  // Extract track dEdx information
  edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;
  iEvent.getByToken(dEdxToken_, dEdxHandle);
  const auto& dEdxData = *dEdxHandle;

  // Select tracks
  reco::TrackCollection selTracks;
  for (size_t iTrk = 0; iTrk < trackHandle->size(); iTrk++) {
    const auto& track = reco::TrackRef(trackHandle, iTrk);

    // Select tracks based on user string
    if (!trackSelection_(*track)) continue;

    const auto p = track->p();
    const auto pT = p / std::cosh(track->eta());
    const auto absEta = std::abs(track->eta());
    const auto tMTDErr = trackInfoMap["SigmaTMTD"][track];
    const bool hasMTD = (tMTDErr > 0) && (absEta < 1.4 ? pT > 0.8 : p > 0.7);

    // Select tracks based on dEdx PID
    if (maxDeDxSignificance_ > 0) {
      const auto dEdx = dEdxData[track].dEdx();
      const auto dEdx_Mean = (absEta < 1.6 ? cutFunction_["BTL"] : cutFunction_["ETL"])["dEdx"]["Mean"].Eval(p);
      const auto dEdx_Error = (absEta < 1.6 ? cutFunction_["BTL"] : cutFunction_["ETL"])["dEdx"]["Error"].Eval(p);
      const auto dEdxSignificance = std::abs(dEdx - dEdx_Mean)/dEdx_Error;
      if (dEdxSignificance > maxDeDxSignificance_) continue;
    }
    else if (maxInvBetaSignificance_ > 0 && !hasMTD) continue;

    // Select tracks based on MTD PID
    if (maxInvBetaSignificance_ > 0 && hasMTD) {
      // Extract MTD information
      const auto tMTD = trackInfoMap["TMTD"][track];
      const auto pathLength = trackInfoMap["PathLength"][track];

      // Associate primary vertex
      reco::Vertex vertex;
      if (vertices.empty()) break;
      else if (vertices.size()==1) vertex = vertices[0];
      else if (maxInvBetaSignificance_ > 0) {
        const auto rsigmazsq = 1./track->dzError()/track->dzError();
        // First try based on association weights, but keep track of closest in z and z-t as well
        reco::Vertex vtxMindz, vtxMinchisq;
        double mindz = 0.1, minchisq = std::numeric_limits<double>::max(), maxW = 0.5;
        for (const auto& vtx : vertices) {
          const auto w = vtx.trackWeight(track);
          if (w > maxW) {
            maxW = w;
            vertex = vtx;
          }
          const auto dz = std::abs(track->dz(vtx.position()));
          if (dz < mindz) {
            mindz = dz;
            vtxMindz = vtx;
          }
          if (vtx.tError() > 0. && vtx.tError() < 0.025) {
            for (const auto& m : MASS) {
              const auto invBeta = std::sqrt(m*m + p*p)/p;
              const auto t0 = tMTD - (pathLength*invBeta/c_cm_ns);
              const auto dt = std::abs(t0 - vtx.t());
              const auto dtsig = dt/tMTDErr;
              const auto chisq = dz*dz*rsigmazsq + dtsig*dtsig;
              if (chisq < minchisq) {
                minchisq = chisq;
                vtxMinchisq = vtx;
              }
            }
          }
        }
        // If no vertex found based on association weights, fall back to closest in z or z-t
        if (!vertex.isValid()) {
          // If closest vertex in z does not have valid time information, just use it, 
          // otherwise use the closest vertex in z-t plane with timing info, with a fallback to the closest in z
          const bool isBadT = (vtxMindz.isValid() && !(vtxMindz.tError() > 0. && vtxMindz.tError() < 0.025));
          if (isBadT) vertex = vtxMindz;
          else if (vtxMinchisq.isValid()) vertex = vtxMinchisq;
          else if (vtxMindz.isValid()) vertex = vtxMindz;
        }
        if (!vertex.isValid()) continue; 
      }
    
      // Compute beta
      const auto t0_PV = vertex.t();
      const auto dt_PV = tMTD - t0_PV;
      const auto beta_PV = (pathLength/dt_PV)*(1./c_cm_ns);
      const auto invBeta = 1/beta_PV;

      // Compute cut significance
      const auto invBeta_Mean = (absEta < 1.6 ? cutFunction_["BTL"] : cutFunction_["ETL"])["invBeta"]["Mean"].Eval(p);
      const auto invBeta_Error = (absEta < 1.6 ? cutFunction_["BTL"] : cutFunction_["ETL"])["invBeta"]["Error"].Eval(p);
      const auto invBetaSignificance = std::abs(invBeta - invBeta_Mean)/invBeta_Error;
      if (invBetaSignificance > maxInvBetaSignificance_) continue;
    }

    // Fill output track collection
    selTracks.push_back(*track);
  }

  // Store output tracks
  auto output = std::make_unique<reco::TrackCollection>(selTracks);
  output->shrink_to_fit();
  iEvent.put(std::move(output));
}


void TrackPIDSelector::beginJob()
{
}


void TrackPIDSelector::endJob()
{
}


// Define as plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(TrackPIDSelector);
