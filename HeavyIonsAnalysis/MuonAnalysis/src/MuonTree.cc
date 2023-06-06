// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

// data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//services and tools
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// root include files
#include "TTree.h"


//
// class declaration
//
class MuonTree : public edm::EDAnalyzer
{
 public:
  explicit MuonTree(const edm::ParameterSet&);
  ~MuonTree();

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  const edm::EDGetTokenT< pat::MuonCollection > tokRecoMu;
  const edm::EDGetTokenT< reco::VertexCollection > tokVtx;
  const edm::EDGetTokenT< reco::GenParticleCollection > tokGenPtl;
  const edm::EDGetTokenT<edm::TriggerResults> tokTrigRes;

  TTree* treeMu;
  edm::Service<TFileService> fs;

  std::map< std::string, int > varI;
  std::map< std::string, float > varF;
  std::map< std::string, std::vector<int> > varIV;
  std::map< std::string, std::vector<float> > varFV;
};


//
// constructors and destructor
//
MuonTree::MuonTree(const edm::ParameterSet& iConfig) :
  tokRecoMu(consumes< pat::MuonCollection >(iConfig.getParameter<edm::InputTag>("muons"))),
  tokVtx(consumes< reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"))),
  tokGenPtl(consumes< reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("genparticles"))),
  tokTrigRes(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"))),
  treeMu(0)
{
}


MuonTree::~MuonTree()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Extract magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  const auto& magField = bFieldHandle.product();

  // Extract gen muons
  std::vector<reco::GenParticleRef> genMuons;
  std::map<reco::GenParticleRef, size_t> genIdx, genToRecIdx;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(tokGenPtl, genParticles);
  if (genParticles.isValid()) {
    for (size_t i=0; i<genParticles->size(); i++) {
      const auto& gp = reco::GenParticleRef(genParticles, i);
      if (abs(gp->pdgId())==13 && gp->status()==1 && gp->isLastCopy()) {
        genMuons.push_back(gp);
        genIdx[gp] = genMuons.size()-1;
      }
    }
    if (!treeMu && genMuons.empty())
      genMuons.emplace_back(reco::GenParticleRef());
  }

  // Fill event information
  varI["event"] = iEvent.id().event();
  varI["run"] = iEvent.id().run();
  varI["lumi"] = iEvent.id().luminosityBlock();

  // Add trigger information
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(tokTrigRes, triggerResults);
  if (triggerResults.isValid()) {
    varI["trigger_fired"] = 0;
    const auto& triggerNames = iEvent.triggerNames(*triggerResults);
    for (ushort trgIdx=0; trgIdx<triggerNames.size(); trgIdx++) {
      if (triggerNames.triggerName(trgIdx).find("HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v")!=std::string::npos &&
          triggerResults->wasrun(trgIdx)) {
        varI.at("trigger_fired") = triggerResults->accept(trgIdx);
        break;
      }
    }
  }

  // Fill vertex information
  math::XYZPoint vertex;
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tokVtx, vertices);
  if (vertices.isValid()) {
    varI["vertex_size"] = vertices->size();
    varF["vertex_x"] = (!vertices->empty() ? vertices->at(0).x() : 99.);
    varF["vertex_y"] = (!vertices->empty() ? vertices->at(0).y() : 99.);
    varF["vertex_z"] = (!vertices->empty() ? vertices->at(0).z() : 99.);
    vertex = math::XYZPoint(varF.at("vertex_x"), varF.at("vertex_y"), varF.at("vertex_z"));
  }

  // Fill reco muon information
  edm::Handle<pat::MuonCollection> muonsH;
  iEvent.getByToken(tokRecoMu, muonsH);
  if (muonsH.isValid()) {
    const auto& muons = (treeMu || !muonsH->empty()) ? *muonsH : pat::MuonCollection({pat::Muon()});
    for (size_t i=0; i<muons.size(); i++) {
      const auto& muon = muons[i];

      // Kinematic information
      varFV["rec_muon_pt"].push_back(muon.pt());
      varFV["rec_muon_eta"].push_back(muon.eta());
      varFV["rec_muon_phi"].push_back(muon.phi());
      varIV["rec_muon_charge"].push_back(muon.charge());

      // Tight ID Muon POG Run 2
      varIV["rec_muon_isGlobal"].push_back(muon.isGlobalMuon());
      varIV["rec_muon_isPF"].push_back(muon.isPFMuon());
      varFV["rec_muon_global_normChi2"].push_back(muon.globalTrack().isNonnull() ? muon.globalTrack()->normalizedChi2() : 99.);
      varIV["rec_muon_global_nMuonHits"].push_back(muon.globalTrack().isNonnull() ? muon.globalTrack()->hitPattern().numberOfValidMuonHits() : -1);
      varIV["rec_muon_nMuonStations"].push_back(muon.numberOfMatchedStations());
      varIV["rec_muon_track_nPixelHits"].push_back(muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().numberOfValidPixelHits() : -1);
      varIV["rec_muon_track_nTrackerLayers"].push_back(muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1);
      varFV["rec_muon_track_dXY"].push_back(muon.innerTrack().isNonnull() && !vertices->empty() ? muon.innerTrack()->dxy(vertex) : 99.);
      varFV["rec_muon_track_dZ"].push_back(muon.innerTrack().isNonnull() && !vertices->empty() ? muon.innerTrack()->dz(vertex) : 99.);
      varIV["rec_muon_isTight"].push_back(varIV.at("rec_muon_isGlobal").back()>0 &&
                                          varIV.at("rec_muon_isPF").back()>0 &&
                                          varFV.at("rec_muon_global_normChi2").back()<10. &&
                                          varIV.at("rec_muon_global_nMuonHits").back()>0 &&
                                          varIV.at("rec_muon_nMuonStations").back()>1 &&
                                          varIV.at("rec_muon_track_nPixelHits").back()>0 &&
                                          varIV.at("rec_muon_track_nTrackerLayers").back()>5 &&
                                          abs(varFV.at("rec_muon_track_dXY").back())<0.2 &&
                                          abs(varFV.at("rec_muon_track_dZ").back())<0.5);

      // Soft ID Muon POG Run 2
      varIV["rec_muon_isTracker"].push_back(muon.isTrackerMuon());
      varIV["rec_muon_oneStationTight"].push_back(muon::isGoodMuon(muon, muon::SelectionType::TMOneStationTight));
      varIV["rec_muon_track_nPixelLayers"].push_back(muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() : -1);
      varIV["rec_muon_track_isHighPurity"].push_back(muon.innerTrack().isNonnull() ? muon.innerTrack()->quality(reco::TrackBase::highPurity) : 0);
      varIV["rec_muon_isSoft"].push_back(varIV.at("rec_muon_isTracker").back()>0 &&
                                         varIV.at("rec_muon_oneStationTight").back()>0 &&
                                         varIV.at("rec_muon_track_nTrackerLayers").back()>5 &&
                                         varIV.at("rec_muon_track_nPixelLayers").back()>0 &&
                                         varIV.at("rec_muon_track_isHighPurity").back()>0 &&
                                         abs(varFV.at("rec_muon_track_dXY").back())<0.3 &&
                                         abs(varFV.at("rec_muon_track_dZ").back())<20.);

      // Hybrid Soft ID HIN PAG Run 2 PbPb
      varIV["rec_muon_isHybridSoft"].push_back(varIV.at("rec_muon_isGlobal").back()>0 &&
                                               varIV.at("rec_muon_track_nTrackerLayers").back()>5 &&
                                               varIV.at("rec_muon_track_nPixelLayers").back()>0 &&
                                               abs(varFV.at("rec_muon_track_dXY").back())<0.3 &&
                                               abs(varFV.at("rec_muon_track_dZ").back())<20.);

      // Muon Trigger Matching
      varIV["rec_muon_fireL1Trigger"].push_back(muon.triggerObjectMatchesByFilter("hltL1sSingleMuOpenNotMBHF2AND").size()>0);

      // Muon Gen Matching
      const auto& gen = muon.genParticleRef();
      varIV["rec_muon_genIdx"].push_back(gen.isNonnull() ? genIdx.at(gen) : -1);
      if (gen.isNonnull())
        genToRecIdx[gen] = i;
    }
  }

  // Fill gen muon information
  for (const auto& gen : genMuons) {
    const auto& muon = gen.isNonnull() ? *gen : reco::GenParticle();
    // Kinematic information
    varFV["gen_muon_pt"].push_back(muon.pt());
    varFV["gen_muon_eta"].push_back(muon.eta());
    varFV["gen_muon_phi"].push_back(muon.phi());
    varIV["gen_muon_charge"].push_back(muon.charge());
    varIV["gen_muon_mass"].push_back(muon.mass());

    // Muon Reco Matching
    const auto& recIt = genToRecIdx.find(gen);
    varIV["gen_muon_recIdx"].push_back(recIt!=genToRecIdx.end() ? recIt->second : -1);
  }

  for (const auto& v : varFV)
    for (const auto& var : v.second)
      std::cout << v.first << " >> " << var << std::endl;
  std::cout << std::endl;

  // Fill reco dimuon information
  if (muonsH.isValid()) {
    const auto& muons = (treeMu || muonsH->size()>1) ? *muonsH : pat::MuonCollection({pat::Muon(), pat::Muon()});
    for (size_t i1=0; i1<muons.size(); i1++) {
      for (size_t i2=i1+1; i2<muons.size(); i2++) {
        const auto& muon1 = muons[i1];
        const auto& muon2 = muons[i2];

        // Kinematic information
        const auto p4 = muon1.p4() + muon2.p4();
        varFV["rec_dimuon_pt"].push_back(p4.pt());
        varFV["rec_dimuon_eta"].push_back(p4.eta());
        varFV["rec_dimuon_rapidity"].push_back(p4.Rapidity());
        varFV["rec_dimuon_phi"].push_back(p4.phi());
        varFV["rec_dimuon_mass"].push_back(p4.mass());
        varIV["rec_dimuon_charge"].push_back(muon1.charge()+muon2.charge());

        // Muon indices
        varIV["rec_dimuon_muon1Idx"].push_back(i1);
        varIV["rec_dimuon_muon2Idx"].push_back(i2);

        // Perform vertex fit
        RefCountedKinematicVertex vtx;
        if (muon1.track().isNonnull() && muon2.track().isNonnull()) {
          float chi = 0.0 , ndf = 0.0, tmp(1E-7);
          std::vector<RefCountedKinematicParticle> muonP;
          KinematicParticleFactoryFromTransientTrack pFactory;
          muonP.push_back(pFactory.particle(reco::TransientTrack(*muon1.track(), magField), muon1.mass(), chi, ndf, tmp));
          muonP.push_back(pFactory.particle(reco::TransientTrack(*muon2.track(), magField), muon2.mass(), chi, ndf, tmp));
          const auto res = KinematicParticleVertexFitter().fit(muonP);
          if (res->isValid())
            vtx = res->currentDecayVertex();
        }
        varFV["rec_dimuon_vertex_prob"].push_back(vtx && vtx->vertexIsValid() ? TMath::Prob(vtx->chiSquared(), vtx->degreesOfFreedom()) : -1.);
        varFV["rec_dimuon_vertex_normChi2"].push_back(vtx && vtx->vertexIsValid() ? vtx->chiSquared()/vtx->degreesOfFreedom() : -1.);
      }
    }
  }

  // Fill tree
  if (!treeMu) {
    treeMu = fs->make<TTree>("MuonTree","MuonTree");
    for (auto& p : varI)
      treeMu->Branch(p.first.c_str(), &(p.second), (p.first+"/I").c_str());
    for (auto& p : varF)
      treeMu->Branch(p.first.c_str(), &(p.second), (p.first+"/F").c_str());
    for (auto& p : varIV)
      treeMu->Branch(p.first.c_str(), &(p.second));
    for (auto& p : varFV)
      treeMu->Branch(p.first.c_str(), &(p.second));
  }
  treeMu->Fill();

  // Clear information
  for (auto& p : varI) p.second = -1;
  for (auto& p : varF) p.second = 99.;
  for (auto& p : varIV) p.second.clear();
  for (auto& p : varFV) p.second.clear();
}


// ------------ method called once each job just before starting event loop  ------------
void
MuonTree::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonTree::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTree);
