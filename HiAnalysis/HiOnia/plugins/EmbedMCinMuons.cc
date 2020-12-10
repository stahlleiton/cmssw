#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

namespace pat {

  class EmbedMCinMuons : public edm::stream::EDProducer<> {
    
  public:
    
    explicit EmbedMCinMuons(const edm::ParameterSet & iConfig):
      muonGenMatchToken_(consumes(iConfig.getParameter<edm::InputTag>("matchedGen"))),//type from CommonTools/UtilAlgos/interface/PhysObjectMatcher.h //not mentionning the type <edm::Association<reco::GenParticleCollection>> to consumes works from 11_2_0_pre6
      muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
    {
      //for (const auto& sel: muonSelectors_) {
      //candidateMuonIDToken_["packedPFCandidate"][sel] = consumes<pat::PackedCandidateRefVector>(edm::InputTag("packedCandidateMuonID", "pfCandidates"+sel));
      //candidateMuonIDToken_["lostTrack"][sel] = consumes<pat::PackedCandidateRefVector>(edm::InputTag("packedCandidateMuonID", "lostTracks"+sel));
      //}
      produces<pat::MuonCollection>();
    }
    ~EmbedMCinMuons() override {};
    
    void produce(edm::Event&, const edm::EventSetup&) override;
    
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    
  private:
    const edm::EDGetTokenT<edm::Association<reco::GenParticleCollection> > muonGenMatchToken_;
    const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    
  };
  
}

void pat::EmbedMCinMuons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // extract MC-muons matching map
  edm::Handle<edm::Association<reco::GenParticleCollection> > matches;
  iEvent.getByToken(muonGenMatchToken_, matches);

  // extract muons from packedCandidate collections
  edm::Handle<edm::View<pat::Muon> > muonCands;
  iEvent.getByToken(muonToken_, muonCands);
  unsigned int ncand = muonCands->size();

  pat::MuonCollection outMuons;

  auto output = std::make_unique<pat::MuonCollection>();
  //loop over muons to grab and add the matched gen muon
  for (size_t i=0; i<ncand; i++) {
    edm::RefToBase<pat::Muon> muonref(muonCands, i);
    reco::GenParticleRef genPRef = (*matches)[muonref]; //PhysicsTools/NanoAOD/plugins/CandMCMatchTableProducer.cc, source mcMap is like muonMatch

    //replacing the associated gen particle
    pat::Muon outmuon = *muonref;
    outmuon.setGenParticleRef(genPRef, true); //2nd argument is for embedding the particle
    
    output->push_back(outmuon);
  }

  iEvent.put(std::move(output));
}
  
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void pat::EmbedMCinMuons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("unpackedMuons"))->setComment("muon input collection");
  desc.add<edm::InputTag>("matchedGen", edm::InputTag("muonMatch"))->setComment("matches with gen muons input collection");
  descriptions.add("unpackedMuonsWithGenMatch", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace pat;
DEFINE_FWK_MODULE(EmbedMCinMuons);
