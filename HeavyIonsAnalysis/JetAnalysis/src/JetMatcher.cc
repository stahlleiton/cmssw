#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

class JetMatcher : public edm::global::EDProducer<> {
public:
  typedef edm::AssociationMap<edm::OneToOne<reco::JetView, reco::JetView> > JetMatchMap;
  explicit JetMatcher(const edm::ParameterSet& iConfig)
      : sourceToken_(consumes<reco::JetView>(iConfig.getParameter<edm::InputTag>("source"))),
        matchedToken_(consumes<reco::JetView>(iConfig.getParameter<edm::InputTag>("matched"))) {
    produces<JetMatchMap>();
  };
  ~JetMatcher() override{};
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

private:
  const edm::EDGetTokenT<reco::JetView> sourceToken_;
  const edm::EDGetTokenT<reco::JetView> matchedToken_;
};

void JetMatcher::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  const auto& source = iEvent.getHandle(sourceToken_);
  const auto& matched = iEvent.getHandle(matchedToken_);
 
  auto matching = std::make_unique<JetMatchMap>(source, matched);
  for (size_t i = 0; i < source->size(); i++) {
    std::pair<float, int> m(999, -1);
    for (size_t j = 0; j < matched->size(); j++) {
      const auto dR = deltaR((*source)[i].p4(), (*matched)[j].p4());
      if (dR < m.first)
        m = {dR, j};
    }
    matching->insert(source->refAt(i), m.second>=0 ? matched->refAt(m.second) : edm::RefToBase<reco::Jet>());
  }
  iEvent.put(std::move(matching));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetMatcher);
