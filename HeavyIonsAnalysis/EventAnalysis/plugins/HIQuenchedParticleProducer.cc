#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TH1F.h"
#include "TF1.h"


class HIQuenchedParticleProducer : public edm::global::EDProducer<> {
  public:
    typedef edm::Ref<edm::View<reco::Jet>> JetRef;
    typedef std::vector<edm::InputTag> VInputTag;
    typedef std::vector<edm::EDGetTokenT<pat::PackedCandidateCollection>> PackedCandGetTokens;
    typedef std::vector<edm::EDPutTokenT<pat::PackedCandidateCollection>> PackedCandPutTokens;

    explicit HIQuenchedParticleProducer(const edm::ParameterSet& iConfig)
        : pvToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
          cenToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBins"))),
          jetToken_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
          packedCandidateGetTokens_(getTokens(iConfig.getParameter<VInputTag>("packedCandidates"), consumesCollector())),
          packedCandidatePutTokens_(putTokens(iConfig.getParameter<VInputTag>("packedCandidates"), producesCollector())),
          quenchingModel_(getQuenchingModel()),
          centralityModel_(getCentralityModel()) {}
    ~HIQuenchedParticleProducer() override{};

    void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    const edm::EDGetTokenT<reco::VertexCollection> pvToken_;
    const edm::EDGetTokenT<int> cenToken_;
    const edm::EDGetTokenT<edm::View<reco::Jet>> jetToken_;
    const PackedCandGetTokens packedCandidateGetTokens_;
    const PackedCandPutTokens packedCandidatePutTokens_;
    
    const std::unique_ptr<TH1F> quenchingModel_;
    const std::unique_ptr<TF1> centralityModel_;

    TH1F* getQuenchingModel() {
      TH1F *quenchingModel = new TH1F("qmh",";Quenching [GeV];AU",200,0,50);
      TF1 quenchingModelFunc("qmf", "[0]/(TMath::Sqrt(2.*TMath::Pi())*0.73*x)*TMath::Exp(-1.*TMath::Power(TMath::Log(x/[0])+1.5,2)/2./0.73/0.73)", 0., 50.);
      quenchingModelFunc.SetParameter(0, 50.); // this sets the omega_c parameter. if we want to make this centrality dependent
      for(int i=0; i<quenchingModel->GetNbinsX(); i++) {
        quenchingModel->SetBinContent(i+1, quenchingModelFunc.Eval(quenchingModel->GetBinCenter(i+1)));
        quenchingModel->SetBinError(i+1, 0);
      }
      quenchingModel->Scale(1./quenchingModel->Integral());
      return quenchingModel;
    };

    TF1* getCentralityModel() {
      TF1 *centralityModel = new TF1("centralityModel", "gaus");
      centralityModel->SetParameter(0,  1.090);
      centralityModel->SetParameter(1, -0.144);
      centralityModel->SetParameter(2,  0.442);
      return centralityModel;
    };

    PackedCandGetTokens getTokens(const VInputTag& tags, edm::ConsumesCollector&& iC) {
      PackedCandGetTokens output;
      for (const auto& tag : tags)
        output.emplace_back(iC.consumes<pat::PackedCandidateCollection>(tag));
      return output;
    }

    PackedCandPutTokens putTokens(const VInputTag& tags, edm::ProducesCollector iP) {
      PackedCandPutTokens output;
      for (const auto& tag : tags)
        output.emplace_back(iP.produces<pat::PackedCandidateCollection>(tag.label()+tag.instance()));
      return output;
    }
};

void HIQuenchedParticleProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // extract input information
  const auto& jetHandle = iEvent.getHandle(jetToken_);
  const auto& pvProd = reco::VertexRefProd(iEvent.getHandle(pvToken_));
  const auto& cenBin = iEvent.get(cenToken_);

  // evaluate centrality suppression
  const auto& centralitySuppression = centralityModel_->Eval(cenBin/200.);

  // extract jets
  std::vector<JetRef> jets;
  jets.reserve(jetHandle->size());
  std::map<JetRef, double> quenchPtRatio;
  for (size_t i=0; i<jetHandle->size(); i++) {
    const auto& jet = jets.emplace_back(jetHandle, i);
    // compute fraction of quenched jet momentum
    const auto dp = centralitySuppression * quenchingModel_->GetRandom();
    quenchPtRatio[jet] = std::max(jet->p() - dp, 0.) / jet->p();
  }
  std::sort(jets.begin(), jets.end(), [&](const JetRef& a, const JetRef& b) { return a->pt() > b->pt(); });

  // loop over packed candidates
  for (size_t i=0; i<packedCandidateGetTokens_.size(); i++) {
    const auto& packedCandidates = iEvent.get(packedCandidateGetTokens_[i]);
    pat::PackedCandidateCollection output;
    output.reserve(packedCandidates.size());
    for (const auto& cand : packedCandidates) {
      // find associated jet
      double quenchF(1.);
      for (const auto& jet : jets)
        if (reco::deltaR(*jet, cand) <= 0.41) {
          quenchF = quenchPtRatio.at(jet);
          break;
        }
      //apply to candidate
      if (quenchF == 1.)
        output.emplace_back(cand);
      else if (cand.pt()*quenchF > 0.1) {
        reco::Candidate::PolarLorentzVector quenchP4(cand.pt()*quenchF, cand.eta(), cand.phi(), cand.mass());
        auto& out = output.emplace_back(quenchP4, cand.vertex(), cand.ptTrk()*quenchF, cand.etaAtVtx(), cand.phiAtVtx(), cand.pdgId(), pvProd, cand.vertexRef().key());
        out.setAssociationQuality(cand.pvAssociationQuality());
        out.setCaloFraction(cand.caloFraction());
        out.setCovarianceVersion(cand.covarianceVersion());
        out.setFirstHit(cand.firstHit());
        out.setGoodEgamma(cand.isGoodEgamma());
        out.setHcalFraction(cand.hcalFraction());
        out.setIsIsolatedChargedHadron(cand.isIsolatedChargedHadron());
        out.setMuonID(cand.isStandAloneMuon(), cand.isGlobalMuon());
        out.setPuppiWeight(cand.puppiWeight(), cand.puppiWeightNoLep());
        out.setRawCaloFraction(cand.rawCaloFraction());
        out.setRawHcalFraction(cand.rawHcalFraction());
        out.setTime(cand.time(), cand.timeError());
        out.setTrackHighPurity(cand.trackHighPurity());
        out.setLostInnerHits(cand.lostInnerHits());
        out.setTrkAlgo(cand.trkAlgo(), cand.trkOriginalAlgo());
        if (cand.hasTrackDetails()) {
          auto track = cand.pseudoTrack();
          const_cast<reco::TrackBase::Vector&>(track.momentum()) = reco::TrackBase::Vector(track.px()*quenchF, track.py()*quenchF, track.pz());
          out.setTrackProperties(track, cand.covarianceSchema(), cand.covarianceVersion());
        }
      }
    }
    iEvent.emplace(packedCandidatePutTokens_[i], std::move(output));
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HIQuenchedParticleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("primaryVertices", {"offlineSlimmedPrimaryVertices"})->setComment("Primary Vertex input collection");
  desc.add<edm::InputTag>("centralityBins", {"centralityBin:HFtowers"})->setComment("Centrality bins");
  desc.add<edm::InputTag>("jets", {"akCs4PFUnquenchedJets"})->setComment("Jet input collection");
  desc.add<VInputTag>("packedCandidates", {edm::InputTag("packedPFCandidates"), edm::InputTag("lostTracks"), edm::InputTag("lostTracks:eleTracks")})->setComment("PF candidate input collections");
  descriptions.add("hiQuenchedParticles", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HIQuenchedParticleProducer);
