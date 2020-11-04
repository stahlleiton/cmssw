#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <iostream>
#include <map>

class HepMCMerger : public edm::stream::EDProducer<>
{
 public:
  explicit HepMCMerger(const edm::ParameterSet& pset);
  ~HepMCMerger() override{};

 private:

  void produce(edm::Event& event, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::HepMCProduct> sigToken_;
  edm::EDGetTokenT<edm::HepMCProduct> bkgToken_;
  std::string label_;
};

HepMCMerger::HepMCMerger(const edm::ParameterSet& pset)
{
  label_ = pset.getUntrackedParameter<std::string>("label", "unsmeared");
  sigToken_ = consumes<edm::HepMCProduct>(edm::InputTag(pset.getParameter<std::string>("signal"), label_));
  bkgToken_ = consumes<edm::HepMCProduct>(edm::InputTag(pset.getParameter<std::string>("background"), label_));

  produces<edm::HepMCProduct>(label_);
};

void HepMCMerger::produce(edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{
  edm::Handle<edm::HepMCProduct> sigHandle;
  iEvent.getByToken(sigToken_, sigHandle);

  edm::Handle<edm::HepMCProduct> bkgHandle;
  iEvent.getByToken(bkgToken_, bkgHandle);

  auto bkgPrd = std::make_unique<edm::HepMCProduct>(*bkgHandle);

  if (label_!="unsmeared") {
    // Step 1: Get signal event vertex
    const auto& sigEvt = sigHandle->GetEvent();
    auto sigVtx = sigEvt->signal_process_vertex();
    for (auto pt=sigEvt->particles_begin(); sigVtx || pt!=sigEvt->particles_end(); ++pt) {
      sigVtx = (*pt)->production_vertex();
    }
    if (!sigVtx) return;
    const auto& sigPos = sigVtx->position();

    // Step 2: Get background event vertex
    const auto& bkgEvt = bkgHandle->GetEvent();
    auto bkgVtx = bkgEvt->signal_process_vertex();
    for (auto pt=bkgEvt->particles_begin(); bkgVtx || pt!=bkgEvt->particles_end(); ++pt) {
      bkgVtx = (*pt)->production_vertex();
    }
    const auto& bkgPos = (bkgVtx ? bkgVtx->position() : HepMC::FourVector());
 
    // Step 3: Shift background event to signal vertex
    const auto& dV = HepMC::FourVector(sigPos.x() - bkgPos.x(), sigPos.y() - bkgPos.y(), sigPos.z() - bkgPos.z(), sigPos.t() - bkgPos.t());
    for (auto vt=bkgPrd->GetEvent()->vertices_begin(); vt!=bkgPrd->GetEvent()->vertices_end(); ++vt) {
      const auto& v = (*vt)->position();
      (*vt)->set_position(HepMC::FourVector(v.x()+dV.x(), v.y()+dV.y(), v.z()+dV.z(), v.t()+dV.t()));
    }
  }

  // Add vertices to output event
  const auto& oEvt = new HepMC::GenEvent(*bkgPrd->GetEvent());
  for (auto vt=sigHandle->GetEvent()->vertices_begin(); vt!=sigHandle->GetEvent()->vertices_end(); ++vt) {
    oEvt->add_vertex(new HepMC::GenVertex(**vt));
  }

  // Store output
  auto output = std::make_unique<edm::HepMCProduct>(*bkgPrd);
  if (output->GetEvent()) delete output->GetEvent();
  output->addHepMCData(oEvt);
  iEvent.put(std::move(output), label_);
};

//define this as a plug-in
DEFINE_FWK_MODULE(HepMCMerger);
