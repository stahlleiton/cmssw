#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEXMLStringProduct.h"

#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include <fstream>
#include <iostream>
#include <map>

class GenParticles2LHEConverter : public edm::one::EDProducer<edm::BeginRunProducer, edm::EndRunProducer>
{
 public:
  explicit GenParticles2LHEConverter(const edm::ParameterSet& pset);
  ~GenParticles2LHEConverter() override{};

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

 private:

  void beginRunProduce(edm::Run&, const edm::EventSetup&) override;
  void endRunProduce(edm::Run&, const edm::EventSetup&) override;
  void produce(edm::Event& event, const edm::EventSetup&) override;

  size_t addParticle(lhef::HEPEUP&, const reco::GenParticleRef&);

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<GenRunInfoProduct> genRunInfoToken_;
  edm::ESHandle<ParticleDataTable> pTable_;

  const double cmEnergy_;
  const int beamID_;

  std::ofstream file;
  std::map<reco::GenParticleRef, size_t> map_;
};

GenParticles2LHEConverter::GenParticles2LHEConverter(const edm::ParameterSet& pset)
    // dummy value to set incident proton pz for particle gun samples
    : cmEnergy_(pset.getUntrackedParameter<double>("cmEnergy", 13000)),
      beamID_(pset.getUntrackedParameter<int>("beamID", 2212))
{
  genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticles"));
  genEventInfoToken_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo"));
  genRunInfoToken_ = consumes<GenRunInfoProduct, edm::InRun>(pset.getParameter<edm::InputTag>("genEventInfo"));

  produces<LHEEventProduct>();
  produces<LHERunInfoProduct, edm::Transition::BeginRun>();

  const auto& lhe_ouputfile = pset.getUntrackedParameter<std::string>("fileName", "");
  if (!lhe_ouputfile.empty()) {
    file.open(lhe_ouputfile, std::fstream::out | std::fstream::trunc);
  }
};

void GenParticles2LHEConverter::produce(edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;
  iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);

  iEventSetup.getData(pTable_);

  lhef::HEPEUP hepeup;
  
  // add gen event information
  hepeup.IDPRUP = genEventInfoHandle->signalProcessID();
  hepeup.XWGTUP = genEventInfoHandle->weight(); 
  hepeup.SCALUP = genEventInfoHandle->qScale();
  hepeup.AQEDUP = genEventInfoHandle->alphaQED();
  hepeup.AQCDUP = genEventInfoHandle->alphaQCD();
  hepeup.XPDWUP = (genEventInfoHandle->pdf() ? genEventInfoHandle->pdf()->xPDF : std::pair<double,double>({}));

  map_.clear();
  for (size_t i=0; i<genParticlesHandle->size(); i++)
    addParticle(hepeup, reco::GenParticleRef(genParticlesHandle, i));

  double originalXWGTUP_ = 1.;
  std::unique_ptr<LHEEventProduct> product(new LHEEventProduct(hepeup, originalXWGTUP_));

  if (file.is_open())
    std::copy(product->begin(), product->end(), std::ostream_iterator<std::string>(file));

  iEvent.put(std::move(product));
};

size_t GenParticles2LHEConverter::addParticle(lhef::HEPEUP& outlhe, const reco::GenParticleRef& par)
{
  if (par.isNull()) return 0;
  if (map_.find(par)!=map_.end()) return map_.at(par);
  size_t parIdx = outlhe.NUP;
  outlhe.resize(outlhe.NUP + 1);

  double mass(0), spin(0), ctau(0), color(0);
  if (pTable_->particle(par->pdgId())) {
    const auto& p = pTable_->particle(par->pdgId());
    mass = p->mass().value();
    spin = p->spin().spin();
    ctau = p->lifetime().value();
    color = p->color();
  }
  else {
    mass = par->mass();
  } 

  outlhe.PUP[parIdx][0] = par->px();
  outlhe.PUP[parIdx][1] = par->py();
  outlhe.PUP[parIdx][2] = par->pz();
  outlhe.PUP[parIdx][3] = par->energy();
  outlhe.PUP[parIdx][4] = mass;
  outlhe.IDUP[parIdx]   = par->pdgId();
  outlhe.SPINUP[parIdx] = spin;
  outlhe.VTIMUP[parIdx] = ctau;
  outlhe.ISTUP[parIdx]  = (par->status()==3 ? 2 : par->status());
  outlhe.ICOLUP[parIdx].first  = color;
  outlhe.ICOLUP[parIdx].second = color;

  size_t firstMother(0), lastMother(0);
  for (size_t i=0; i<par->numberOfMothers(); i++) {
    const auto& idx = addParticle(outlhe, par);
    if (idx==0) continue;
    firstMother = (firstMother>0 ? std::min(firstMother, idx) : idx);
    lastMother = (lastMother>0 ? std::max(lastMother, idx) : idx);
  }
  outlhe.MOTHUP[parIdx].first = firstMother;
  outlhe.MOTHUP[parIdx].second = lastMother;

  map_[par] = parIdx;

  return parIdx;
};

void GenParticles2LHEConverter::beginRunProduce(edm::Run& iRun, const edm::EventSetup&)
{
  edm::Handle<GenRunInfoProduct> genRunInfoHandle;
  iRun.getByToken(genRunInfoToken_, genRunInfoHandle);

  // fill HEPRUP common block and store in edm::Run
  lhef::HEPRUP heprup;

  // set number of processes: 1 for Z to tau tau
  heprup.resize(1);

  //Process independent information

  //beam particles ID
  //heprup.IDBMUP.first = beamID_;
  //heprup.IDBMUP.second = beamID_;

  //beam particles energies
  //heprup.EBMUP.first = cmEnergy_/2.;
  //heprup.EBMUP.second = cmEnergy_/2.;

  //take default pdf group for both beamparticles
  //heprup.PDFGUP.first = -1;
  //heprup.PDFGUP.second = -1;

  //take certan pdf set ID (same as in officially produced DYJets LHE files)
  //heprup.PDFSUP.first = -1;
  //heprup.PDFSUP.second = -1;

  heprup.IDWTUP = 3;
  auto xsec = genRunInfoHandle->internalXSec().value();
  heprup.XSECUP[0] = (xsec==0. ? 1. : xsec);
  heprup.XERRUP[0] = genRunInfoHandle->internalXSec().error();
  heprup.XMAXUP[0] = 1;
  heprup.LPRUP[0] = 1;

  std::unique_ptr<LHERunInfoProduct> runInfo(new LHERunInfoProduct(heprup));

  if (file.is_open())
    std::copy(runInfo->begin(), runInfo->end(), std::ostream_iterator<std::string>(file));
  iRun.put(std::move(runInfo));
};

void GenParticles2LHEConverter::endRunProduce(edm::Run& run, const edm::EventSetup& es)
{
  if (file.is_open()) {
    file << LHERunInfoProduct::endOfFile();
    file.close();
  }
};

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenParticles2LHEConverter::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
};

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticles2LHEConverter);
