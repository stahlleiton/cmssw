// -*- C++ -*-
//
// Package:    VNAnalyzer
// Class:      VNAnalyzer
// 

// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/VNAnalysis/interface/TrackEfficiency.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>
	
#include <vector>
using std::vector;
using std::rand;
using namespace std;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
using namespace hi;
using namespace edm;


static const int ncentbins = 13;
static const  double centbins[]={0, 5, 10, 15, 20, 25, 30, 35, 40,  50, 60, 70, 80, 100};

static const int nptbins = 28;
static const float ptbins[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
			      2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 10., 12.0, 14.0, 16.0,  20.0, 26.0, 35.0, 45.0, 60.0, 80.0, 100., 200.};

static const int netabinsDefault = 12;
static const float etabinsDefault[]={-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
static const int nanals = 3;
enum AnalType {
        N2,          N3,         N4       
};
string AnalNames[]={
       "N2",          "N3",           "N4" 
};


//
// class declaration
//

class VNAnalyzer : public edm::EDAnalyzer {
public:
  explicit VNAnalyzer(const edm::ParameterSet&);
  ~VNAnalyzer();
      
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  // ----------member data ---------------------------
  int eporder_;
  int mx;
  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;
  edm::InputTag centralityBinTag_;
  edm::EDGetTokenT<int> centralityBinToken;
  edm::Handle<int> cbin_;
  edm::EDGetTokenT<int> tag_;

  edm::InputTag bTag_;
  edm::EDGetTokenT<double> bToken;
  edm::Handle<double> b_;
  double b;

  edm::InputTag RPTag_;
  edm::EDGetTokenT<double> RPToken;
  edm::Handle<double> RP_;
  double rp;

  edm::InputTag centralityTag_;
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;

  edm::InputTag pfTag_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfToken_;

  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::Service<TFileService> fs;

  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hb;
  TH1D * hRP;
  TH1D * hvtx;
  TH1D * hvtxRaw;
  TH1D * hcentbins;
  TH1D * hcentres;
  TH1D * hptNtrk;
  TH1D * hptNtrkGood;
  TH1I * hNtrkRet;
  //  TH2D * hEff[ncentbins];
  double centval;
  int ntrkval;
  double vtx;

  double reso_;
  int nvtx_;
  double minvz_;
  double maxvz_;
  double dzdzerror_pix_;
  double chi2_;

  double dzdzerror_;
  double d0d0error_;
  double pterror_;

  Double_t epang[NumEPNames];
  Double_t eporig[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];

  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t epmult[NumEPNames];
  Double_t sumw[NumEPNames];
  Double_t sumw2[NumEPNames];
  Double_t vn[NumEPNames];

  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  TH1D * hPsi[NumEPNames];
  TH1D * hPsiOffset[NumEPNames];
  TH1D * hPsiFlat[NumEPNames];


  unsigned int runno_;

  TH1D * hNtrk;
  int nEtaBins;
  string rpnames[NumEPNames];
  string effTable_;
  TTree * tree;
  TrackEfficiency *teff;
  int FlatOrder_;
  int NumFlatBins_;
  int CentBinCompression_;
  int flatnvtxbins_;
  double flatminvtx_;
  double flatdelvtx_;
  TH2D * qxtrk[7];
  TH2D * qytrk[7];
  TH2D * qcnt;
  TH2D * avpt;
  TH2D * res[7][ncentbins];
  TH2D * resw[7][ncentbins];
  TH2D * resep[7][ncentbins];
  TH2D * rescnt[7][ncentbins];
  
  HiEvtPlaneFlatten * flat[NumEPNames];
  bool loadDB_;
  bool useNcentbins_; 
  bool bypassCentrality_;
  bool FirstEvent_;
  bool makeTree_;
  int nCentBins_ = 1;
  int ntrack;


  //==============  Harmonics ============
  TH2D * ptav[ncentbins];
  TH2D * ptcnt[ncentbins];
  TH2D * qA[ncentbins][11];
  TH2D * qB[ncentbins][11];
  TH2D * v1MC[ncentbins];
  TH2D * v2MC[ncentbins];
  TH2D * v3MC[ncentbins];
  TH2D * vnCntMC[ncentbins];
  TH1D * qres;
  TH1D * thA[ncentbins][11];
  TH1D * thB[ncentbins][11];
  TH1D * thC[ncentbins][11];
  TH1D * thA2[ncentbins][11];
  TH1D * thB2[ncentbins][11];
  TH1D * thC2[ncentbins][11];
  TH1D * thn[ncentbins][11];

  TH2D * wnA[ncentbins][11];
  TH2D * wnB[ncentbins][11];
  TH2D * hTemplate;
  TH2D * qxt = 0;
  TH2D * qyt = 0;
  TH2D * qct = 0;
  TH2D * qxy = 0;
  TH2D * qxxy = 0;
  TH2D * qxyy = 0;
  TH2D * qx2y3 = 0;
  TH2D * qy2x3 = 0;
  TH2D * qx2x2x3;
  TH2D * qx2x3y2;
  TH2D * qx3y2y2;
  TH2D * qx2x2y3;
  TH2D * qx2y2y3;
  TH2D * qy2y2y3;
  TH2D * qcnt3;

  struct qvec {
    TH2D * qA[ncentbins][11];
    TH2D * qB[ncentbins][11];
    TH2D * wnA[ncentbins][11];
    TH2D * wnB[ncentbins][11];
    TH1D * thA[ncentbins][11];
    TH1D * thB[ncentbins][11];
    TH1D * thC[ncentbins][11];
    TH1D * thA2[ncentbins][11];
    TH1D * thB2[ncentbins][11];
    TH1D * thC2[ncentbins][11];
    TH1D * thn[ncentbins][11];
    TH1D * qBA[ncentbins][11];
    TH1D * qCA[ncentbins][11];
    TH1D * qCB[ncentbins][11];
    TH1D * qBAcnt[ncentbins][11];
    TH1D * qCAcnt[ncentbins][11];
    TH1D * qCBcnt[ncentbins][11];
    TH1D * qBA2[ncentbins][11];
    TH1D * qCA2[ncentbins][11];
    TH1D * qCB2[ncentbins][11];
    TH1D * qBAcnt2[ncentbins][11];
    TH1D * qCAcnt2[ncentbins][11];
    TH1D * qCBcnt2[ncentbins][11];
  } qanal[nanals];

  //===================================

  
  TH2D * Eff_0_5;
  TH2D * Eff_5_10;
  TH2D * Eff_10_30;
  TH2D * Eff_30_50;
  TH2D * Eff_50_100;
  enum    TrackCut {trackUndefine = 0, ppReco = 1, HIReco, Pixel, GenMC};
  TrackCut sTrackQuality;
  bool TrackQuality_ppReco(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);
  bool TrackQuality_HIReco(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);
  bool TrackQuality_Pixel(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);
  bool TrackQuality_GenMC(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);


  TRandom * ran;
#include "HeavyIonsAnalysis/VNAnalysis/interface/Harmonics.h"
  
  //===================================
  
  
  
  int fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup, int bin)
  {
    int Ntrk = 0;
    using namespace edm;
    using namespace reco;
    for(int i = 0; i<7; i++) {
      qxtrk[i]->Reset();
      qytrk[i]->Reset();
    }
    qcnt->Reset();
    avpt->Reset();
    iEvent.getByToken(vertexToken,vertex_);
    VertexCollection recoVertices = *vertex_;
    int primaryvtx = 0;
    if(sTrackQuality == ppReco) {
      if ( recoVertices.size() > 100 ) return -1;
      sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	  if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
	  return a.tracksSize() > b.tracksSize();
	});
    } else {
      sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	  return a.tracksSize() > b.tracksSize();
	});    
    }
          
    vtx = recoVertices[primaryvtx].z();
    hvtxRaw->Fill(vtx);

    if (fabs(vtx) < minvz_ || fabs(vtx) > maxvz_) {
      return -1;
    }        
    
    iEvent.getByLabel(trackTag_,trackCollection_);
    int k = (vtx-flatminvtx_)/flatdelvtx_;
    for(TrackCollection::const_iterator itTrack = trackCollection_->begin(); itTrack != trackCollection_->end(); ++itTrack) { 
      if ( sTrackQuality == HIReco and not TrackQuality_HIReco(itTrack, recoVertices) ) continue;
      else if ( sTrackQuality == ppReco and not TrackQuality_ppReco(itTrack, recoVertices) ) continue;
      else if ( sTrackQuality == Pixel  and not TrackQuality_Pixel (itTrack, recoVertices) ) continue;
      else if ( sTrackQuality == GenMC  and not TrackQuality_GenMC (itTrack, recoVertices) ) continue;
      double eff = 1.;
      if(k>=0&&k<flatnvtxbins_) {
	for(int iorder = 1; iorder <=7; iorder++) {
	  qxtrk[iorder-1]->Fill(itTrack->pt(), itTrack->eta(), eff*(TMath::Cos(iorder*itTrack->phi())));
	  qytrk[iorder-1]->Fill(itTrack->pt(), itTrack->eta(), eff*(TMath::Sin(iorder*itTrack->phi())));
	}
	qcnt->Fill(itTrack->pt(), itTrack->eta(), eff);
	avpt->Fill(itTrack->pt(), itTrack->eta(), eff*itTrack->pt());
      }
      if( itTrack->pt() < 0.2 ) continue;
      ++Ntrk;
    }
    return Ntrk;
  }
};


//
// constructors and destructor
//
VNAnalyzer::VNAnalyzer(const edm::ParameterSet& iConfig):runno_(0)  
{
  ran = new TRandom();
  ran->SetSeed(0);
  runno_ = 0;
  loadDB_ = kTRUE;
  FirstEvent_ = kTRUE;
  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    eporig[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    epmult[i] = 0;
    rescor[i] = 0;
    rescorErr[i] = 0;
  }
  
  centralityVariable_ = iConfig.getParameter<std::string>("centralityVariable");
  if(iConfig.exists("nonDefaultGlauberModel")){
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_+centralityMC_;
  
  centralityBinTag_ = iConfig.getParameter<edm::InputTag>("centralityBinTag_");
  centralityBinToken = consumes<int>(centralityBinTag_);
  
  bTag_ = iConfig.getParameter<edm::InputTag>("bTag_");
  bToken = consumes<double>(bTag_);
  
  RPTag_ = iConfig.getParameter<edm::InputTag>("RPTag_");
  RPToken = consumes<double>(RPTag_);
  
  centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);
  if(centralityToken.isUninitialized()) {
    std::cout<<"centralityToken is uninitialized."<<std::endl;
  }
  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);
  if(vertexToken.isUninitialized()) {
    std::cout<<"vertexToken is uninitialized."<<std::endl;
  }
  
  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);
  if(trackToken.isUninitialized()) {
    std::cout<<"trackToken is uninitialized."<<std::endl;
  }
  if ( trackTag_.label() == "hiGeneralTracks" ) {
    sTrackQuality = HIReco;
    cout<<"hiGeneralTracks"<<endl;
  } else if ( trackTag_.label() == "generalTracks" ) {
    sTrackQuality = ppReco;
    cout<<"generalTracks"<<endl;
  } else if ( trackTag_.label() == "hiGeneralAndPixelTracks" ) {
    sTrackQuality = Pixel;
    cout<<"hiGeneralAndPixelTracks"<<endl;
  } else if (trackTag_.label() == "mcEvtPlane") {
    sTrackQuality = GenMC;
    cout<<"mcEvtPlane"<<endl;
  } else {
    sTrackQuality = trackUndefine;
    std::cout<<"UNDEFINED: "<<trackTag_.label()<<std::endl;
  }
  
  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);
  if(inputPlanesToken.isUninitialized()) {
    std::cout<<"inputPlanesToken is uninitialized."<<std::endl;
  }
  tag_ = consumes<int>(iConfig.getParameter<edm::InputTag>("BinLabel"));
  
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",20);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  flatnvtxbins_ = iConfig.getParameter<int>("flatnvtxbins") ;
  flatminvtx_ = iConfig.getParameter<double>("flatminvtx") ;
  flatdelvtx_ = iConfig.getParameter<double>("flatdelvtx") ;
  
  effTable_ = iConfig.getParameter<std::string>("effTable_");
  makeTree_ = iConfig.getUntrackedParameter<bool>("makeTree_",false);
  
  nvtx_ = iConfig.getUntrackedParameter<int>("nvtx_", 100);
  dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 3.);
  d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
  pterror_ = iConfig.getUntrackedParameter<double>("pterror_",0.1);
  dzdzerror_pix_ = iConfig.getUntrackedParameter<double>("dzdzerror_pix_") ;
  chi2_  = iConfig.getUntrackedParameter<double>("chi2_") ;
  teff = 0;
  if(effTable_!="NULL") teff = new TrackEfficiency(effTable_.data());
  minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
  maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
  mx = ncentbins;
  
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"centralityBinTag_           "<<centralityBinTag_.encode()<<std::endl;
  std::cout<<"centralityTag_              "<<centralityTag_.encode()<<std::endl;
  std::cout<<"vertexTag_                  "<<vertexTag_.encode()<<std::endl;
  std::cout<<"trackTag_                   "<<trackTag_.encode()<<std::endl;
  std::cout<<"inputPlanesTag_             "<<inputPlanesTag_.encode()<<std::endl;
  std::cout<<"FlatOrder_                  "<<FlatOrder_<<std::endl;
  std::cout<<"NumFlatBins_                "<<NumFlatBins_<<std::endl;
  std::cout<<"flatnvtxbins_               "<<flatnvtxbins_<<std::endl;
  std::cout<<"flatminvtx_                 "<<flatminvtx_<<std::endl;
  std::cout<<"flatdelvtx_                 "<<flatdelvtx_<<std::endl;
  std::cout<<"CentBinCompression_         "<<CentBinCompression_<<std::endl;
  std::cout<<"effTable_                   "<<effTable_<<std::endl;
  std::cout<<"dzerror_                    "<<dzdzerror_<<endl;
  std::cout<<"d0error_                    "<<d0d0error_<<endl;
  std::cout<<"pterror_                    "<<pterror_<<endl;
  std::cout<<"nvtx_                       "<<nvtx_<<endl;
  std::cout<<"dzdzerror_pix_               "<<dzdzerror_pix_<<std::endl;
  std::cout<<"chi2_                        "<<chi2_<<std::endl;
  std::cout<<"==============================================="<<std::endl;
  TDirectory * save = gDirectory;
  TFileDirectory conddir = fs->mkdir("Conditions");
  conddir.make<TH1I>(centralityBinTag_.label().data(),centralityBinTag_.label().data(),1,0,1);
  conddir.make<TH1I>(centralityTag_.label().data(), centralityTag_.label().data(),1,0,1);
  conddir.make<TH1I>(vertexTag_.label().data(), vertexTag_.label().data(),1,0,1);
  conddir.make<TH1I>(trackTag_.label().data(), trackTag_.label().data(),1,0,1);
  conddir.make<TH1I>(inputPlanesTag_.label().data(), inputPlanesTag_.label().data(),1,0,1);
  string etable = Form("EffTable_%s",effTable_.data());
  conddir.make<TH1I>(etable.data(), etable.data(),1,0,1);
  string note_FlatOrder = Form("FlatOrder_%d",FlatOrder_);
  conddir.make<TH1I>(note_FlatOrder.data(), note_FlatOrder.data(),1,0,1);
  string note_NumFlatBins = Form("NumFlatBins_%d",NumFlatBins_);
  conddir.make<TH1I>(note_NumFlatBins.data(), note_NumFlatBins.data(),1,0,1);
  string note_flatnvtxbins = Form("flatnvtxbins_%d",flatnvtxbins_);
  conddir.make<TH1I>(note_flatnvtxbins.data(), note_flatnvtxbins.data(),1,0,1);
  string note_flatminvtx = Form("flatminvtx_%07.2f",flatminvtx_);
  conddir.make<TH1I>(note_flatminvtx.data(), note_flatminvtx.data(),1,0,1);
  string note_flatdelvtx = Form("flatdelvtx_%07.2f",flatdelvtx_);
  conddir.make<TH1I>(note_flatdelvtx.data(), note_flatdelvtx.data(),1,0,1);
  string note_dzdzerror = Form("dzdzerror_%07.2f",dzdzerror_);
  conddir.make<TH1I>(note_dzdzerror.data(), note_dzdzerror.data(),1,0,1);
  string note_d0d0error = Form("d0d0error_%07.2f",d0d0error_);
  conddir.make<TH1I>(note_d0d0error.data(), note_d0d0error.data(),1,0,1);
  string note_dzdzerror_pix = Form("dzdzerror_pix_%07.2f",dzdzerror_pix_);
  conddir.make<TH1I>(note_dzdzerror_pix.data(), note_dzdzerror_pix.data(),1,0,1);
  string note_chi2 = Form("chi2_%07.2f",chi2_);
  conddir.make<TH1I>(note_chi2.data(), note_chi2.data(),1,0,1);
  string note_vtx_range = Form("vtx_%5.1f_%5.1f",minvz_,maxvz_);
  conddir.make<TH1I>(note_vtx_range.data(), note_vtx_range.data(),1,0,1);
  string note_nvtx = Form("nvtx_%d",nvtx_);
  conddir.make<TH1I>(note_nvtx.data(), note_nvtx.data(),1,0,1);
  
  save->cd();
  hNtrk = fs->make<TH1D>("Ntrk","Ntrk",1001,0,3000);
  int npt = nptbins;
  for(int iorder = 1; iorder<=7; iorder++) {
    qxtrk[iorder-1] = fs->make<TH2D>(Form("qxtrk%d",iorder),Form("qxtrk%d",iorder),npt,ptbins, netabinsDefault, etabinsDefault);
    qytrk[iorder-1] = fs->make<TH2D>(Form("qytrk%d",iorder),Form("qytrk%d",iorder),npt,ptbins, netabinsDefault, etabinsDefault);
    qxtrk[iorder-1]->SetOption("colz");
    qytrk[iorder-1]->SetOption("colz");
    qxtrk[iorder-1]->Sumw2();
    qytrk[iorder-1]->Sumw2();
    qxtrk[iorder-1]->SetXTitle("p_{T} (GeV/c");
    qxtrk[iorder-1]->SetYTitle(Form("#eta (n=%d)",iorder));
    qytrk[iorder-1]->SetXTitle("p_{T} (GeV/c");
    qytrk[iorder-1]->SetYTitle(Form("#eta (n=%d)",iorder));
  }
  qcnt =  fs->make<TH2D>("qcnt", "qcnt",npt,ptbins, netabinsDefault, etabinsDefault);
  qcnt->SetXTitle("p_{T} (GeV/c");
  qcnt->SetYTitle("#eta");
  avpt =  fs->make<TH2D>("avpt","avpt",npt,ptbins, netabinsDefault, etabinsDefault);
  qcnt->SetOption("colz");
  avpt->SetOption("colz");
  qcnt->Sumw2();
  avpt->Sumw2();
  hTemplate = (TH2D *) qcnt->Clone("hTemplate");
  hTemplate->SetDirectory(0);
  hTemplate->Reset();
  hcent = fs->make<TH1D>("cent","cent",220,-10,110);
  hvtx = fs->make<TH1D>("vtx","vtx",600,-30,30);
  hvtxRaw = fs->make<TH1D>("vtxRaw","vtxRaw",600,-30,30);
  hb = fs->make<TH1D>("b","b",10000,0,10000);
  hRP = fs->make<TH1D>("rp","rp",1000,-4,4);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  hcentres = fs->make<TH1D>("centres","centres",ncentbins,centbins);
  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  NumFlatBins_ = ncentbins;
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,EPNames[i],EPOrder[i]);
    Double_t psirange = 4;
    if(EPOrder[i]==1 ) psirange = 3.5;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;
    if(EPOrder[i]==7) psirange = 0.5;
    
    hPsi[i] = subdir.make<TH1D>("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    hPsiOffset[i] = subdir.make<TH1D>("psiOffset","psiOffset",800,-psirange,psirange);
    hPsiOffset[i]->SetXTitle("#Psi");
    hPsiOffset[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    
    hPsiFlat[i] = subdir.make<TH1D>("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
  }
  //=====================
  TFileDirectory hardir = fs->mkdir("Harmonics");
  int nanalbins = ncentbins;
  //==============   Resolution terms  ========
  TFileDirectory resdir = fs->mkdir("Resolutions");
  for(int i = 0; i<nanalbins; i++) {
    TFileDirectory ressubdir;
    ressubdir = resdir.mkdir(Form("%d_%d",(int)centbins[i],(int)centbins[i+1]));
    for(int iorder = 1; iorder<=7; iorder++) {
      res[iorder-1][i] = ressubdir.make<TH2D>(Form("res%d",iorder),Form("res%d",iorder),50,0,50,50,0,50);
      resw[iorder-1][i] = ressubdir.make<TH2D>(Form("resw%d",iorder),Form("resw%d",iorder),50,0,50,50,0,50);
      resep[iorder-1][i] = ressubdir.make<TH2D>(Form("resep%d",iorder),Form("resep%d",iorder),50,0,50,50,0,50);
      rescnt[iorder-1][i] = ressubdir.make<TH2D>(Form("rescnt%d",iorder),Form("rescnt%d",iorder),50,0,50,50,0,50);
      res[iorder-1][i]->Reset();
      res[iorder-1][i]->Sumw2();
      res[iorder-1][i]->SetOption("colz");
      resw[iorder-1][i]->Reset();
      resw[iorder-1][i]->Sumw2();
      resw[iorder-1][i]->SetOption("colz");
      resep[iorder-1][i]->Reset();
      resep[iorder-1][i]->Sumw2();
      resep[iorder-1][i]->SetOption("colz");
      rescnt[iorder-1][i]->Reset();
      rescnt[iorder-1][i]->Sumw2();
      rescnt[iorder-1][i]->SetOption("colz");
    }
  }
  
  for(int i = 0; i<nanalbins; i++) {

    TFileDirectory subdir;
    subdir = hardir.mkdir(Form("%d_%d",(int)centbins[i],(int)centbins[i+1]));
    ptav[i] = subdir.make<TH2D>("ptav","ptav",npt,ptbins, netabinsDefault, etabinsDefault);
    ptcnt[i] = subdir.make<TH2D>("ptcnt","ptcnt",npt,ptbins, netabinsDefault, etabinsDefault);
    ptav[i]->Sumw2();
    ptcnt[i]->Sumw2();

    ptav[i]->SetOption("colz");
    ptcnt[i]->SetOption("colz");
    for(int ian = 0; ian<nanals; ian++) {
      TFileDirectory andir = subdir.mkdir(AnalNames[ian].data());
      
      qanal[ian].qA[i][0] = andir.make<TH2D>("qA","qA",npt,ptbins, netabinsDefault, etabinsDefault);
      qanal[ian].qB[i][0] = andir.make<TH2D>("qB","qB",npt,ptbins, netabinsDefault, etabinsDefault);
      qanal[ian].wnA[i][0] = andir.make<TH2D>("wnA","wnA",npt,ptbins, netabinsDefault, etabinsDefault);
      qanal[ian].wnB[i][0] = andir.make<TH2D>("wnB","wnB",npt,ptbins, netabinsDefault, etabinsDefault);
      qanal[ian].qA[i][0]->Sumw2();
      qanal[ian].qB[i][0]->Sumw2();
      qanal[ian].wnA[i][0]->Sumw2();
      qanal[ian].wnB[i][0]->Sumw2();
      qanal[ian].qA[i][0]->SetOption("colz");
      qanal[ian].qB[i][0]->SetOption("colz");
      qanal[ian].wnA[i][0]->SetOption("colz");
      qanal[ian].wnB[i][0]->SetOption("colz");
      
      qanal[ian].thA[i][0]  = andir.make<TH1D>("thA","thA",160,-4,4);
      qanal[ian].thB[i][0]  = andir.make<TH1D>("thB","thB",160,-4,4);
      qanal[ian].thC[i][0]  = andir.make<TH1D>("thC","thC",160,-4,4);
      qanal[ian].thA2[i][0] = andir.make<TH1D>("thA2","thA2",160,-4,4);
      qanal[ian].thB2[i][0] = andir.make<TH1D>("thB2","thB2",160,-4,4);
      qanal[ian].thC2[i][0] = andir.make<TH1D>("thC2","thC2",160,-4,4);
      qanal[ian].thn[i][0]  = andir.make<TH1D>("thn","thn",160,-4,4);

      qanal[ian].qBA[i][0] = andir.make<TH1D>("qBA","qBA",1,0,1);
      qanal[ian].qCA[i][0] = andir.make<TH1D>("qCA","qCA",1,0,1);
      qanal[ian].qCB[i][0] = andir.make<TH1D>("qCB","qCB",1,0,1);
      qanal[ian].qBAcnt[i][0] = andir.make<TH1D>("qBAcnt","qBAcnt",1,0,1);
      qanal[ian].qCAcnt[i][0] = andir.make<TH1D>("qCAcnt","qCAcnt",1,0,1);
      qanal[ian].qCBcnt[i][0] = andir.make<TH1D>("qCBcnt","qCBcnt",1,0,1);
 

      qanal[ian].qBA2[i][0] = andir.make<TH1D>("qBA2","qBA2",1,0,1);
      qanal[ian].qCA2[i][0] = andir.make<TH1D>("qCA2","qCA2",1,0,1);
      qanal[ian].qCB2[i][0] = andir.make<TH1D>("qCB2","qCB2",1,0,1);
      qanal[ian].qBAcnt2[i][0] = andir.make<TH1D>("qBAcnt2","qBAcnt2",1,0,1);
      qanal[ian].qCAcnt2[i][0] = andir.make<TH1D>("qCAcnt2","qCAcnt2",1,0,1);
      qanal[ian].qCBcnt2[i][0] = andir.make<TH1D>("qCBcnt2","qCBcnt2",1,0,1);
     
      TFileDirectory subev = andir.mkdir("SubEvents");
      for(int j = 1; j<=10; j++) {	  
	qanal[ian].qA[i][j] = subev.make<TH2D>(Form("qA_%d",j),Form("qA_%d",j),npt,ptbins, netabinsDefault, etabinsDefault); 
	qanal[ian].qB[i][j] = subev.make<TH2D>(Form("qB_%d",j),Form("qB_%d",j),npt,ptbins, netabinsDefault, etabinsDefault);
	qanal[ian].wnA[i][j] = subev.make<TH2D>(Form("wnA_%d",j),Form("wnA_%d",j),npt,ptbins, netabinsDefault, etabinsDefault);
	qanal[ian].wnB[i][j] = subev.make<TH2D>(Form("wnB_%d",j),Form("wnB_%d",j),npt,ptbins, netabinsDefault, etabinsDefault);
	qanal[ian].qA[i][j]->Sumw2();
	qanal[ian].qB[i][j]->Sumw2();
	qanal[ian].wnA[i][j]->Sumw2();
	qanal[ian].wnB[i][j]->Sumw2();
	qanal[ian].qA[i][j]->SetOption("colz");
	qanal[ian].qB[i][j]->SetOption("colz");
	qanal[ian].wnA[i][j]->SetOption("colz");
	qanal[ian].wnB[i][j]->SetOption("colz");

	qanal[ian].thA[i][j]  = subev.make<TH1D>(Form("thA_%d",j),Form("thA_%d",j),160,-4,4);
	qanal[ian].thB[i][j]  = subev.make<TH1D>(Form("thB_%d",j),Form("thB_%d",j),160,-4,4);
	qanal[ian].thC[i][j]  = subev.make<TH1D>(Form("thC_%d",j),Form("thC_%d",j),160,-4,4);
	qanal[ian].thA2[i][j]  = subev.make<TH1D>(Form("thA2_%d",j),Form("thA2_%d",j),160,-4,4);
	qanal[ian].thB2[i][j]  = subev.make<TH1D>(Form("thB2_%d",j),Form("thB2_%d",j),160,-4,4);
	qanal[ian].thC2[i][j]  = subev.make<TH1D>(Form("thC2_%d",j),Form("thC2_%d",j),160,-4,4);
	qanal[ian].thn[i][j]  = subev.make<TH1D>(Form("thn_%d",j),Form("thn_%d",j),160,-4,4);
	
	qanal[ian].qBA[i][j] = subev.make<TH1D>(Form("qBA_%d",j),Form("qBA_%d",j),1,0,1);
	qanal[ian].qCA[i][j] = subev.make<TH1D>(Form("qCA_%d",j),Form("qCA_%d",j),1,0,1);      
	qanal[ian].qCB[i][j] = subev.make<TH1D>(Form("qCB_%d",j),Form("qCB_%d",j),1,0,1);
	
	qanal[ian].qBAcnt[i][j] = subev.make<TH1D>(Form("qBAcnt_%d",j),Form("qBAcnt_%d",j),1,0,1);
	qanal[ian].qCAcnt[i][j] = subev.make<TH1D>(Form("qCAcnt_%d",j),Form("qCAcnt_%d",j),1,0,1);      
	qanal[ian].qCBcnt[i][j] = subev.make<TH1D>(Form("qCBcnt_%d",j),Form("qCBcnt_%d",j),1,0,1);
	
	qanal[ian].qBA2[i][j] = subev.make<TH1D>(Form("qBA2_%d",j),Form("qBA2_%d",j),1,0,1);
	qanal[ian].qCA2[i][j] = subev.make<TH1D>(Form("qCA2_%d",j),Form("qCA2_%d",j),1,0,1);      
	qanal[ian].qCB2[i][j] = subev.make<TH1D>(Form("qCB2_%d",j),Form("qCB2_%d",j),1,0,1);
	
	qanal[ian].qBAcnt2[i][j] = subev.make<TH1D>(Form("qBAcnt2_%d",j),Form("qBAcnt2_%d",j),1,0,1);
	qanal[ian].qCAcnt2[i][j] = subev.make<TH1D>(Form("qCAcnt2_%d",j),Form("qCAcnt2_%d",j),1,0,1);      
	qanal[ian].qCBcnt2[i][j] = subev.make<TH1D>(Form("qCBcnt2_%d",j),Form("qCBcnt2_%d",j),1,0,1);
	
      }
      
    }
  }
  //==============================
  if(makeTree_) {
    tree = fs->make<TTree>("tree","EP tree");
    tree->Branch("Cent",&centval,"cent/D");
    tree->Branch("ntrkflat",&ntrkval,"nofftrak/I");
    tree->Branch("Vtx",&vtx,"vtx/D");
    tree->Branch("epang",&epang, epnames.Data());
    tree->Branch("eporig",&eporig, epnames.Data());
    tree->Branch("qx",      &qx,       epnames.Data());
    tree->Branch("qy",      &qy,       epnames.Data());
    tree->Branch("q",       &q,       epnames.Data());
    tree->Branch("vn", &vn, epnames.Data());
    tree->Branch("mult",    &epmult,  epnames.Data());
    tree->Branch("sumw",    &sumw,  epnames.Data());
    tree->Branch("sumw2",    &sumw2,  epnames.Data());
    tree->Branch("Run",     &runno_,   "run/i");
    tree->Branch("Rescor",  &rescor,   epnames.Data());
    tree->Branch("RescorErr",  &rescorErr,   epnames.Data());
    tree->Branch("qcnt",    "TH2D",  &qcnt, 128000, 0);
    tree->Branch("avpt",    "TH2D",  &avpt, 128000, 0);
  }
}


VNAnalyzer::~VNAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
    //
    //Get Size of Centrality Table
    //
    nCentBins_ = 200;
    
    edm::ESHandle<CentralityTable> centDB_;
    iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
    nCentBins_ = (int) centDB_->m_table.size();
    for(int i = 0; i<NumEPNames; i++) {
      flat[i]->setCaloCentRefBins(-1,-1);
    }
    
    //
    // Get flattening parameter file.  
    //
    
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
    if(!db->IsSuccess()) {
      std::cout<<"Flattening db inconsistancy, will continue without: "<<std::endl;
      loadDB_ = kFALSE;
    }        
    
  } //First event
  
  //
  // Get Centrality
  //

  int cbin = 0;
  iEvent.getByToken(centralityBinToken, cbin_);
  cbin = *cbin_;
  double cscale = 100./nCentBins_;
  centval = cscale*cbin;
  hcentres->Fill(centval);
  
  int ibin = hcentres->FindBin(centval)-1;
  hcent->Fill(centval);
  hcentbins->Fill(cbin);
  int ntrkval=fillTracks(iEvent, iSetup, ibin);
  if(ntrkval<=0) return;
  hvtx->Fill(vtx);
  hNtrk->Fill(ntrkval);
  //
  // Get Event Planes
  //
  iEvent.getByToken(inputPlanesToken,inputPlanes_);
  
  if(!inputPlanes_.isValid()){
    cout << "Error! Can't get hiEvtPlaneFlat product!" << endl;
    return ;
  }
  
  Int_t indx = 0;
  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    vn[i] = 0;
    epmult[i] = 0;
    sumw[i] = 0;
    sumw2[i] = 0;
  }
  for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
    if(indx != rp->indx() ) {
      cout<<"EP inconsistency found. Abort."<<endl;
      return;
    }
    if(rp->sumSin()!=0 || rp->sumCos()!=0) {
      if(rp->mult()>3 && fabs(vtx)<15) {
	epang[indx]=rp->angle();
	eporig[indx]=rp->angle(0);
	epsin[indx] = rp->sumSin();
	epcos[indx] = rp->sumCos();
	hPsi[indx]->Fill(rp->angle(0));
	hPsiOffset[indx]->Fill(rp->angle(1));
	hPsiFlat[indx]->Fill(rp->angle(2));
	qx[indx]  = rp->qx(); 
	qy[indx]  = rp->qy();
	q[indx]   = rp->q();
	vn[indx] = rp->vn(0);
	epmult[indx] = (double) rp->mult();
	sumw[indx] = rp->sumw();
	sumw2[indx] = rp->sumw2();
	rescor[indx] = flat[indx]->getCentRes1((int) centval);
	rescorErr[indx] = flat[indx]->getCentResErr1((int) centval);
      } 
    }
    ++indx; 
  }
  ptav[ibin]->Add(avpt);
  ptcnt[ibin]->Add(qcnt);
  for(int ian = 0; ian<nanals; ian++) {
     if(ian==N2) Fill_N( N2, 2., ibin, qxtrk[1], qytrk[1], qcnt, qx[HFp2], qy[HFp2], qx[HFm2], qy[HFm2], qx[trackmid2], qy[trackmid2], sumw[HFp2], sumw[HFm2], sumw[trackmid2]);
     if(ian==N3) Fill_N( N3, 3., ibin, qxtrk[2], qytrk[2], qcnt, qx[HFp3], qy[HFp3], qx[HFm3], qy[HFm3], qx[trackmid3], qy[trackmid3], sumw[HFp3], sumw[HFm3], sumw[trackmid3]);
     if(ian==N4) Fill_N( N4, 4., ibin, qxtrk[3], qytrk[3], qcnt, qx[HFp4], qy[HFp4], qx[HFm4], qy[HFm4], qx[trackmid4], qy[trackmid4], sumw[HFp4], sumw[HFm4], sumw[trackmid4]);
 }
  for(int iorder = 1; iorder <=4; iorder++) {
    int epmin = 0;
    int epmax = 0;
    if(iorder == 1 ) {
      epmin = HFm1;
      epmax = Castor1;
    }else if (iorder == 2) {
      epmin = HFm2;
      epmax = Castor2;
    }else if (iorder == 3) {
      epmin = HFm3;
      epmax = trackp3;
    }else if (iorder == 4) {
      epmin = HFm4;
      epmax = trackp4;
    }
    for(int i = epmin; i<= epmax; i++) {
      for(int j = i; j<=epmax; j++) {
	double ii = i-epmin;
	double jj = j-epmin;
	double w =(pow(qx[i],2)+pow(qy[i],2))*(pow(qx[j],2)+pow(qy[j],2));
	if(w>0) {
	  res[iorder-1][ibin]->Fill(ii,jj,   qx[i]*qx[j]+qy[i]*qy[j]);
	  resw[iorder-1][ibin]->Fill(ii,jj,  sumw[i]*sumw[j]);
	  resep[iorder-1][ibin]->Fill(ii,jj, (qx[i]*qx[j]+qy[i]*qy[j])/sqrt(w));
	  rescnt[iorder-1][ibin]->Fill(ii,jj);
	}
      }
    }
    
  }
  
  if(makeTree_) tree->Fill(); 
}



// ------------ method called once each job just before starting event loop  ------------
void 
VNAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VNAnalyzer::endJob() {
}
///
bool
VNAnalyzer::TrackQuality_ppReco(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
  if ( itTrack->charge() == 0 ) return false;
  if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
  if ( itTrack->ptError()/itTrack->pt() > pterror_ ) return false;
  int primaryvtx = 0;
  math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
  double vxError = recoVertices[primaryvtx].xError();
  double vyError = recoVertices[primaryvtx].yError();
  double vzError = recoVertices[primaryvtx].zError();
  double d0 = -1.* itTrack->dxy(v1);
  double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
  if ( fabs( d0/derror ) > d0d0error_ ) {
    return false;
  }
  double dz=itTrack->dz(v1);
  double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
  if ( fabs( dz/dzerror ) > dzdzerror_ ) {
    return false;
  }
  return true;
}

///
bool
VNAnalyzer::TrackQuality_HIReco(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
  if ( itTrack->charge() == 0 ) return false;
  if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
  if ( itTrack->numberOfValidHits() < 11 ) return false;
  if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
    return false;
  }
  if ( itTrack->ptError()/itTrack->pt() > pterror_ ) {
    return false;
  }
  if (
      itTrack->originalAlgo() != 4 and
      itTrack->originalAlgo() != 5 and
      itTrack->originalAlgo() != 6 and
      itTrack->originalAlgo() != 7
      ) {
    return false;
  }
  
  int primaryvtx = 0;
  math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
  double vxError = recoVertices[primaryvtx].xError();
  double vyError = recoVertices[primaryvtx].yError();
  double vzError = recoVertices[primaryvtx].zError();
  double d0 = -1.* itTrack->dxy(v1);
  double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
  if ( fabs( d0/derror ) > d0d0error_ ) {
    return false;
  }
  
  double dz=itTrack->dz(v1);
  double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
  if ( fabs( dz/dzerror ) > dzdzerror_ ) {
    return false;
  }
  return true;
}

///
bool
VNAnalyzer::TrackQuality_Pixel(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
  if ( itTrack->charge() == 0 ) return false;
  if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
  bool bPix = false;
  int nHits = itTrack->numberOfValidHits();
  
  int primaryvtx = 0;
  math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
  double vxError = recoVertices[primaryvtx].xError();
  double vyError = recoVertices[primaryvtx].yError();
  double vzError = recoVertices[primaryvtx].zError();
  double d0 = -1.* itTrack->dxy(v1);
  
  double dz=itTrack->dz(v1);
  double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
  if ( itTrack->pt() < 2.4 and (nHits==3 or nHits==4 or nHits==5 or nHits==6) ) bPix = true;
  if ( not bPix ) {
    if ( nHits < 11 ) return false;
    if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
      return false;
    }
    if ( itTrack->ptError()/itTrack->pt() > pterror_ ) {
      return false;
    }
    if (
	itTrack->pt() > 2.4 and
	itTrack->originalAlgo() != 4 and
	itTrack->originalAlgo() != 5 and
	itTrack->originalAlgo() != 6 and
	itTrack->originalAlgo() != 7
	) {
      return false;
    }
    
    double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
    if ( fabs( d0/derror ) > d0d0error_ ) {
      return false;
    }
    
    if ( fabs( dz/dzerror ) > dzdzerror_ ) {
      return false;
    }
  } else {
    if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > chi2_ ) return false;
    if ( fabs( dz/dzerror ) > dzdzerror_pix_ ) {
      return false;
    }
  }
  return true;
}
///
bool
VNAnalyzer::TrackQuality_GenMC(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
	if ( itTrack->charge() == 0 ) return false;
	if( fabs(itTrack->eta())>2.4  )return false;
	return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(VNAnalyzer);

