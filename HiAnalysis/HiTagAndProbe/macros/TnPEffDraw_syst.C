#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include "TSystem.h"
#include "TTree.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TPave.h"
#include "TText.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFitResult.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"
#include "TMath.h"

#include "toyStudy.h"

using namespace RooFit;
using namespace std;

/////////////////////////////////////////////////
//              P A R A M E T E R S            //
/////////////////////////////////////////////////

// Choose the efficiency type.
// Possible values: MUIDTRG, MUID, TRG, TRK, STA
#define MUIDTRG

// pp or PbPb?
bool isPbPb = false; // if true, will compute the centrality dependence
TString collTag = "PbPb"; // isPbPb ? "PbPb" : "pp";

// do the toy study for the correction factors? (only applies if MUIDTRG)
bool doToys = false;

// Location of the files
const int nCentBins = 7;
const char* fCentMCNames[nCentBins] = 
{
   "fitmc/MuIDTrg_cent_20150909/res_Cent0005_/MuIDTrg_cent_20150909_Cent0005_.root",
   "fitmc/MuIDTrg_cent_20150909/res_Cent0510_/MuIDTrg_cent_20150909_Cent0510_.root",
   "fitmc/MuIDTrg_cent_20150909/res_Cent1020_/MuIDTrg_cent_20150909_Cent1020_.root",
   "fitmc/MuIDTrg_cent_20150909/res_Cent2040_/MuIDTrg_cent_20150909_Cent2040_.root",
   "fitmc/MuIDTrg_cent_20150909/res_Cent4050_/MuIDTrg_cent_20150909_Cent4050_.root",
   "fitmc/MuIDTrg_cent_20150909/res_Cent5060_/MuIDTrg_cent_20150909_Cent5060_.root",
   "fitmc/MuIDTrg_cent_20150909/res_Cent60100_/MuIDTrg_cent_20150909_Cent60100_.root",
};
const char* fCentDataNames[nCentBins] = 
{
   "fitdata/MuIdTrg_cent_20150911/res_Cent0005_/muidtrg_Cent0005_.root",
   "fitdata/MuIdTrg_cent_20150911/res_Cent0510_/muidtrg_Cent0510_.root",
   "fitdata/MuIdTrg_cent_20150911/res_Cent1020_/muidtrg_Cent1020_.root",
   "fitdata/MuIdTrg_cent_20150911/res_Cent2040_/muidtrg_Cent2040_.root",
   "fitdata/MuIdTrg_cent_20150911/res_Cent4050_/muidtrg_Cent4050_.root",
   "fitdata/MuIdTrg_cent_20150911/res_Cent5060_/muidtrg_Cent5060_.root",
   "fitdata/MuIdTrg_cent_20150911/res_Cent60100_/muidtrg_Cent60100_.root",
};

const int nSyst = 5;
// the first file is for the nominal case, the following ones are for the systematics
const char* fDataName[nSyst] = {
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/nominal/pbpb_data/tnp_Ana_RD_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/massrange/pbpb_data/tnp_Ana_RD_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/signalfcn/pbpb_data/tnp_Ana_RD_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/bkgdfcn/pbpb_data/tnp_Ana_RD_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/tagsel/pbpb_data/tnp_Ana_RD_PbPb_MuonIDTrg_AllMB.root",
};
const char* fMCName[nSyst] = {
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/nominal/pbpb_mc/tnp_Ana_MC_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/massrange/pbpb_mc/tnp_Ana_MC_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/signalfcn/pbpb_mc/tnp_Ana_MC_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/bkgdfcn/pbpb_mc/tnp_Ana_MC_PbPb_MuonIDTrg_AllMB.root",
   "/home/emilien/Documents/Postdoc_LLR/TagAndProbe/systs/tagsel/pbpb_mc/tnp_Ana_MC_PbPb_MuonIDTrg_AllMB.root",
};

// names for systematics
const char* systName[nSyst] = {
   "nominal",
   "massrange",
   "signalfcn",
   "bkgfcn",
   "tagsel"
};



//////////////////////////////////////////////////////////////////////////

// Other parameters
const double tx[nCentBins]={381.23,330.35,261.37,158.7,86.3,53.5,16.2};
#ifdef STA
TString etaTag("PassingSTA_eta");
TString ptTag("PassingSTA_pt");
TString allTag("PassingSTA_1bin");
#else
TString etaTag("PassingGlb_eta");
TString ptTag("PassingGlb_pt");
TString allTag("PassingGlb_1bin");
#endif // ifdef STA

// Automatic parameters (shouldn't require modification by the user)
#if defined MUIDTRG || defined STA
TString absetaVar("abseta");
ofstream file_sfs("correction_functions.txt");
#ifdef MUIDTRG
const int nAbsEtaBins = 4;
TString absetaTag("PassingGlb_abseta");
TString cutTag("MuonIDTrg"); 
TString cutLegend("Muon ID + trigger");
#else
const int nAbsEtaBins = 4;//2;
TString absetaTag("PassingSTA_abseta");
TString cutTag = isPbPb ? TString("fitTkFromSta") : TString("tpTree"); 
TString cutLegend("STA muon");
#endif // ifdef MUIDTRG
#else
TString absetaVar("eta");
TString absetaTag("PassingGlb_1bin");
const int nAbsEtaBins = 1;
#ifdef TRG
TString cutTag("MuonTrg"); 
TString cutLegend("Muon trigger");
#else
#ifdef TRK
TString cutTag("MuonTrk"); 
TString cutLegend("Inner tracking");
#else
TString cutTag("MuonIDTrg"); 
TString cutLegend("Muon ID");
#endif // ifdef TRK
#endif // ifdef TRG
#endif // ifdef MUIDTRG

#ifdef TRK
const double effmin = 0.8;
const double sfrange = 0.08;
#else
#ifdef MUID
const double effmin = 0.8;
const double sfrange = 0.08;
#else
const double effmin = 0.;
const double sfrange = 0.55;
#endif // ifdef MUID
#endif // ifdef TRK

// Function Define
TH2F *plotEff2D(RooDataSet *a, TString b);
vector<TGraphAsymmErrors*> plotEff_Nbins(RooDataSet *a, int aa, const char* varx, const char* var2);
TGraphAsymmErrors *plotEff_1bin(RooDataSet *a, int aa, const char* varx, int rebin=1, double* tntot=NULL);
TGraphAsymmErrors *plotEffCent(RooDataSet **a1,int aa);
void formatTH1F(TH1* a, int b, int c, int d);
void formatTGraph(TGraph* a, int b, int c, int d);
void formatTLeg(TLegend* a);
void CalEffErr(TGraph *a, double *b);
void CalEffErr(vector<TGraphAsymmErrors*> a, double **b);
void plotSysts(TGraphAsymmErrors *graphs[nSyst], TCanvas *c1, TPad *p1, TH1F *h1, TPad *pr, TH1F *hr, TString header, TString tag);

// From here you need to set up your environments.
void TnPEffDraw_syst() {

  // gROOT->Macro("~/logon.C");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleYOffset(1.0);


  //centrallity root files
  TFile* fCentMC[nCentBins];
  TFile* fCentData[nCentBins];

  if (isPbPb) {
     for (int i=0; i<nCentBins; i++)
     {
        fCentMC[i] = new TFile(fCentMCNames[i]);
        fCentData[i] = new TFile(fCentDataNames[i]);
     }
  }
  
  //data and MC root files as well as single bin for integrated efficiency
  TFile* f9[nSyst];
  TFile* f10[nSyst];

  for (int i=0; i<nSyst; i++) {
     f9[i] = new TFile(fMCName[i]);
     f10[i] = new TFile(fDataName[i]);
  }

  TCanvas *c1 = new TCanvas("c1","",700,600);

  vector<RooDataSet*> daPtData0[nSyst], daPtData1[nSyst];

  for (int k=0; k<nSyst; k++) {
     if (nAbsEtaBins==1)
     {
        daPtData0[k].push_back((RooDataSet*)f9[k]->Get(cutTag + "/" + ptTag + "/fit_eff")); 
        daPtData1[k].push_back((RooDataSet*)f10[k]->Get(cutTag + "/" + ptTag + "/fit_eff")); 
     }
     else
     {
        for (int i=0; i<nAbsEtaBins; i++)
        {
           daPtData0[k].push_back((RooDataSet*)f9[k]->Get(cutTag + "/" + ptTag + Form("%i/fit_eff",i+1))); 
           daPtData1[k].push_back((RooDataSet*)f10[k]->Get(cutTag + "/" + ptTag + Form("%i/fit_eff",i+1))); 
        }
     }
  }

  vector<TGraphAsymmErrors*> ComPt0[nSyst], ComPt1[nSyst];
  
  for (int k=0; k<nSyst; k++) {
     for (unsigned int i=0; i<daPtData0[k].size(); i++)
     {
        ComPt0[k].push_back(plotEff_1bin(daPtData0[k][i], 1, "pt"));
        ComPt1[k].push_back(plotEff_1bin(daPtData1[k][i], 1, "pt"));
     }
  }


  RooDataSet* daEtaData0[nSyst]; 
  RooDataSet* daEtaData1[nSyst];

  for (int i=0; i<nSyst; i++) {
     daEtaData0[i] = (RooDataSet*)f9[i]->Get(cutTag + "/" + etaTag + "/fit_eff"); 
     daEtaData1[i] = (RooDataSet*)f10[i]->Get(cutTag + "/" + etaTag + "/fit_eff"); 
  }

  TGraphAsymmErrors* ComEta0[nSyst];
  TGraphAsymmErrors* ComEta1[nSyst];

  for (int i=0; i<nSyst; i++) {
     ComEta0[i] = plotEff_1bin(daEtaData0[i], 1, "eta");
     ComEta1[i] = plotEff_1bin(daEtaData1[i], 1, "eta");
  }

  RooDataSet *daCentMC1Bin[nCentBins];
  RooDataSet *daCentData1Bin[nCentBins];
  if (isPbPb) {
     for (int i=0; i<nCentBins; i++)
     {
        daCentMC1Bin[i] = (RooDataSet*) fCentMC[i]->Get(cutTag + "/" + allTag + "/fit_eff");
        daCentData1Bin[i] = (RooDataSet*) fCentData[i]->Get(cutTag + "/" + allTag + "/fit_eff");
     }
  }

  RooDataSet* daPtMC1Bin0[nSyst];
  RooDataSet* daPtData1Bin0[nSyst];
  RooDataSet* daAbsEtaMC1[nSyst];
  RooDataSet* daAbsEtaData1[nSyst];

  for (int i=0; i<nSyst; i++) {
     daPtMC1Bin0[i] = (RooDataSet*)f9[i]->Get(cutTag + "/" + allTag + "/fit_eff"); 
     daPtData1Bin0[i] = (RooDataSet*)f10[i]->Get(cutTag + "/" + allTag + "/fit_eff"); 
     daAbsEtaMC1[i] = (RooDataSet*)f9[i]->Get(cutTag + "/" + absetaTag + "/fit_eff"); 
     daAbsEtaData1[i] = (RooDataSet*)f10[i]->Get(cutTag + "/" + absetaTag + "/fit_eff"); 
  }

  TGraphAsymmErrors* Com0Pt0[nSyst];
  TGraphAsymmErrors* Com0Pt1[nSyst];
  vector<TGraphAsymmErrors*> Com0AbsEta0[nSyst];
  vector<TGraphAsymmErrors*> Com0AbsEta1[nSyst];

  for (int k=0; k<nSyst; k++) {
     Com0Pt0[k] = plotEff_1bin(daPtMC1Bin0[k],0,"eta");
     Com0Pt1[k] = plotEff_1bin(daPtData1Bin0[k],0,"eta");
     Com0AbsEta0[k] = plotEff_Nbins(daAbsEtaMC1[k],0,"pt",absetaVar);
     Com0AbsEta1[k] = plotEff_Nbins(daAbsEtaData1[k],0,"pt",absetaVar);
  }

  TGraphAsymmErrors *effCentMC = isPbPb ? plotEffCent(daCentMC1Bin, 1) : NULL;
  TGraphAsymmErrors *effCentData = isPbPb ? plotEffCent(daCentData1Bin, 1) : NULL;

  if (isPbPb) {
     effCentMC->SetMarkerStyle(20);
     effCentMC->SetMarkerSize(1.4);
     effCentMC->SetMarkerColor(kRed+1);
     effCentMC->SetLineColor(kRed+1);
     effCentData->SetMarkerStyle(25);
     effCentData->SetMarkerSize(1.4);
     effCentData->SetMarkerColor(kBlue+1);
     effCentData->SetLineColor(kBlue+1);
  }

  int nbins_abseta = ComPt0[0].size();
  for (int k=0; k<nSyst; k++)
  {
     for (int i=0; i<nbins_abseta; i++)
     {
        ComPt0[k][i]->SetMarkerStyle(20);
        ComPt0[k][i]->SetMarkerSize(1.4);
        ComPt0[k][i]->SetMarkerColor(kRed+1);
        ComPt0[k][i]->SetLineColor(kRed+1);
        ComPt1[k][i]->SetMarkerStyle(25);
        ComPt1[k][i]->SetMarkerSize(1.4);
        ComPt1[k][i]->SetMarkerColor(kBlue+1);
        ComPt1[k][i]->SetLineColor(kBlue+1);
     }

     ComEta0[k]->SetMarkerStyle(20);
     ComEta0[k]->SetMarkerSize(1.4);
     ComEta0[k]->SetMarkerColor(kRed+1);
     ComEta0[k]->SetLineColor(kRed+1);
     ComEta1[k]->SetMarkerStyle(25);
     ComEta1[k]->SetMarkerSize(1.4);
     ComEta1[k]->SetMarkerColor(kBlue+1);
     ComEta1[k]->SetLineColor(kBlue+1);
  }


  // 1 Bin
  double Trk0[nSyst][4];
  double Trk1[nSyst][4];
  double*** TrkAbsEta0 = new double**[nSyst];
  double*** TrkAbsEta1 = new double**[nSyst];
  for (int k=0; k<nSyst; k++) {
     TrkAbsEta0[k] = new double*[nbins_abseta];
     TrkAbsEta1[k] = new double*[nbins_abseta];
     for (int i=0; i<nbins_abseta; i++) 
     {
        TrkAbsEta0[k][i] = new double[4];
        TrkAbsEta1[k][i] = new double[4];
     }
     CalEffErr(Com0Pt0[k], Trk0[k]);
     CalEffErr(Com0Pt1[k], Trk1[k]);
     CalEffErr(Com0AbsEta0[k], TrkAbsEta0[k]);
     CalEffErr(Com0AbsEta1[k], TrkAbsEta1[k]);
  }


  c1->cd();
  TPad *pad1 = new TPad("pad1","pad1",0, 0.26 ,1,1);
  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3*0.96);
  pad1->SetBottomMargin(0.04);
  pad2->SetTopMargin(0);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  pad2->SetBottomMargin(gStyle->GetPadBottomMargin()/0.3);
  pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
  // pad2->SetGridy();
  pad1->Draw();
  pad1->cd();

  // // adapt text size
  // lTextSize *= 1./0.7;


  TH1F *hPad = new TH1F("hPad",";p^{#mu}_{T} [GeV/c];Single #mu Efficiency",5,0,30);
  TH1F *hPad1 = new TH1F("hPad1",";#eta^{#mu} ;Single #mu Efficiency",5,-2.4,2.4);
  TH1F *hPad2 = new TH1F("hPad2",";N_{part} ;Single #mu Efficiency",5,0,400);
  hPad->GetXaxis()->CenterTitle();
  hPad1->GetXaxis()->CenterTitle();
  hPad2->GetXaxis()->CenterTitle();
  hPad->GetXaxis()->SetLabelSize(0.);
  hPad->GetXaxis()->SetTitleSize(0.);
  hPad->GetXaxis()->SetTitleOffset(1.2);
  hPad->GetYaxis()->SetLabelSize(0.05);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetTitleOffset(1.);
  hPad1->GetXaxis()->SetLabelSize(0.);
  hPad1->GetXaxis()->SetTitleSize(0.);
  hPad1->GetXaxis()->SetTitleOffset(1.2);
  hPad1->GetYaxis()->SetLabelSize(0.05);
  hPad1->GetYaxis()->SetTitleSize(0.05);
  hPad1->GetYaxis()->SetTitleOffset(1.);
  hPad2->GetXaxis()->SetLabelSize(0.);
  hPad2->GetXaxis()->SetTitleSize(0.);
  hPad2->GetXaxis()->SetTitleOffset(1.2);
  hPad2->GetYaxis()->SetLabelSize(0.05);
  hPad2->GetYaxis()->SetTitleSize(0.05);
  hPad2->GetYaxis()->SetTitleOffset(1.);


  hPad->GetYaxis()->SetRangeUser(effmin,1.05);
  hPad1->GetYaxis()->SetRangeUser(effmin,1.05);
  hPad2->GetYaxis()->SetRangeUser(effmin,1.05);

  pad2->cd();
  pad2->SetGridy();
  double tsize = (1./0.36)*hPad->GetYaxis()->GetTitleSize(); // 1./0.36
  TH1F *hPadr = (TH1F*) hPad->Clone("hPadr"); hPadr->GetYaxis()->SetRangeUser(1.-sfrange,1.+sfrange);
  hPadr->GetYaxis()->SetTitle("Scale Factor");
  hPadr->GetXaxis()->SetTitleSize(tsize);
  hPadr->GetXaxis()->SetLabelSize(tsize);
  hPadr->GetYaxis()->SetTitleSize(tsize);
  hPadr->GetYaxis()->SetLabelSize(tsize);
  hPadr->GetYaxis()->SetNdivisions(504,kTRUE);
  TH1F *hPadr_syst = (TH1F*) hPadr->Clone("hPadr_syst");hPadr_syst->GetYaxis()->SetRangeUser(1.-.1,1.+.1);
  TH1F *hPad1r = (TH1F*) hPad1->Clone("hPad1r"); hPad1r->GetYaxis()->SetRangeUser(1.-sfrange,1.+sfrange);
  hPad1r->GetYaxis()->SetTitle("Scale Factor");
  hPad1r->GetXaxis()->SetTitleSize(tsize);
  hPad1r->GetXaxis()->SetLabelSize(tsize);
  hPad1r->GetYaxis()->SetTitleSize(tsize);
  hPad1r->GetYaxis()->SetLabelSize(tsize);
  hPad1r->GetYaxis()->SetNdivisions(504,kTRUE);
  TH1F *hPad1r_syst = (TH1F*) hPad1r->Clone("hPad1r_syst");hPad1r_syst->GetYaxis()->SetRangeUser(1.-.1,1.+.1);
  TH1F *hPad2r = (TH1F*) hPad2->Clone("hPad2r"); hPad2r->GetYaxis()->SetRangeUser(1.-sfrange,1.+sfrange);
  hPad2r->GetYaxis()->SetTitle("Scale Factor");
  hPad2r->GetXaxis()->SetTitleSize(tsize);
  hPad2r->GetXaxis()->SetLabelSize(tsize);
  hPad2r->GetYaxis()->SetTitleSize(tsize);
  hPad2r->GetYaxis()->SetLabelSize(tsize);
  hPad2r->GetYaxis()->SetNdivisions(504,kTRUE);

  pad1->cd();

  TString header;
  for (int i=0; i<nbins_abseta; i++)
  {
     for (int k=0; k<nSyst; k++)
     {
        pad1->cd();
        hPad->Draw();
        TLatex *lt1 = new TLatex();
        lt1->SetNDC();

        char legs[512];
        TLegend *leg1 = new TLegend(0.43,0.05,0.66,0.43);
        leg1->SetFillStyle(0);
        leg1->SetFillColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.035);
        double ptmin = ((RooRealVar*) daPtData0[k][i]->get()->find("pt"))->getBinning().binLow(0);
        double etamin, etamax;
        if (daPtData0[k][i]->get()->find("abseta"))
        {
           etamin = ((RooRealVar*) daPtData0[k][i]->get()->find("abseta"))->getBinning().binLow(0);
           etamax = ((RooRealVar*) daPtData0[k][i]->get()->find("abseta"))->getBinning().binHigh(0);
           header = TString("#splitline{") + cutLegend + Form(" Efficiency}{(p^{#mu}_{T}>%.1f, |#eta| #in [%.1f, %.1f])}",ptmin,etamin,etamax);
        }
        else
        {
           etamin = ((RooRealVar*) daPtData0[k][i]->get()->find("eta"))->getBinning().binLow(0);
           etamax = ((RooRealVar*) daPtData0[k][i]->get()->find("eta"))->getBinning().binHigh(0);
           header = TString("#splitline{") + cutLegend + Form(" Efficiency}{(p^{#mu}_{T}>%.1f, #eta #in [%.1f, %.1f])}",ptmin,etamin,etamax);
        }
        leg1->SetHeader(header);
        sprintf(legs,"MC PYTHIA+EvtGen: %.4f^{ + %.3f}_{ - %.3f}", TrkAbsEta0[k][i][0], TrkAbsEta0[k][i][1], TrkAbsEta0[k][i][2]);
        leg1->AddEntry(ComPt0[k][i],legs,"pl");
        sprintf(legs,"Data: %.4f^{ + %.3f}_{ - %.3f}", TrkAbsEta1[k][i][0], TrkAbsEta1[k][i][1], TrkAbsEta1[k][i][2]);
        leg1->AddEntry(ComPt1[k][i],legs,"pl");
        leg1->Draw("same");

        ComPt0[k][i]->Draw("pz same");
        ComPt1[k][i]->Draw("pz same");

        lt1->SetTextSize(0.05);
        lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
        //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
        lt1->DrawLatex(0.43,0.54,collTag + "  #sqrt{s_{NN}} = 2.76 TeV");

        // now take care of the data/mc ratio panel
        c1->cd();
        // pad2->SetFrameFillStyle(4000);
        pad2->Draw();
        pad2->cd();
        hPadr->Draw();

        int nbins_mc = ComPt0[k][i]->GetN();
        int nbins = ComPt1[k][i]->GetN();
        double *xr = new double[nbins];
        double *yr = new double[nbins];
        double *xrlo = new double[nbins];
        double *yrlo = new double[nbins];
        double *xrhi = new double[nbins];
        double *yrhi = new double[nbins];

        TGraphAsymmErrors* ComPt0_forRatio = NULL;
        if (nbins_mc == nbins) {
           ComPt0_forRatio = ComPt0[k][i];
        }
        else {
           double* tntot = new double[nbins_mc];
           TString tag = absetaVar + "_bin0__pt_bin";
           for (int j=0; j<nbins_mc; j++) {
              RooFitResult *fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbGausPlusExpo/fitresults",j));
              if (!fitres) fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbPlusPoly1/fitresults",j));
              if (!fitres) fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbPlusPoly2/fitresults",j));
              if (!fitres) fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbGausPlusPoly/fitresults",j));
              if (!fitres) fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbPlusPoly2nd/fitresults",j));
              if (!fitres) fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbPlusPoly3rd/fitresults",j));
              if (!fitres) fitres = (RooFitResult*) f9[0]->Get(cutTag + "/" + ptTag + Form("%i/",i+1) + tag + Form("%i__cbPlusExpo/fitresults",j));
              if (!fitres) cerr << "ERROR I couldn't find the fit results! Expect a crash soon..." << endl;
              tntot[j] = ((RooRealVar*) fitres->floatParsFinal().find("numSignalAll"))->getVal();
              delete fitres;
           }
           ComPt0_forRatio = plotEff_1bin(daPtData0[k][i], 1, "pt", nbins_mc/nbins, tntot);
           delete[] tntot;
        }

        // here we assume that the mc uncertainty is negligible compared to the data one: simply scale everything by the central value.
        for (int j=0; j<nbins; j++)
        {
           xr[j] = ComPt1[k][i]->GetX()[j];
           xrlo[j] = ComPt1[k][i]->GetErrorXlow(j);
           xrhi[j] = ComPt1[k][i]->GetErrorXhigh(j);
           yr[j] = ComPt1[k][i]->GetY()[j]/ComPt0_forRatio->GetY()[j];
           yrlo[j] = ComPt1[k][i]->GetErrorYlow(j)/ComPt0_forRatio->GetY()[j];
           yrhi[j] = ComPt1[k][i]->GetErrorYhigh(j)/ComPt0_forRatio->GetY()[j];
        }
        TGraphAsymmErrors *gratio = new TGraphAsymmErrors(nbins,xr,yr,xrlo,xrhi,yrlo,yrhi);
        gratio->SetMarkerStyle(20);
        gratio->SetMarkerColor(kBlack);
        gratio->SetMarkerSize(1.0);
        gratio->SetLineColor(kBlack);
        gratio->SetLineWidth(1);
        gratio->Draw("pz same");

        // save (nominal only)
        if (k==0) {
           c1->SaveAs(cutTag + Form("Eff%i_",i) + collTag + "_RD_MC_PT.root");
           c1->SaveAs(cutTag + Form("Eff%i_",i) + collTag + "_RD_MC_PT.pdf");
        }

        // in case we are looking at muon Id + trigger: get the scale factor at the same time
#if defined MUIDTRG || defined STA
        pad1->cd();
        double ptmax = ((RooRealVar*) daPtData0[k][i]->get()->find("pt"))->getMax();
        TLatex tchi; tchi.SetNDC();
        tchi.SetTextSize(0.035);
        double chi2, pval; int dof;

        // fit data
        TF1 *fdata = new TF1("fdata","[0]*TMath::Erf((x-[1])/[2])",ptmin,ptmax);
        fdata->SetParNames("eff0","x0","m");
        // Initialize the normalization to the efficiency in the last point
        fdata->SetParameters(ComPt1[k][i]->GetX()[ComPt1[k][i]->GetN()-1],0.1,1.0);
        fdata->SetParLimits(0,0,1);
        fdata->SetParLimits(1,0.,10.);
        fdata->SetParLimits(2,0,10.);
        fdata->SetLineWidth(2);
        fdata->SetLineColor(kBlue);
        ComPt1[k][i]->Fit(fdata,"WRME");
        leg1->AddEntry(fdata,Form("%0.2f*TMath::Erf((x-%0.2f)/%0.2f)",fdata->GetParameter(0),fdata->GetParameter(1),fdata->GetParameter(2)),"pl");

        chi2 = ComPt1[k][i]->Chisquare(fdata);
        dof = ComPt1[k][i]->GetN() - fdata->GetNpar();
        pval = TMath::Prob(chi2,dof);
        tchi.SetTextColor(kBlue);
        tchi.DrawLatex(0.6,0.92,Form("#chi^{2}/dof = %.1f/%d (p-value: %.2f)",chi2,dof,pval));

        // fit mc
        TF1 *fmc = (TF1*) fdata->Clone("fmc");;
        // Initialize the normalization to the efficiency in the last point
        if (isPbPb) fmc->SetParameters(ComPt0[k][i]->GetX()[ComPt0[k][i]->GetN()-1],0.5,2.5); 
        else fmc->SetParameters(ComPt0[k][i]->GetX()[ComPt0[k][i]->GetN()-1],2.2,1.5); 
        fmc->SetLineColor(kRed);
        ComPt0[k][i]->Fit(fmc,"WRME");
        leg1->AddEntry(fmc,Form("%0.2f*TMath::Erf((x-%0.2f)/%0.2f)",fmc->GetParameter(0),fmc->GetParameter(1),fmc->GetParameter(2)),"pl");

        chi2 = ComPt0[k][i]->Chisquare(fmc);
        dof = ComPt0[k][i]->GetN() - fmc->GetNpar();
        pval = TMath::Prob(chi2,dof);
        tchi.SetTextColor(kRed);
        tchi.DrawLatex(0.6,0.88,Form("#chi^{2}/dof = %.1f/%d (p-value: %.2f)",chi2,dof,pval));

        leg1->Draw();

        // now the bottom panel
        pad2->cd();
        // hPadr->Draw();
        TF1 *fratio = new TF1("fratio","[0]*TMath::Erf((x-[1])/[2])/([3]*TMath::Erf((x-[4])/[5]))",ptmin,ptmax);
        fratio->SetParameters(
              fdata->GetParameter(0),fdata->GetParameter(1),fdata->GetParameter(2),
              fmc->GetParameter(0),fmc->GetParameter(1),fmc->GetParameter(2)
              );
        fratio->Draw("same");

        chi2 = gratio->Chisquare(fratio);
        dof = gratio->GetN() - fratio->GetNpar();
        pval = TMath::Prob(chi2,dof);
        tchi.SetTextColor(kBlack);
        tchi.SetTextSize(0.035*0.7/0.3);
        tchi.DrawLatex(0.6,0.8,Form("#chi^{2}/dof = %.1f/%d (p-value: %.2f)",chi2,dof,pval));

        // save (nominal only)
        if (k==0) {
           c1->SaveAs(cutTag + Form("SF%i_",i) + collTag + "_RD_MC_PT.root");
           c1->SaveAs(cutTag + Form("SF%i_",i) + collTag + "_RD_MC_PT.pdf");
        }

        // print the fit results to file
        file_sfs << "Data " << etamin << " " << etamax << endl;
        file_sfs << Form("%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",fdata->GetParameter(0),fdata->GetParameter(1),fdata->GetParameter(2)) << endl;
        file_sfs << "MC " << etamin << " " << etamax << endl;
        file_sfs << Form("%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",fmc->GetParameter(0),fmc->GetParameter(1),fmc->GetParameter(2)) << endl;
        file_sfs << endl;
     }

     // plot systematics
     TGraphAsymmErrors *graphssyst[nSyst];
     //data
     for (int k=0; k<nSyst; k++) graphssyst[k] = ComPt0[k][i];
     plotSysts(graphssyst, c1, pad1, hPad, pad2, hPadr_syst, header, Form("syst_data_pt_%i",i));
     //mc
     for (int k=0; k<nSyst; k++) graphssyst[k] = ComPt1[k][i];
     plotSysts(graphssyst, c1, pad1, hPad, pad2, hPadr_syst, header, Form("syst_mc_pt_%i",i));

     // toys study ///FIXME
     // if (doToys) toyStudy(ComPt1[i], ComPt0[i], fdata, fmc, cutTag + Form("toys%i_",i) + collTag + "_RD_MC_PT");
#endif // ifdef MUIDTRG or STA
  }

   //---------- This is for eta dependence
  TLegend *leg1 = new TLegend(0.43,0.21,0.66,0.43);
  TLatex *lt1 = new TLatex();
  for (int k=0; k<nSyst; k++)
  {
     pad1->cd();
     hPad1->Draw();

     ComEta0[k]->Draw("pz same");
     ComEta1[k]->Draw("pz same");

     lt1->SetNDC();
     char legs[512];
     leg1->SetFillStyle(0);
     leg1->SetFillColor(0);
     leg1->SetBorderSize(0);
     leg1->SetTextSize(0.035);
     double ptmin = ((RooRealVar*) daEtaData0[k]->get()->find("pt"))->getBinning().binLow(0);
     leg1->SetHeader(TString("#splitline{") + cutLegend + Form(" Efficiency}{(p^{#mu}_{T}>%.1f)}",ptmin));
     sprintf(legs,"MC PYTHIA+EvtGen: %.4f^{ + %.3f}_{ - %.3f}", Trk0[k][0], Trk0[k][1], Trk0[k][2]);
     leg1->AddEntry(ComPt0[k][0],legs,"pl");
     sprintf(legs,"Data: %.4f^{ + %.3f}_{ - %.3f}", Trk1[k][0], Trk1[k][1], Trk1[k][2]);
     leg1->AddEntry(ComPt1[k][0],legs,"pl");
     leg1->Draw("same");
     leg1->Draw("same");

     lt1->SetTextSize(0.05);
     lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
     //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
     lt1->DrawLatex(0.43,0.54,collTag + "  #sqrt{s_{NN}} = 2.76 TeV");

     // now take care of the data/mc ratio panel
     c1->cd();
     // pad2->SetFrameFillStyle(4000);
     pad2->Draw();
     pad2->cd();
     hPad1r->Draw();

     int nbins = ComEta0[k]->GetN();
     double *xr = new double[nbins];
     double *yr = new double[nbins];
     double *xrlo = new double[nbins];
     double *yrlo = new double[nbins];
     double *xrhi = new double[nbins];
     double *yrhi = new double[nbins];

     // here we assume that the mc uncertainty is negligible compared to the data one: simply scale everything by the central value.
     for (int j=0; j<nbins; j++)
     {
        xr[j] = ComEta1[k]->GetX()[j];
        xrlo[j] = ComEta1[k]->GetErrorXlow(j);
        xrhi[j] = ComEta1[k]->GetErrorXhigh(j);
        yr[j] = ComEta1[k]->GetY()[j]/ComEta0[k]->GetY()[j];
        yrlo[j] = ComEta1[k]->GetErrorYlow(j)/ComEta0[k]->GetY()[j];
        yrhi[j] = ComEta1[k]->GetErrorYhigh(j)/ComEta0[k]->GetY()[j];
     }
     TGraphAsymmErrors *gratio1 = new TGraphAsymmErrors(nbins,xr,yr,xrlo,xrhi,yrlo,yrhi);
     gratio1->SetMarkerStyle(20);
     gratio1->SetMarkerColor(kBlack);
     gratio1->SetMarkerSize(1.0);
     gratio1->SetLineColor(kBlack);
     gratio1->SetLineWidth(1);
     gratio1->Draw("pz same");

     if (k==0) {
        c1->SaveAs(cutTag + "Eff_" + collTag + "_RD_MC_Eta.root");
        c1->SaveAs(cutTag + "Eff_" + collTag + "_RD_MC_Eta.pdf");
     }
  }

  // plot systematics
  //data
  plotSysts(ComEta0, c1, pad1, hPad1, pad2, hPad1r_syst, header, "syst_data_eta");
  //mc
  plotSysts(ComEta1, c1, pad1, hPad1, pad2, hPad1r_syst, header, "syst_mc_eta");
  
  //-------- This is for centrality dependence
  if (isPbPb) {
     pad1->cd();
     hPad2->Draw();

     effCentMC->Draw("pz same");
     effCentData->Draw("pz same");

     leg1->Draw("same");

     lt1->SetTextSize(0.05);
     lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
     //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
     lt1->DrawLatex(0.43,0.54,collTag + "  #sqrt{s_{NN}} = 2.76 TeV");

     // now take care of the data/mc ratio panel
     c1->cd();
     // pad2->SetFrameFillStyle(4000);
     pad2->Draw();
     pad2->cd();
     hPad2r->Draw();

     int nbins2 = effCentMC->GetN();
     double* xr2 = new double[nbins2];
     double* yr2 = new double[nbins2];
     double* xr2lo = new double[nbins2];
     double* yr2lo = new double[nbins2];
     double* xr2hi = new double[nbins2];
     double* yr2hi = new double[nbins2];

     // here we assume that the mc uncertainty is negligible compared to the data one: simply scale everything by the central value.
     for (int j=0; j<nbins2; j++)
     {
        xr2[j] = effCentData->GetX()[j];
        xr2lo[j] = effCentData->GetErrorXlow(j);
        xr2hi[j] = effCentData->GetErrorXhigh(j);
        yr2[j] = effCentData->GetY()[j]/effCentMC->GetY()[j];
        yr2lo[j] = effCentData->GetErrorYlow(j)/effCentMC->GetY()[j];
        yr2hi[j] = effCentData->GetErrorYhigh(j)/effCentMC->GetY()[j];
     }
     TGraphAsymmErrors *gratio2 = new TGraphAsymmErrors(nbins2,xr2,yr2,xr2lo,xr2hi,yr2lo,yr2hi);
     gratio2->SetMarkerStyle(20);
     gratio2->SetMarkerColor(kBlack);
     gratio2->SetMarkerSize(1.0);
     gratio2->SetLineColor(kBlack);
     gratio2->SetLineWidth(1);
     gratio2->Draw("pz same");

     c1->SaveAs(cutTag + "Eff_" + collTag + "_RD_MC_Cent.root");
     c1->SaveAs(cutTag + "Eff_" + collTag + "_RD_MC_Cent.pdf");
  }


  // close the file for correction functions
#if defined MUIDTRG || defined STA
  file_sfs.close();
#endif
}

void formatTH1F(TH1* a, int b, int c, int d){
  a->SetLineWidth(2);
  a->SetLineStyle(c);
  a->SetMarkerSize(2);
  a->SetLineColor(b);
  a->SetMarkerColor(b);
  a->GetYaxis()->SetTitle("Efficiency");
  if(d == 1){
    a->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  }else if(d == 2){
    a->GetXaxis()->SetTitle("#eta");
  }else if(d == 3){
    a->GetXaxis()->SetTitle("rapidity");
  }else if(d == 4){
    a->GetXaxis()->SetTitle("Centrality");
  }
  a->GetXaxis()->CenterTitle(true);
  a->GetXaxis()->SetLabelSize(0.05);
  a->GetXaxis()->SetTitleSize(0.05);
  a->GetXaxis()->SetTitleOffset(0.9);
  a->GetYaxis()->SetLabelSize(0.05);
  a->GetYaxis()->SetTitleSize(0.05);
  a->GetYaxis()->SetTitleOffset(0.8);

}

void formatTLeg(TLegend* a){

  a->SetFillStyle(0);
  a->SetFillColor(0);
  a->SetBorderSize(0);
  a->SetTextSize(0.03);

}


void formatTGraph(TGraph* a, int b, int c, int d)
{

  a->SetMarkerStyle(c);
  a->SetMarkerColor(b);
  a->SetMarkerSize(1.0);
  a->SetLineColor(b);
  a->SetLineWidth(1);
  a->GetXaxis()->SetLabelSize(0.05);
  a->GetXaxis()->SetTitleSize(0.06);
  a->GetXaxis()->SetTitleOffset(0.9);
  a->GetYaxis()->SetTitle("Efficiency");
  a->GetXaxis()->CenterTitle();
  if(d == 1){ 
    a->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  }else if(d == 2){ 
    a->GetXaxis()->SetTitle("eta"); 
  }else if(d == 3){ 
    a->GetXaxis()->SetTitle("rapidity"); 
  }else if(d == 4){
    a->GetXaxis()->SetTitle("Centrality");
  }   
  a->GetYaxis()->SetRangeUser(0,1);
  a->GetXaxis()->SetRangeUser(0,10);
  a->GetYaxis()->SetLabelSize(0.05);
  a->GetYaxis()->SetTitleSize(0.05);
  a->GetYaxis()->SetTitleOffset(0.9);


}

TGraphAsymmErrors *plotEff_1bin(RooDataSet *a, int aa, const char* varx, int rebin, double *tntot){
  const RooArgSet *set = a->get();
  RooRealVar *xAx = (RooRealVar*)set->find(varx);
  RooRealVar *eff = (RooRealVar*)set->find("efficiency");

  int nbins = xAx->getBinning().numBins();

  double tx[nbins], txhi[nbins], txlo[nbins];
  double ty[nbins], tyhi[nbins], tylo[nbins];

  for (int i=0; i<nbins; i++) {
    a->get(i);
    ty[i] = eff->getVal();
    tx[i] = xAx->getVal();
    txhi[i] = fabs(xAx->getErrorHi());
    txlo[i] = fabs(xAx->getErrorLo()); 
    tyhi[i] = fabs(eff->getErrorHi());
    tylo[i] = fabs(eff->getErrorLo()); 
  }

  // do the rebinning here if necessary
  if (rebin>1)
  {
     if (nbins % rebin != 0) cout << "ERROR YOU ARE TRYING TO REBIN " << nbins << " BY " << rebin << " HOW IS THIS SUPPOSED TO WORK ???" << endl;

     nbins = nbins/rebin;
     for (int i=0; i<nbins; i++) {
        double ntot1 = tntot[2*i];
        double ntot2 = tntot[2*i+1];
        ty[i] = (ty[2*i]*ntot1 + ty[2*i+1]*ntot2) / (ntot1 + ntot2);
        tx[i] = (tx[2*i]*ntot1 + tx[2*i+1]*ntot2) / (ntot1 + ntot2);
        cout << "..." << ntot1 << " " << tx[2*i] << " " << ntot2 << " " << tx[2*i+1] << " " << tx[i] << endl;
        txhi[i] = (txhi[2*i]*ntot1 + txhi[2*i+1]*ntot2) / (ntot1 + ntot2);
        txlo[i] = (txlo[2*i]*ntot1 + txlo[2*i+1]*ntot2) / (ntot1 + ntot2); 
        tyhi[i] = (tyhi[2*i]*ntot1 + tyhi[2*i+1]*ntot2) / (ntot1 + ntot2);
        tylo[i] = (tylo[2*i]*ntot1 + tylo[2*i+1]*ntot2) / (ntot1 + ntot2); 
     }
  }


  cout<<"NBins : "<<nbins<<endl;

  const double *x = tx;
  const double *xhi = txhi;
  const double *xlo = txlo;
  const double *y = ty;
  const double *yhi = tyhi;
  const double *ylo = tylo;


  TGraphAsymmErrors *b = new TGraphAsymmErrors();
  if(aa == 1) {*b = TGraphAsymmErrors(nbins,x,y,xlo,xhi,ylo,yhi);}
  if(aa == 0) {*b = TGraphAsymmErrors(nbins,x,y,0,0,ylo,yhi);}
  b->SetMaximum(1.1);
  b->SetMinimum(0.0);
  b->SetMarkerStyle(20);
  b->SetMarkerColor(kRed+2);
  b->SetMarkerSize(1.0);
  b->SetTitle("");
  b->GetXaxis()->SetTitleSize(0.1);
  b->GetYaxis()->SetTitleSize(0.1);
  b->GetXaxis()->SetTitle("#eta");
  b->GetYaxis()->SetTitle("Efficiency");
  b->GetXaxis()->CenterTitle();
  //b->Draw("apz");

  for (int i=0; i<nbins; i++) {
    cout << x[i] << " " << y[i] << " " << yhi[i] << " " << ylo[i] << endl;
  }

  return b;

}

vector<TGraphAsymmErrors*> plotEff_Nbins(RooDataSet *a, int aa, const char* varx, const char* var2){
  const RooArgSet *set = a->get();
  RooRealVar *xAx = (RooRealVar*)set->find(varx);
  RooRealVar *abseta = (RooRealVar*)set->find(var2);
  RooRealVar *eff = (RooRealVar*)set->find("efficiency");

  int nbins = xAx->getBinning().numBins();
  const int nbins2 = abseta->getBinning().numBins();
  // cout << nbins << " " << nbins2 << endl;

  double **tx = new double*[nbins2], **txhi = new double*[nbins2], **txlo = new double*[nbins2];
  double **ty = new double*[nbins2], **tyhi = new double*[nbins2], **tylo = new double*[nbins2];
  for (int i=0; i<nbins2; i++)
  {
     tx[i] = new double[nbins];
     txlo[i] = new double[nbins];
     txhi[i] = new double[nbins];
     ty[i] = new double[nbins];
     tylo[i] = new double[nbins];
     tyhi[i] = new double[nbins];
  }

  for (int i=0; i<nbins*nbins2; i++) {
    a->get(i);
    int bin1 = xAx->getBinning().binNumber(xAx->getVal());
    int bin2 = abseta->getBinning().binNumber(abseta->getVal());
    ty[bin2][bin1] = eff->getVal();
    tx[bin2][bin1] = xAx->getVal();
    txhi[bin2][bin1] = fabs(xAx->getErrorHi());
    txlo[bin2][bin1] = fabs(xAx->getErrorLo()); 
    tyhi[bin2][bin1] = fabs(eff->getErrorHi());
    tylo[bin2][bin1] = fabs(eff->getErrorLo()); 
  }

  // cout<<"NBins pt: "<<nbins<<endl;

  vector<TGraphAsymmErrors*> result;
  for (int i=0; i<nbins2; i++)
  {
     const double *x = tx[i];
     const double *xhi = txhi[i];
     const double *xlo = txlo[i];
     const double *y = ty[i];
     const double *yhi = tyhi[i];
     const double *ylo = tylo[i];
     TGraphAsymmErrors *b = new TGraphAsymmErrors();
     if(aa == 1) {*b = TGraphAsymmErrors(nbins,x,y,xlo,xhi,ylo,yhi);}
     if(aa == 0) {*b = TGraphAsymmErrors(nbins,x,y,0,0,ylo,yhi);}
     b->SetMaximum(1.1);
     b->SetMinimum(0.0);
     b->SetMarkerStyle(20);
     b->SetMarkerColor(kRed+2);
     b->SetMarkerSize(1.0);
     b->SetTitle("");
     b->GetXaxis()->SetTitleSize(0.05);
     b->GetYaxis()->SetTitleSize(0.05);
     b->GetXaxis()->SetTitle("p_{T} [GeV/c]");
     b->GetYaxis()->SetTitle("Efficiency");
     b->GetXaxis()->CenterTitle();
     //b->Draw("apz");
     result.push_back(b);
  }

  // for (int i=0; i<nbins; i++) {
  //   cout << x[i] << " " << y[i] << " " << yhi[i] << " " << ylo[i] << endl;
  // }

  return result;

}

TGraphAsymmErrors *plotEffCent(RooDataSet **a1,int aa){
  const RooArgSet *setn[nCentBins];
  RooRealVar *effn[nCentBins];
  for (int i=0; i<nCentBins; i++)
  {
     setn[i] = a1[i]->get();
     effn[i] = (RooRealVar*)setn[i]->find("efficiency");
  }

  const int nbins = nCentBins;

  double ty[nCentBins], tyhi[nCentBins], tylo[nCentBins];

  for (int i=0; i<nCentBins; i++)
  {
     ty[i] = effn[i]->getVal();
     tyhi[i] = fabs(effn[i]->getErrorHi());
     tylo[i] = fabs(effn[i]->getErrorLo());
  }

  const double *x = tx;
  const double *y = ty;
  const double *yhi = tyhi;
  const double *ylo = tylo;

  TGraphAsymmErrors *b = new TGraphAsymmErrors();
  if(aa == 1) {*b = TGraphAsymmErrors(nbins,x,y,0,0,ylo,yhi);}
  if(aa == 0) {*b = TGraphAsymmErrors(nbins,x,y,0,0,ylo,yhi);}
  b->SetMaximum(1.1);
  b->SetMinimum(0.0);
  b->SetMarkerStyle(20);
  b->SetMarkerColor(kRed+2);
  b->SetMarkerSize(1.0);
  b->SetTitle("");
  b->GetXaxis()->SetTitleSize(0.05);
  b->GetYaxis()->SetTitleSize(0.05);
  b->GetXaxis()->SetTitle("centrality");
  b->GetYaxis()->SetTitle("efficiency");
  b->GetXaxis()->CenterTitle();
  //b->Draw("apz");

  for (int i=0; i<nbins; i++) {
    cout << x[i] << " " << y[i] << " " << yhi[i] << " " << ylo[i] << endl;
  }

  return b;

}



TH2F *plotEff2D(RooDataSet *a, TString b){
  const RooArgSet *set = a->get();
  RooRealVar *yAx = (RooRealVar*)set->find("pt");
  RooRealVar *xAx = (RooRealVar*)set->find("eta");
  RooRealVar *eff = (RooRealVar*)set->find("efficiency");

  //const int xnbins = xAx->getBinning().numBins();
  //const int ynbins = yAx->getBinning().numBins();

  const double *xvbins = xAx->getBinning().array();
  const double *yvbins = yAx->getBinning().array();

  TH2F* h = new TH2F(b, "", xAx->getBinning().numBins(), xvbins, yAx->getBinning().numBins(), yvbins);

  gStyle->SetPaintTextFormat("5.2f");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPalette(1);
  h->SetOption("colztexte");
  h->GetZaxis()->SetRangeUser(-0.001,1.001);
  h->SetStats(kFALSE);
  h->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h->GetXaxis()->SetTitle("#eta");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetXaxis()->SetTitleOffset(0.9);
  for(int i=0; i<a->numEntries(); i++){
    a->get(i);
    h->SetBinContent(h->FindBin(xAx->getVal(), yAx->getVal()), eff->getVal());
    h->SetBinError(h->FindBin(xAx->getVal(), yAx->getVal()), (eff->getErrorHi()-eff->getErrorLo())/2.);
  }

  return h;

}


void CalEffErr(TGraph *a, double *b){
    const int nbins = 100;
    double x[nbins], y[nbins];
    double sum = 0, errHighSum = 0, errLowSum = 0, sqSumHigh = 0, sqSumLow = 0;
    //double b[3] = 0;
    
    int nBins = a->GetN();
    for(int i=0;i<a->GetN();i++){
        a->GetPoint(i,x[i],y[i]);
        //cout<<"Eff x = "<<x[i]<<" y = "<<y[i]<<endl;
        double eHigh = a->GetErrorYhigh(i);
        double eLow = a->GetErrorYlow(i);
        //cout<<"Err high = "<<eHigh<<" low = "<<eLow<<endl;
        sum += y[i];
        errHighSum += eHigh;
        sqSumHigh += eHigh*eHigh;
        errLowSum += eLow;
        sqSumLow += eLow*eLow;
    }
    b[0] = sum/nBins;
    b[1] = sqrt(sqSumHigh)/nBins;
    b[2] = sqrt(sqSumLow)/nBins;
    //cout<<"b1 : "<<b[0]<<", b2 : "<<b[1]<<", b3 : "<<b[2]<<endl;
    
    cout<<b[0]<<"^{"<<b[1]<<"}_{"<<b[2]<<"}"<<endl;
    //return b[3];
}

void CalEffErr(vector<TGraphAsymmErrors*> a, double **b){
    const int nbins = 100;
    const int vsize = a.size();

    for (int vbin=0; vbin<vsize; vbin++)
    {
       double x[nbins], y[nbins];
       double sum = 0, errHighSum = 0, errLowSum = 0, sqSumHigh = 0, sqSumLow = 0;
       //double b[3] = 0;

       int nBins = a[vbin]->GetN();
       for(int i=0;i<a[vbin]->GetN();i++){
          a[vbin]->GetPoint(i,x[i],y[i]);
          //cout<<"Eff x = "<<x[i]<<" y = "<<y[i]<<endl;
          double eHigh = a[vbin]->GetErrorYhigh(i);
          double eLow = a[vbin]->GetErrorYlow(i);
          //cout<<"Err high = "<<eHigh<<" low = "<<eLow<<endl;
          sum += y[i];
          errHighSum += eHigh;
          sqSumHigh += eHigh*eHigh;
          errLowSum += eLow;
          sqSumLow += eLow*eLow;
       }
       b[vbin][0] = sum/nBins;
       b[vbin][1] = sqrt(sqSumHigh)/nBins;
       b[vbin][2] = sqrt(sqSumLow)/nBins;
       //cout<<"b1 : "<<b[0]<<", b2 : "<<b[1]<<", b3 : "<<b[2]<<endl;

       cout<<b[vbin][0]<<"^{"<<b[vbin][1]<<"}_{"<<b[vbin][2]<<"}"<<endl;
    }
    //return b[3];
}

void plotSysts(TGraphAsymmErrors *graphs[nSyst], TCanvas *c1, TPad *p1, TH1F *h1, TPad *pr, TH1F *hr, TString header, TString tag){
        p1->cd();
        h1->Draw();
        TLatex *lt1 = new TLatex();
        lt1->SetNDC();

        char legs[512];
        TLegend *leg1 = new TLegend(0.43,0.05,0.66,0.43);
        leg1->SetFillStyle(0);
        leg1->SetFillColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.035);
        leg1->SetHeader(header);
        for (int k=0; k<nSyst; k++) {
           graphs[k]->SetMarkerColor(k<4 ? k+1 : k+2);
           graphs[k]->SetMarkerStyle(k+20);
           graphs[k]->SetLineColor(k<4 ? k+1 : k+2);
           leg1->AddEntry(graphs[k],systName[k],"pl");
           graphs[k]->Draw("pz same");
        }
        leg1->Draw("same");

        // lt1->SetTextSize(0.05);
        // lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
        // //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
        // lt1->DrawLatex(0.43,0.54,collTag + "  #sqrt{s_{NN}} = 2.76 TeV");

        // now take care of the data/mc ratio panel
        c1->cd();
        // pad2->SetFrameFillStyle(4000);
        pr->Draw();
        pr->cd();
        hr->Draw();

        int nbins_mc = graphs[0]->GetN();
        int nbins = graphs[0]->GetN();
        double **xr = new double*[nSyst];
        double **yr = new double*[nSyst];
        double **xrlo = new double*[nSyst];
        double **yrlo = new double*[nSyst];
        double **xrhi = new double*[nSyst];
        double **yrhi = new double*[nSyst];
        for (int j=0; j<nSyst; j++) {
           xr[j] = new double[nbins];
           yr[j] = new double[nbins];
           xrlo[j] = new double[nbins];
           yrlo[j] = new double[nbins];
           xrhi[j] = new double[nbins];
           yrhi[j] = new double[nbins];
        }

        TGraphAsymmErrors *gratio[nSyst];

        TGraphAsymmErrors* ComPt0_forRatio = NULL;
        ComPt0_forRatio = graphs[0];

        // here we assume that the mc uncertainty is negligible compared to the data one: simply scale everything by the central value.
        for (int k=1; k<nSyst; k++) {
           for (int j=0; j<nbins; j++)
           {
              xr[k][j] = graphs[k]->GetX()[j];
              xrlo[k][j] = graphs[k]->GetErrorXlow(j);
              xrhi[k][j] = graphs[k]->GetErrorXhigh(j);
              yr[k][j] = graphs[k]->GetY()[j]/ComPt0_forRatio->GetY()[j];
              yrlo[k][j] = graphs[k]->GetErrorYlow(j)/ComPt0_forRatio->GetY()[j];
              yrhi[k][j] = graphs[k]->GetErrorYhigh(j)/ComPt0_forRatio->GetY()[j];
           }
           gratio[k] = new TGraphAsymmErrors(nbins,xr[k],yr[k],xrlo[k],xrhi[k],yrlo[k],yrhi[k]);
           gratio[k]->SetMarkerStyle(20+k);
           gratio[k]->SetMarkerColor(k<4 ? k+1 : k+2);
           gratio[k]->SetMarkerSize(1.0);
           gratio[k]->SetLineColor(k<4 ? k+1 : k+2);
           gratio[k]->SetLineWidth(1);
           gratio[k]->Draw("pz same");
        }

        // save
        c1->SaveAs(tag + ".root");
        c1->SaveAs(tag + ".pdf");
}
