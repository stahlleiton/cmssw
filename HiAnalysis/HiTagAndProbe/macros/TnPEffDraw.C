#include <iostream>
#include "TSystem.h"
#include "TTree.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TPave.h"
#include "TText.h"
#include <sstream>
#include <string.h>

#include "TROOT.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"

#include "TStyle.h"
#include "TLatex.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TPostScript.h"


using namespace RooFit;
using namespace std;

const int nCentBins = 7;
const int nAbsEtaBins = 1;
const double tx[nCentBins]={381.23,330.35,261.37,158.7,86.3,53.5,16.2};
TString cutTag("MuonTrg"); // MuonIDTrg
TString cutLegend("Muon ID");
TString etaTag("PassingGlb_eta");
TString absetaTag("PassingGlb_1bin"); // "PassingGlb_abseta"
TString absetaVar("eta");
TString ptTag("PassingGlb_pt");
TString allTag("PassingGlb_1bin");
const char* fCentMCNames[nCentBins] = 
{
   "fitdata/Trg_cent_20150907/res_Cent0005_/data_Trg_Cent0005_.root",
   "fitdata/Trg_cent_20150907/res_Cent0510_/data_Trg_Cent0510_.root",
   "fitdata/Trg_cent_20150907/res_Cent1020_/data_Trg_Cent1020_.root",
   "fitdata/Trg_cent_20150907/res_Cent2040_/data_Trg_Cent2040_.root",
   "fitdata/Trg_cent_20150907/res_Cent4050_/data_Trg_Cent4050_.root",
   "fitdata/Trg_cent_20150907/res_Cent5060_/data_Trg_Cent5060_.root",
   "fitdata/Trg_cent_20150907/res_Cent60100_/data_Trg_Cent60100_.root",
};
const char* fCentDataNames[nCentBins] = 
{
   "fitdata/Trg_cent_20150907/res_Cent0005_/data_Trg_Cent0005_.root",
   "fitdata/Trg_cent_20150907/res_Cent0510_/data_Trg_Cent0510_.root",
   "fitdata/Trg_cent_20150907/res_Cent1020_/data_Trg_Cent1020_.root",
   "fitdata/Trg_cent_20150907/res_Cent2040_/data_Trg_Cent2040_.root",
   "fitdata/Trg_cent_20150907/res_Cent4050_/data_Trg_Cent4050_.root",
   "fitdata/Trg_cent_20150907/res_Cent5060_/data_Trg_Cent5060_.root",
   "fitdata/Trg_cent_20150907/res_Cent60100_/data_Trg_Cent60100_.root",
};
const char* fDataName = "fitdata/Trg_4etabins_20150908/tnp_Ana_Data_PbPb_MuonTrg_AllMB.root";
const char* fMCName = "fitmc/Trg_20150908/tnp_Ana_MC_PbPb_MuonTrg_AllMB.root";

// Function Define
TH2F *plotEff2D(RooDataSet *a, TString b);
vector<TGraphAsymmErrors*> plotEff_Nbins(RooDataSet *a, int aa, const char* varx, const char* var2);
TGraphAsymmErrors *plotEff_1bin(RooDataSet *a, int aa, const char* varx);
TGraphAsymmErrors *plotEffCent(RooDataSet **a1,int aa);
void formatTH1F(TH1* a, int b, int c, int d);
void formatTGraph(TGraph* a, int b, int c, int d);
void formatTLeg(TLegend* a);
void CalEffErr(TGraph *a, double *b);
void CalEffErr(vector<TGraphAsymmErrors*> a, double **b);

// From here you need to set up your environments.
void TnPEffDraw() {

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

  for (int i=0; i<nCentBins; i++)
  {
     fCentMC[i] = new TFile(fCentMCNames[i]);
     fCentData[i] = new TFile(fCentDataNames[i]);
  }
  
  //data and MC root files as well as single bin for integrated efficiency
  TFile *f9 = new TFile(fDataName);
  TFile *f10 = new TFile(fMCName);

  TCanvas *c1 = new TCanvas("c1","",700,600);

  vector<RooDataSet*> daPtData0, daPtData1;

  if (nAbsEtaBins==1)
  {
     daPtData0.push_back((RooDataSet*)f9->Get(cutTag + "/" + ptTag + "/fit_eff")); 
     daPtData1.push_back((RooDataSet*)f10->Get(cutTag + "/" + ptTag + "/fit_eff")); 
  }
  else
  {
     for (int i=0; i<nAbsEtaBins; i++)
     {
        daPtData0.push_back((RooDataSet*)f9->Get(cutTag + "/" + ptTag + Form("%i/fit_eff",i+1))); 
        daPtData1.push_back((RooDataSet*)f10->Get(cutTag + "/" + ptTag + Form("%i/fit_eff",i+1))); 
     }
  }

  vector<TGraphAsymmErrors*> ComPt0, ComPt1;
  
  for (unsigned int i=0; i<daPtData0.size(); i++)
  {
  cout << daPtData0[i] << endl;
     ComPt0.push_back(plotEff_1bin(daPtData0[i], 1, "pt"));
     ComPt1.push_back(plotEff_1bin(daPtData1[i], 1, "pt"));
  }


  RooDataSet *daEtaData0 = (RooDataSet*)f9->Get(cutTag + "/" + etaTag + "/fit_eff"); 
  RooDataSet *daEtaData1 = (RooDataSet*)f10->Get(cutTag + "/" + etaTag + "/fit_eff"); 

  TGraphAsymmErrors *ComEta0 = plotEff_1bin(daEtaData0, 1, "eta");
  TGraphAsymmErrors *ComEta1 = plotEff_1bin(daEtaData1, 1, "eta");

  RooDataSet *daCentMC1Bin[nCentBins];
  RooDataSet *daCentData1Bin[nCentBins];
  for (int i=0; i<nCentBins; i++)
  {
     daCentMC1Bin[i] = (RooDataSet*) fCentMC[i]->Get(cutTag + "/" + allTag + "/fit_eff");
     daCentData1Bin[i] = (RooDataSet*) fCentData[i]->Get(cutTag + "/" + allTag + "/fit_eff");
  }

  RooDataSet *daPtMC1Bin0 = (RooDataSet*)f9->Get(cutTag + "/" + allTag + "/fit_eff"); 
  RooDataSet *daPtData1Bin0 = (RooDataSet*)f10->Get(cutTag + "/" + allTag + "/fit_eff"); 
  RooDataSet *daAbsEtaMC1 = (RooDataSet*)f9->Get(cutTag + "/" + absetaTag + "/fit_eff"); 
  RooDataSet *daAbsEtaData1 = (RooDataSet*)f10->Get(cutTag + "/" + absetaTag + "/fit_eff"); 

  TGraphAsymmErrors* Com0Pt0 = plotEff_1bin(daPtMC1Bin0,0,"eta");
  TGraphAsymmErrors* Com0Pt1 = plotEff_1bin(daPtData1Bin0,0,"eta");
  vector<TGraphAsymmErrors*> Com0AbsEta0 = plotEff_Nbins(daAbsEtaMC1,0,"pt",absetaVar);
  vector<TGraphAsymmErrors*> Com0AbsEta1 = plotEff_Nbins(daAbsEtaData1,0,"pt",absetaVar);

  TGraphAsymmErrors *effCentMC = plotEffCent(daCentMC1Bin, 1);
  TGraphAsymmErrors *effCentData = plotEffCent(daCentData1Bin, 1);

  effCentMC->SetMarkerStyle(20);
  effCentMC->SetMarkerSize(1.4);
  effCentMC->SetMarkerColor(kRed+1);
  effCentMC->SetLineColor(kRed+1);
  effCentData->SetMarkerStyle(25);
  effCentData->SetMarkerSize(1.4);
  effCentData->SetMarkerColor(kBlue+1);
  effCentData->SetLineColor(kBlue+1);

  int nbins_abseta = ComPt0.size();
  for (int i=0; i<nbins_abseta; i++)
  {
     ComPt0[i]->SetMarkerStyle(20);
     ComPt0[i]->SetMarkerSize(1.4);
     ComPt0[i]->SetMarkerColor(kRed+1);
     ComPt0[i]->SetLineColor(kRed+1);
     ComPt1[i]->SetMarkerStyle(25);
     ComPt1[i]->SetMarkerSize(1.4);
     ComPt1[i]->SetMarkerColor(kBlue+1);
     ComPt1[i]->SetLineColor(kBlue+1);
  }

  ComEta0->SetMarkerStyle(20);
  ComEta0->SetMarkerSize(1.4);
  ComEta0->SetMarkerColor(kRed+1);
  ComEta0->SetLineColor(kRed+1);
  ComEta1->SetMarkerStyle(25);
  ComEta1->SetMarkerSize(1.4);
  ComEta1->SetMarkerColor(kBlue+1);
  ComEta1->SetLineColor(kBlue+1);


  // 1 Bin
  double Trk0[4];
  double Trk1[4];
  double** TrkAbsEta0 = new double*[nbins_abseta];
  double** TrkAbsEta1 = new double*[nbins_abseta];
  for (int i=0; i<nbins_abseta; i++) 
  {
     TrkAbsEta0[i] = new double[4];
     TrkAbsEta1[i] = new double[4];
  }
  CalEffErr(Com0Pt0, Trk0);
  CalEffErr(Com0Pt1, Trk1);
  CalEffErr(Com0AbsEta0, TrkAbsEta0);
  CalEffErr(Com0AbsEta1, TrkAbsEta1);




  TH1F *hPad = new TH1F("hPad",";p^{#mu}_{T} [GeV/c];Single #mu Efficiency",5,0,30);
  TH1F *hPad1 = new TH1F("hPad1",";#eta^{#mu} ;Single #mu Efficiency",5,-2.4,2.4);
  TH1F *hPad2 = new TH1F("hPad2",";N_{part} ;Single #mu Efficiency",5,0,400);
  hPad->GetXaxis()->CenterTitle();
  hPad1->GetXaxis()->CenterTitle();
  hPad2->GetXaxis()->CenterTitle();
  hPad->GetXaxis()->SetLabelSize(0.05);
  hPad->GetXaxis()->SetTitleSize(0.05);
  hPad->GetXaxis()->SetTitleOffset(1.2);
  hPad->GetYaxis()->SetLabelSize(0.05);
  hPad->GetYaxis()->SetTitleSize(0.05);
  hPad->GetYaxis()->SetTitleOffset(1.);
  hPad1->GetXaxis()->SetLabelSize(0.05);
  hPad1->GetXaxis()->SetTitleSize(0.05);
  hPad1->GetXaxis()->SetTitleOffset(1.2);
  hPad1->GetYaxis()->SetLabelSize(0.05);
  hPad1->GetYaxis()->SetTitleSize(0.05);
  hPad1->GetYaxis()->SetTitleOffset(1.);
  hPad2->GetXaxis()->SetLabelSize(0.05);
  hPad2->GetXaxis()->SetTitleSize(0.05);
  hPad2->GetXaxis()->SetTitleOffset(1.2);
  hPad2->GetYaxis()->SetLabelSize(0.05);
  hPad2->GetYaxis()->SetTitleSize(0.05);
  hPad2->GetYaxis()->SetTitleOffset(1.);


  hPad->GetYaxis()->SetRangeUser(0.,1.05);
  hPad1->GetYaxis()->SetRangeUser(0.,1.05);
  hPad2->GetYaxis()->SetRangeUser(0.,1.05);


  for (int i=0; i<nbins_abseta; i++)
  {
     hPad->Draw();
     TLatex *lt1 = new TLatex();
     lt1->SetNDC();

     char legs[512];
     TLegend *leg1 = new TLegend(0.43,0.21,0.66,0.43);
     leg1->SetFillStyle(0);
     leg1->SetFillColor(0);
     leg1->SetBorderSize(0);
     leg1->SetTextSize(0.035);
     double ptmin = ((RooRealVar*) daPtData0[i]->get()->find("pt"))->getBinning().binLow(0);
     double etamin = ((RooRealVar*) daPtData0[i]->get()->find("abseta"))->getBinning().binLow(0);
     double etamax = ((RooRealVar*) daPtData0[i]->get()->find("abseta"))->getBinning().binHigh(0);
     leg1->SetHeader(TString("#splitline{") + cutLegend + Form(" Efficiency}{(p^{#mu}_{T}>%.1f, |#eta| #in [%.1f, %.1f])}",ptmin,etamin,etamax));
     sprintf(legs,"MC PYTHIA+EvtGen: %.4f^{ + %.3f}_{ - %.3f}", TrkAbsEta0[i][0], TrkAbsEta0[i][1], TrkAbsEta0[i][2]);
     leg1->AddEntry(ComPt0[i],legs,"pl");
     sprintf(legs,"Data: %.4f^{ + %.3f}_{ - %.3f}", TrkAbsEta1[i][0], TrkAbsEta1[i][1], TrkAbsEta1[i][2]);
     leg1->AddEntry(ComPt1[i],legs,"pl");
     leg1->Draw("same");

     ComPt0[i]->Draw("pz same");
     ComPt1[i]->Draw("pz same");

     lt1->SetTextSize(0.05);
     lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
     //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
     lt1->DrawLatex(0.43,0.54,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");

     c1->SaveAs(cutTag + Form("Eff%i_PbPb_RD_MC_PT.png",i));
     c1->SaveAs(cutTag + Form("Eff%i_PbPb_RD_MC_PT.pdf",i));
  }

   //---------- This is for eta dependence
  hPad1->Draw();

  ComEta0->Draw("pz same");
  ComEta1->Draw("pz same");

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  char legs[512];
  TLegend *leg1 = new TLegend(0.43,0.21,0.66,0.43);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  double ptmin = ((RooRealVar*) daEtaData0->get()->find("pt"))->getBinning().binLow(0);
  leg1->SetHeader(TString("#splitline{") + cutLegend + Form(" Efficiency}{(p^{#mu}_{T}>%.1f)}",ptmin));
  sprintf(legs,"MC PYTHIA+EvtGen: %.4f^{ + %.3f}_{ - %.3f}", Trk0[0], Trk0[1], Trk0[2]);
  leg1->AddEntry(ComPt0[0],legs,"pl");
  sprintf(legs,"Data: %.4f^{ + %.3f}_{ - %.3f}", Trk1[0], Trk1[1], Trk1[2]);
  leg1->AddEntry(ComPt1[0],legs,"pl");
  leg1->Draw("same");
  leg1->Draw("same");

  lt1->SetTextSize(0.05);
  lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
  //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
  lt1->DrawLatex(0.43,0.54,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");

  c1->SaveAs(cutTag + "Eff_PbPb_RD_MC_Eta.png");
  c1->SaveAs(cutTag + "Eff_PbPb_RD_MC_Eta.pdf");
  
  //-------- This is for cent dependence
  hPad2->Draw();

  effCentMC->Draw("pz same");
  effCentData->Draw("pz same");

  leg1->Draw("same");

  lt1->SetTextSize(0.05);
  lt1->DrawLatex(0.43,0.60,"CMS Preliminary");
  //lt1->DrawLatex(0.43,0.54,"pp  #sqrt{s} = 2.76 TeV");
  lt1->DrawLatex(0.43,0.54,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");

  c1->SaveAs(cutTag + "Eff_PbPb_RD_MC_Cent.png");
  c1->SaveAs(cutTag + "Eff_PbPb_RD_MC_Cent.pdf");


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

TGraphAsymmErrors *plotEff_1bin(RooDataSet *a, int aa, const char* varx){
  const RooArgSet *set = a->get();
  RooRealVar *xAx = (RooRealVar*)set->find(varx);
  RooRealVar *eff = (RooRealVar*)set->find("efficiency");

  const int nbins = xAx->getBinning().numBins();

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

  const int nbins = xAx->getBinning().numBins();
  const int nbins2 = abseta->getBinning().numBins();
  cout << nbins << " " << nbins2 << endl;

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

  cout<<"NBins pt: "<<nbins<<endl;


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

