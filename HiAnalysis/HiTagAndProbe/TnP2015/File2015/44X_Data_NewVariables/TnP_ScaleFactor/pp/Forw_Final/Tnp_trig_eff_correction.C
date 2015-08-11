// This is macro to calculate scale factor using T&P method
// Editor : Dong Ho Moon

#include <iostream>
#include <string>

#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

using namespace std;

void Tnp_trig_eff_correction(string dataFile="Jpsi_Regit_RD_2014_Trg2_Eff.root", string mcFile="Jpsi_Regit_MC_2014_Trg2_Eff.root", bool savePlots=true)
{

  int etaRegion = 2; // 1 : mid, 2 : forward

  gROOT->Macro("~/rootlogon.C");
  gStyle->SetOptFit(0);
  TGraphAsymmErrors *grEffRD;
  TGraphAsymmErrors *grEffMC;
  dataFile="tnp_pp_RD_Forw.root";
  mcFile="tnp_pp_MC_Forw.root";

  TFile *tnpRD = new TFile(dataFile.c_str(),"READ");
  grEffRD = (TGraphAsymmErrors*) tnpRD->Get("gRD_pT")->Clone("grEffRD");
  //  tnpRD->Close();
  grEffRD->SetMarkerColor(kBlue);
  grEffRD->SetMarkerStyle(21);
  TFile *tnpMC = new TFile(mcFile.c_str(),"READ");
  grEffMC = (TGraphAsymmErrors*) tnpMC->Get("gMC_pT")->Clone("grEffMC");
  //grEffMC = (TGraphAsymmErrors*) tnpMC->Get("gTrg_pT_MC")->Clone("grEffMC");
  //  tnpMC->Close();
  grEffMC->SetMarkerColor(kRed);
  grEffMC->SetMarkerStyle(20);

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",0.0,20);

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();

  f1->SetParNames("eff0","x0","m");
  f1->SetParameters(0.8,0.1,1.0); // 0.8, 0.2, 2.5
  f1->SetLineWidth(2);

  TH1F *hPad = new TH1F("hPad",";p_{T} GeV/c;Trigger X Muon Id Efficiency",10,0,20);
  TH1F *hPad1 = new TH1F("hPad1",";p_{T} GeV/c;Trigger X Muon Id Efficiency",10,0,20);
  grEffRD->Fit(f1,"RME");
  f1->SetParameters(0.9,0.5,2.5);
  grEffMC->Fit(f1,"WRME");

  grEffRD->GetFunction("f1")->SetLineColor(kBlue);
  grEffMC->GetFunction("f1")->SetLineColor(kRed);

  cout<<"%%% Fit Summary for Data %%%"<<endl;
  cout<<"eff0 : "<<grEffRD->GetFunction("f1")->GetParameter(0)<<endl;
  cout<<"x0 : "<<grEffRD->GetFunction("f1")->GetParameter(1)<<endl;
  cout<<"m : "<<grEffRD->GetFunction("f1")->GetParameter(2)<<endl;
  char fnName1[512];
  sprintf(fnName1,"%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",grEffRD->GetFunction("f1")->GetParameter(0),grEffRD->GetFunction("f1")->GetParameter(1),grEffRD->GetFunction("f1")->GetParameter(2));
  cout<<fnName1<<endl;

  cout<<"%%% Fit Summary for MC %%%"<<endl;
  cout<<"eff0 : "<<grEffMC->GetFunction("f1")->GetParameter(0)<<endl;
  cout<<"x0 : "<<grEffMC->GetFunction("f1")->GetParameter(1)<<endl;
  cout<<"m : "<<grEffMC->GetFunction("f1")->GetParameter(2)<<endl;

  char fnName2[512];
  sprintf(fnName2,"%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",grEffMC->GetFunction("f1")->GetParameter(0),grEffMC->GetFunction("f1")->GetParameter(1),grEffMC->GetFunction("f1")->GetParameter(2));
  cout<<fnName2<<endl;

  cout<<"%%%% Final Correction Function : "<<"("<<fnName1<<")"<<"/"<<"("<<fnName2<<")"<<endl;



  TCanvas *c1 = new TCanvas("c1","c1");
  hPad->Draw();
  grEffRD->Draw("PZ");
  grEffMC->Draw("PZ");

  grEffRD->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  grEffRD->GetYaxis()->SetTitle("Trigger Efficiency");

  TLegend *leg1 = new TLegend(0.6694631,0.7947295,0.8691275,0.8942059);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(grEffRD,"Data","PL");
  leg1->AddEntry(grEffMC,"MC","PL");

  leg1->Draw();

  if (savePlots) {
    if(etaRegion == 1){
      c1->cd();
      lt1->DrawLatex(0.5,0.5,"|#eta^{#mu}| < 1.6");
      c1->SaveAs("tnp_trig_eff_mid.png");
      c1->SaveAs("tnp_trig_eff_mid.pdf");
    }
    if(etaRegion == 2){
      c1->cd();
      lt1->DrawLatex(0.615,0.2,"pp  #sqrt{s} = 2.76");
      lt1->DrawLatex(0.5,0.5,"1.6 < |#eta^{#mu}| < 2.4");
      lt1->DrawLatex(0.5,0.44,"6.5 < p_{T}^{J/#psi} < 30 GeV/c");
      //lt1->DrawLatex(0.6,0.44,"6.5 < p_{T} < 30 GeV/c");
      c1->SaveAs("tnp_trig_eff_forward.png");
      c1->SaveAs("tnp_trig_eff_forward.pdf");
    }
  }

  return;
}
