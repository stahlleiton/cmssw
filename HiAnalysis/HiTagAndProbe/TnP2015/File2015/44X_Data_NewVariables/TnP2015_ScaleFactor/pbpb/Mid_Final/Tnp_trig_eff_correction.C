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

using namespace std;

void Tnp_trig_eff_correction(string dataFile="Jpsi_Regit_RD_2014_Trg2_Eff.root", string mcFile="Jpsi_Regit_MC_2014_Trg2_Eff.root", bool savePlots=true)
{

  int etaRegion = 1; // 1 : mid, 2 : forward

  gROOT->Macro("/Users/dmoon/rootlogon.C");
  gStyle->SetOptFit(0);
  TGraphAsymmErrors *grEffRD;
  TGraphAsymmErrors *grEffMC;
  dataFile="tnp_pbpb_RD_Mid.root";
  mcFile="tnp_pbpb_MC_Mid.root";

  TFile *tnpRD = new TFile(dataFile.c_str(),"READ");
  grEffRD = (TGraphAsymmErrors*) tnpRD->Get("gRD_pT")->Clone("grEffRD");
  //  tnpRD->Close();
  grEffRD->SetMarkerColor(kBlue);
  grEffRD->SetMarkerStyle(21);
  TFile *tnpMC = new TFile(mcFile.c_str(),"READ");
  grEffMC = (TGraphAsymmErrors*) tnpMC->Get("gMC_pT")->Clone("grEffMC");
  //  tnpMC->Close();
  grEffMC->SetMarkerColor(kRed);
  grEffMC->SetMarkerStyle(20);

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",0.0,20);

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();

  f1->SetParNames("eff0","x0","m");
  f1->SetParameters(0.8,0.1,1.0); // 0.8, 0.2, 2.5
  f1->SetLineWidth(2);

  TH1F *hPad = new TH1F("hPad",";p_{T} GeV/c;Trigger x Muon Id Efficiency",10,0,20);
  TH1F *hPad1 = new TH1F("hPad1",";p_{T} GeV/c;Trigger X Muon Id Efficiency",10,0,20);
  cout<<"%%% Fit to Data curve %%%"<<endl;
  grEffRD->Fit(f1,"RME");
  cout<<"%%% Fit to MC curve %%%"<<endl;
  f1->SetParameters(0.9,0.5,2.5);
  grEffMC->Fit(f1,"WRME");


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

  cout<<"Final Correction Function : "<<"("<<fnName1<<")"<<"/"<<"("<<fnName2<<")"<<endl;

  grEffRD->GetFunction("f1")->SetLineColor(kBlue);
  grEffMC->GetFunction("f1")->SetLineColor(kRed);


  TCanvas *c1 = new TCanvas("c1","c1");
  hPad->SetMaximum(1.0);
  hPad->SetMinimum(0.0);
  hPad->Draw();
  grEffRD->Draw("PZ same");
  grEffMC->Draw("PZ same");

  TLegend *leg1 = new TLegend(0.533557,0.4737762,0.7936242,0.5646853);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(grEffRD,"Data","PL");
  leg1->AddEntry(grEffMC,"MC","PL");

  leg1->Draw("same");
  lt1->DrawLatex(0.60,0.22,"PbPb  #sqrt{s_{NN}} = 2.76");

  grEffRD->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  grEffRD->GetYaxis()->SetTitle("Trigger Efficiency");


  if (savePlots) {
    if(etaRegion == 1){
      c1->cd();
      lt1->DrawLatex(0.60,0.16,"|#eta^{#mu}| < 1.6");
      //lt1->DrawLatex(0.16,0.18,"PbPb  #sqrt{s_{NN}} = 2.76");
      //lt1->DrawLatex(0.5,0.5,"|#eta^{#mu}| < 1.6");
      //lt1->DrawLatex(0.5,0.44,"6.5 < p_{T}^{J/#psi} < 30 GeV/c");
      c1->SaveAs("tnp_trig_eff_mid_pbpb.png");
      c1->SaveAs("tnp_trig_eff_mid_pbpb.pdf");
    }
    if(etaRegion == 2){
      c1->cd();
      lt1->DrawLatex(0.65,0.16,"1.6 < |#eta^{#mu}| < 2.4");
      lt1->DrawLatex(0.16,0.18,"PbPb  #sqrt{s_{NN}} = 2.76");
      lt1->DrawLatex(0.6,0.5,"1.6 < |#eta^{#mu}| < 2.4");
      lt1->DrawLatex(0.6,0.44,"6.5 < p_{T}^{J/#psi} < 30 GeV/c");
      c1->SaveAs("tnp_trig_eff_forward_pbpb.png");
      c1->SaveAs("tnp_trig_eff_forward_pbpb.pdf");
    }
  }

  return;
}
