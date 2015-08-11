// This is macro to get 100 scale factor functions for systematic uncertainty study using Toy model
// Editor : Dong Ho Moon, if you have problem, please send e-mail to dmoon@cern.ch

#include <iostream>
#include <string>

#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"

using namespace std;

double assym_gaus(double *x, double *par)
{
  double mean = par[0];
  double sigmah = par[1];
  double sigmal = par[2];

  if ((x[0] < mean))
    return TMath::Exp(-TMath::Power((x[0]-mean)/par[1],2));
  else
    return TMath::Exp(-TMath::Power((x[0]-mean)/par[2],2));
}


void Tnp_trig_eff_correction_ToyModel(string dataFile="Jpsi_Regit_RD_2014_Trg2_Eff.root", string mcFile="Jpsi_Regit_MC_2014_Trg2_Eff.root", bool savePlots=true)
{

  int etaRegion = 2; // 1 : mid, 2 : forward

  gROOT->Macro("~/rootlogon.C");
  gStyle->SetOptFit(0);
  TGraphAsymmErrors *grEffRD;
  TGraphAsymmErrors *grEffMC;
  dataFile="tnp_pbpb_RD_Forw.root";
  mcFile="tnp_pbpb_MC_Forw.root";

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

  TH1F *hPad = new TH1F("hPad",";p_{T} GeV/c;Trigger X Muon Id Efficiency",10,0,20);
  TH1F *hPad1 = new TH1F("hPad1",";p_{T} GeV/c;Trigger X Muon Id Efficiency",10,0,20);
  grEffRD->Fit(f1,"RME");
  f1->SetParameters(0.9,0.5,2.5);
  grEffMC->Fit(f1,"WRME");

  grEffRD->GetFunction("f1")->SetLineColor(kBlue);
  grEffMC->GetFunction("f1")->SetLineColor(kRed);


  char fnName2[512];
  sprintf(fnName2,"%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",grEffMC->GetFunction("f1")->GetParameter(0),grEffMC->GetFunction("f1")->GetParameter(1),grEffMC->GetFunction("f1")->GetParameter(2));
  //cout<<fnName2<<endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  hPad->Draw();
  grEffRD->Draw("PZ");
  grEffMC->Draw("PZ");

  TLegend *leg1 = new TLegend(0.533557,0.4737762,0.7936242,0.5646853);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(grEffRD,"Data","PL");
  leg1->AddEntry(grEffMC,"MC","PL");

  leg1->Draw("same");
  lt1->DrawLatex(0.16,0.9,"PbPb  #sqrt{s_{NN}} = 2.76");

  grEffRD->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  grEffRD->GetYaxis()->SetTitle("Trigger Efficiency");
  
  // Study with Toy Model


  const int nBin = 4;
  double effRD[nBin] = {0.0};
  double errRDh[nBin] = {0.0};
  double errRDl[nBin] = {0.0};
  double xaxisRD[nBin] = {0.0};
  double effMC[nBin] = {0.0};
  double errMCh[nBin] = {0.0};
  double errMCl[nBin] = {0.0};
  double xaxisMC[nBin] = {0.0};
  for(int i = 0; i < nBin; i++){
    grEffRD->GetPoint(i, xaxisRD[i], effRD[i]);
    errRDh[i] = grEffRD->GetErrorYhigh(i);
    errRDl[i] = grEffRD->GetErrorYlow(i);
    cout<<"Data : "<<xaxisRD[i]<<", "<<effRD[i]<<", "<<errRDh[i]<<", "<<errRDl[i]<<endl;

    grEffMC->GetPoint(i, xaxisMC[i], effMC[i]);
    errMCh[i] = grEffMC->GetErrorYhigh(i);
    errMCl[i] = grEffMC->GetErrorYlow(i);
    cout<<"MC : "<<xaxisMC[i]<<", "<<effMC[i]<<", "<<errMCh[i]<<", "<<errMCl[i]<<endl;
  }

  

  TF1 *fEffRDPt[nBin];
  char tmp[512];
  for(int i = 0; i < nBin; i++){
    sprintf(tmp,"fEffRDPt_%d",i);
    //cout<<"name of fn : "<<tmp<<endl;
    fEffRDPt[i] = new TF1(tmp,assym_gaus,0,1.0,3);
    ////fEffRDPt[i] = new TF1(tmp,"[0]*exp(-(x-[1])*(x-[1])/[2]*[2])",0,1.0);
    //fEffRDPt[i] = new TF1(tmp,"gaus",0,1.0);
    fEffRDPt[i]->SetParameter(0, effRD[i]);
    fEffRDPt[i]->SetParameter(1, errRDl[i]);
    fEffRDPt[i]->SetParameter(2, errRDh[i]);
  }

  TCanvas *cAs = new TCanvas();
  cAs->cd();
  TF1 *tfAsyGauss = new TF1("tfAsy",assym_gaus,0,1.0,3);
  tfAsyGauss->SetParameter(0, 0.5);
  tfAsyGauss->SetParameter(1, 0.08);
  tfAsyGauss->SetParameter(2, 0.12);
  //tfAsyGauss->SetParameter(0, effRD[0]);
  //tfAsyGauss->SetParameter(1, errRDh[0]);
  //tfAsyGauss->SetParameter(2, errRDl[0]);
  tfAsyGauss->Draw();
  cAs->SaveAs("tmp_gauss_asym.png");
  cAs->SaveAs("tmp_gauss_asym.pdf");


  // pN = number of toy

  TCanvas *c2 = new TCanvas();
  const int pN = 100; 
  double effToy[pN][nBin] ={{0.0}};
  double xErrh[nBin] = {0.0};
  double xErrl[nBin] = {0.0};
  TH1F *h1 = new TH1F("h1","",100,0,1.0);
  for(int itoy = 0; itoy < pN; itoy++){
    for(int i = 0; i < nBin; i++){
      effToy[itoy][i]= fEffRDPt[i]->GetRandom();
      //cout<<"eff Toy : "<<effToy[itoy][i]<<endl;
      //if(i == 3) h1->Fill(effToy[itoy][i]);
    }
  }
  //h1->Draw();
  char tmp1[512];
  TGraphAsymmErrors *grEffRDToy[pN];
  for(int itoy = 0; itoy < pN; itoy++){
    grEffRDToy[itoy] = new TGraphAsymmErrors(nBin, xaxisRD, effToy[itoy], xErrh, xErrl, errRDh, errRDl);
    hPad->Draw();
    grEffRDToy[itoy]->Fit(f1,"RME");
    grEffRDToy[itoy]->Draw("pz");
    sprintf(tmp1,"figs_eff/eff_toy_%d.pdf",itoy);
    c2->SaveAs(tmp1);
  }

  TCanvas *c3 = new TCanvas("c3","",1670,990);
  c3->Divide(10,10);
  char fnName1[512];
  double pr1[pN] = {0.0};
  double pr2[pN] = {0.0};
  double pr3[pN] = {0.0};
  for(int i = 0; i < pN; i++){
    c3->cd(i+1); 
    hPad->Draw(); 
    grEffRDToy[i]->Fit(f1,"RME");
    pr1[i]=grEffRDToy[i]->GetFunction("f1")->GetParameter(0);
    pr2[i]=grEffRDToy[i]->GetFunction("f1")->GetParameter(1);
    pr3[i]=grEffRDToy[i]->GetFunction("f1")->GetParameter(2);
    //sprintf(fnName1,"%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",pr1,pr2,pr3);
    //cout<<fnName1<<endl;
    //if(i == 5 || i == 43) cout<<"Final Correction Function : "<<"("<<fnName1<<")"<<"/"<<"("<<fnName2<<")"<<endl;
    grEffRDToy[i]->Draw("pz");
  }
  c3->SaveAs("eff_comp.pdf");

  TCanvas *c4 = new TCanvas("c3","",660,600);


  /*
  for(int i = 0; i < pN+1; i++){
    if(i == 0){c3->cd(i+1); hPad->Draw(); grEffRD->Draw("pz");}
    if(i != 0){
      c3->cd(i+1); 
      hPad->Draw(); 
      grEffRDToy[i-1]->Fit(f1,"RME");
      pr1[i-1]=grEffRDToy[i-1]->GetFunction("f1")->GetParameter(0);
      pr2[i-1]=grEffRDToy[i-1]->GetFunction("f1")->GetParameter(1);
      pr3[i-1]=grEffRDToy[i-1]->GetFunction("f1")->GetParameter(2);
      //sprintf(fnName1,"%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",pr1,pr2,pr3);
      //cout<<fnName1<<endl;
      //cout<<"Final Correction Function : "<<"("<<fnName1<<")"<<"/"<<"("<<fnName2<<")"<<endl;
      grEffRDToy[i-1]->Draw("pz");
    }
  }
  c3->SaveAs("eff_comp.pdf");
  */

  for(int i = 0; i < pN; i++){
    sprintf(fnName1,"%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",pr1[i],pr2[i],pr3[i]);
    //cout<<fnName1<<endl;
    cout<<"Final Correction Function : "<<"("<<fnName1<<")"<<"/"<<"("<<fnName2<<")"<<endl;
  }


  if (savePlots) {
    if(etaRegion == 1){
      c1->cd();
      lt1->DrawLatex(0.16,0.85,"|#eta^{#mu}| < 1.6");
      c1->SaveAs("tnp_trig_eff_mid_pbpb.png");
      c1->SaveAs("tnp_trig_eff_mid_pbpb.pdf");
    }
    if(etaRegion == 2){
      c1->cd();
      lt1->DrawLatex(0.16,0.84,"1.6 < |#eta^{#mu}| < 2.4");
      c1->SaveAs("tnp_trig_eff_forward_pbpb.png");
      c1->SaveAs("tnp_trig_eff_forward_pbpb.pdf");
    }
  }

  /*
  TF1 *test = new TF1("test",assym_gaus,0,1.0,3);
  test->SetParameter(0, effRD[0]);
  test->SetParameter(1, errRDh[0]);
  test->SetParameter(2, errRDl[0]);
  TCanvas *cc2 = new TCanvas("cc2","",500,550);
  cc2->cd();
  test->Draw();
  */
  //return;
}
