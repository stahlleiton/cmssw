#ifndef toyStudy_h
#define toyStudy_h

#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TArrow.h"

#include <fstream>

const int nToys = 100;

double assym_gaus(double *x, double *par) {
   double mean = par[0];
   double sigmal = par[1];
   double sigmah = par[2];

   // beware, sometimes one of the two parameters is 0... that's bad
   if (sigmal<1e-3) sigmal = max(1e-3, sigmah);
   if (sigmah<1e-3) sigmah = max(1e-3, sigmal);

   if ((x[0] < mean))
      return TMath::Gaus(x[0],mean,sigmal,kTRUE);
   else
      return TMath::Gaus(x[0],mean,sigmah,kTRUE);
};

double* systerr(TGraphAsymmErrors **graph, int nsyst) {
   int nbins = graph[0]->GetN();
   double *ans = new double[nbins];

   for (int i=0; i<nbins; i++) {
      double ynom = graph[0]->GetY()[i];
      double eylnom = graph[0]->GetEYlow()[i];
      double eyhnom = graph[0]->GetEYhigh()[i];
      ans[i] = 0;
      for (int k=1; k<nsyst; k++) {
         double ysyst = graph[k]->GetY()[i];
         double eylsyst = graph[k]->GetEYlow()[i];
         double eyhsyst = graph[k]->GetEYhigh()[i];
         if (ysyst-ynom>0) { // upward fluctuation
            double errc = sqrt(fabs(pow(eyhsyst,2)-pow(eyhnom,2)));
            if (ysyst-ynom>errc && ysyst-ynom>ans[i]) ans[i] = ysyst-ynom;
            // if (ysyst-ynom<errc && ysyst-ynom>1e-2) cout << "AAAA " << i << " " << k << " " << ysyst-ynom << " " << errc << endl;
         }
         if (ysyst-ynom<0) { // downward fluctuation
            double errc = sqrt(fabs(pow(eylsyst,2)-pow(eylnom,2)));
            if (-ysyst+ynom>errc && -ysyst+ynom>ans[i]) ans[i] = -ysyst+ynom;
            // if (-ysyst+ynom<errc && -ysyst+ynom>1e-2) cout << "AAAA " << i << " " << k << " " << ysyst-ynom << " " << errc << endl;
         }
      }
   }

   return ans;
}

void toyStudy(int nsyst, TGraphAsymmErrors **gdata, TGraphAsymmErrors **gmc, TF1 *fdata, TF1 *fmc, TString outputname, int systmode=0) {
   const int nbinsdata = gdata[0]->GetN();
   const double *xdata = gdata[0]->GetX();
   const double *exldata = gdata[0]->GetEXlow();
   const double *exhdata = gdata[0]->GetEXhigh();
   const double *ydata = gdata[0]->GetY();
   const double *eyldata = gdata[0]->GetEYlow();
   const double *eyhdata = gdata[0]->GetEYhigh();
   double *eyldata_syst = systerr(gdata,nsyst);
   double *eyhdata_syst = eyldata_syst;
   const int nbinsmc = gmc[0]->GetN();
   const double *xmc = gmc[0]->GetX();
   const double *exlmc = gmc[0]->GetEXlow();
   const double *exhmc = gmc[0]->GetEXhigh();
   const double *ymc = gmc[0]->GetY();
   const double *eylmc = gmc[0]->GetEYlow();
   const double *eyhmc = gmc[0]->GetEYhigh();
   double *eylmc_syst = systerr(gmc,nsyst);
   double *eyhmc_syst = eyldata_syst;

   TF1 *tfAsyGauss = new TF1("tfAsy",assym_gaus,0,1.0,3);
   // tfAsyGauss->SetNpx(1000);

   ofstream fileout(TString("correction_functions_toys_") + outputname + ".txt");
   TFile *tfileout = new TFile(TString("correction_functions_toys_") + outputname + ".root","RECREATE");
   tfileout->cd();

   float *effdata = new float[nbinsdata];
   float *effdata_fit = new float[nbinsdata];
   float *effmc = new float[nbinsmc];
   float *effmc_fit = new float[nbinsmc];

   TTree *tr = new TTree("tr","tr");
   tr->Branch("effdata",effdata,Form("effdata[%i]/F",nbinsdata));
   tr->Branch("effdata_fit",effdata_fit,Form("effdata_fit[%i]/F",nbinsdata));
   tr->Branch("effmc",effmc,Form("effmc[%i]/F",nbinsmc));
   tr->Branch("effmc_fit",effmc_fit,Form("effmc_fit[%i]/F",nbinsmc));

   TCanvas *cdata = new TCanvas("cdata");
   TCanvas *cmc = new TCanvas("cmc");
   TH1F *hPad = new TH1F("hPad",";p^{#mu}_{T} [GeV/c];Single #mu Efficiency",5,0,30);
   hPad->GetXaxis()->CenterTitle();
   hPad->GetXaxis()->SetLabelSize(0.05);
   hPad->GetXaxis()->SetTitleSize(0.05);
   hPad->GetXaxis()->SetTitleOffset(1.2);
   hPad->GetYaxis()->SetLabelSize(0.05);
   hPad->GetYaxis()->SetTitleSize(0.05);
   hPad->GetYaxis()->SetTitleOffset(1.);
   hPad->GetYaxis()->SetRangeUser(0.,1.05);
   cdata->cd(); hPad->Draw();
   cmc->cd(); hPad->Draw();

   int systmodemin=systmode;
   int systmodemax=systmode;
   if (systmode==2) {
      systmodemin=0;
      systmodemax=1;
   }

   for (int systmode2=systmodemax; systmode2>=systmodemin; systmode2--) {
      for (int i=0; i<nToys; i++) {
         // generate the toys
         // data
         TGraphAsymmErrors *gtoydata = new TGraphAsymmErrors(nbinsdata,xdata,ydata,exldata,exhdata,eyldata,eyhdata);
         for (int j=0; j<nbinsdata; j++) {
            if (systmode2==0) tfAsyGauss->SetParameters(ydata[j],eyldata[j],eyhdata[j]);
            else tfAsyGauss->SetParameters(ydata[j],sqrt(pow(eyldata[j],2)+pow(eyldata_syst[j],2)),sqrt(pow(eyhdata[j],2)+pow(eyhdata_syst[j],2)));
            effdata[j] = tfAsyGauss->GetRandom();
            gtoydata->SetPoint(j,xdata[j],effdata[j]);
         }

         // mc
         TGraphAsymmErrors *gtoymc = new TGraphAsymmErrors(nbinsmc,xmc,ymc,exlmc,exhmc,eylmc,eyhmc);
         for (int j=0; j<nbinsmc; j++) {
            if (systmode2==0) tfAsyGauss->SetParameters(ymc[j],eylmc[j],eyhmc[j]);
            else tfAsyGauss->SetParameters(ymc[j],sqrt(pow(eylmc[j],2)+pow(eylmc_syst[j],2)),sqrt(pow(eyhmc[j],2)+pow(eyhmc_syst[j],2)));
            effmc[j] = tfAsyGauss->GetRandom();
            gtoymc->SetPoint(j,xmc[j],effmc[j]);
         }

         // fit the graphs
         TF1 *ftoydata = (TF1*) fdata->Clone(Form("ftoydata%i",i));
         ftoydata->SetLineColor(systmode2==1 ? kOrange+3 : kOrange+7);
         // ftoydata->SetParameter(2,gtoydata->GetX()[nbinsdata]);
         TF1 *ftoymc = (TF1*) fmc->Clone(Form("ftoymc%i",i));
         ftoymc->SetLineColor(systmode2==1 ? kGray : kBlack);
         // ftoymc->SetParameter(2,gtoymc->GetX()[nbinsmc]);

         gtoydata->Fit(ftoydata,"WRM");
         gtoymc->Fit(ftoymc,"WRM");

         // compute the efficiency from fit, for the tree
         for (int j=0; j<nbinsdata; j++) {
            effdata_fit[j] = ftoydata->Eval(xdata[j]);
         }
         for (int j=0; j<nbinsmc; j++) {
            effmc_fit[j] = ftoymc->Eval(xmc[j]);
         }

         // draw
         cdata->cd();
         ftoydata->Draw("l same");
         cmc->cd();
         ftoymc->Draw("l same");

         // print results to file
         if (systmode2==systmode) {
            fileout << Form("%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",ftoydata->GetParameter(0),ftoydata->GetParameter(1),ftoydata->GetParameter(2)) << endl;
            // fileout << Form("%0.4f*TMath::Erf((x-%0.4f)/%0.4f)",ftoymc->GetParameter(0),ftoymc->GetParameter(1),ftoymc->GetParameter(2)) << endl;
            tr->Fill();
         }
      }
   }


   // superimpose nominal results
   cdata->cd();
   double eyldatastatsyst[nbinsdata];
   double eyhdatastatsyst[nbinsdata];
   for (int i=0; i<nbinsdata; i++) {
      eyldatastatsyst[i] = sqrt(pow(eyldata[i],2)+pow(eyldata_syst[i],2));
      eyhdatastatsyst[i] = sqrt(pow(eyhdata[i],2)+pow(eyhdata_syst[i],2));
   }
   TGraphAsymmErrors *gdatasyst = new TGraphAsymmErrors(nbinsdata,xdata,ydata,exldata,exhdata,eyldatastatsyst,eyhdatastatsyst);
   gdatasyst->Draw("pz same");
   gdata[0]->Draw("p|| same");
   fdata->Draw("l same");
   TLegend *tlegdata = new TLegend(0.6,0.2,0.9,0.4);
   tlegdata->AddEntry(gdata[0],"Data","lp");
   tlegdata->AddEntry(fdata,"Nominal fit","lp");
   tlegdata->Draw();
   cmc->cd();
   double eylmcstatsyst[nbinsmc];
   double eyhmcstatsyst[nbinsmc];
   for (int i=0; i<nbinsmc; i++) {
      eylmcstatsyst[i] = sqrt(pow(eylmc[i],2)+pow(eylmc_syst[i],2));
      eyhmcstatsyst[i] = sqrt(pow(eyhmc[i],2)+pow(eyhmc_syst[i],2));
   }
   TGraphAsymmErrors *gmcsyst = new TGraphAsymmErrors(nbinsmc,xmc,ymc,exlmc,exhmc,eylmcstatsyst,eyhmcstatsyst);
   gmcsyst->Draw("pz same");
   gmc[0]->Draw("p|| same");
   fmc->Draw("l same");
   TLegend *tlegmc = new TLegend(0.6,0.2,0.9,0.4);
   tlegmc->AddEntry(gmc[0],"MC","lp");
   tlegmc->AddEntry(fmc,"Nominal fit","lp");
   tlegmc->Draw();

   cdata->SaveAs(TString("toysdata_") + outputname + ".pdf");
   cmc->SaveAs(TString("toysmc_") + outputname + ".pdf");

   // plot the efficiency distribution for each bin
   TCanvas *cdata_dists = new TCanvas("cdata_dists","cdata_dists",1670,560);
   cdata_dists->Divide(nbinsdata/2,2);
   for (int j=0; j<nbinsdata; j++) {
      cdata_dists->cd(j+1);
      TString hname = Form("distdata%i",j);
      TString hname_fit = Form("distdata%i_fit",j);
      double xmin = min(ydata[j]-5.*eyldata[j], fdata->Eval(xdata[j])-5.*eyldata[j]);
      double xmax = max(ydata[j]+5.*eyldata[j], fdata->Eval(xdata[j])+5.*eyldata[j]);
      TH1F *hdist = new TH1F(hname,hname,10,xmin,xmax);
      hdist->GetXaxis()->SetTitle("Efficiency"); hdist->GetYaxis()->SetTitle("Counts");
      hdist->SetLineColor(kBlue); hdist->SetMarkerColor(kBlue);
      TH1F *hdist_fit = new TH1F(hname_fit,hname_fit,10,ydata[j]-5.*eyldata[j],ydata[j]+5.*eyhdata[j]);
      hdist_fit->GetXaxis()->SetTitle("Efficiency"); hdist_fit->GetYaxis()->SetTitle("Counts");
      hdist_fit->SetLineColor(kCyan); hdist_fit->SetMarkerColor(kCyan);
      tr->Draw(Form("effdata[%i]>>",j) + hname);
      tr->Draw(Form("effdata_fit[%i]>>",j) + hname_fit);
      hdist->Draw(); hdist->GetYaxis()->SetRangeUser(0,100); 
      hdist_fit->Draw("same");
      if (j==0) {
         TLegend *tleg = new TLegend(0.7,0.7,0.9,0.9);
         tleg->AddEntry(hdist,"Toys","lp");
         tleg->AddEntry(hdist_fit,"Fit","lp");
         tleg->Draw();
      }
      TLatex *lt1 = new TLatex();
      lt1->SetNDC();
      lt1->DrawLatex(0.43,0.95,Form("p_{T} #in [%.1f,%.1f]",xdata[j]-exldata[j],xdata[j]+exhdata[j]));

      // represent the input for these distribitions
      TArrow *arrlow = new TArrow(ydata[j],70.,ydata[j]-eyldata[j],70.,0.01,"|->");
      arrlow->SetLineColor(kBlue); arrlow->Draw();
      TArrow *arrhigh = new TArrow(ydata[j],70.,ydata[j]+eyhdata[j],70.,0.01,"|->");
      arrhigh->SetLineColor(kBlue); arrhigh->Draw();
      TArrow *arrfit = new TArrow(fdata->Eval(xdata[j]),50.,fdata->Eval(xdata[j]),0.,0.01,"->");
      arrfit->SetLineColor(kCyan); arrfit->Draw();
   }
   cdata_dists->SaveAs(TString("toysdata_dists_") + outputname + ".pdf");

   TCanvas *cmc_dists = new TCanvas("cmc_dists","cmc_dists",1670,560);
   cmc_dists->Divide(nbinsmc/2,2);
   for (int j=0; j<nbinsmc; j++) {
      cmc_dists->cd(j+1);
      TString hname = Form("distmc%i",j);
      TString hname_fit = Form("distmc%i_fit",j);
      double xmin = min(ymc[j]-5.*eylmc[j], fmc->Eval(xmc[j])-5.*eylmc[j]);
      double xmax = max(ymc[j]+5.*eylmc[j], fmc->Eval(xmc[j])+5.*eylmc[j]);
      TH1F *hdist = new TH1F(hname,hname,10,xmin,xmax);
      hdist->GetXaxis()->SetTitle("Efficiency"); hdist->GetYaxis()->SetTitle("Counts");
      hdist->SetLineColor(kRed); hdist->SetMarkerColor(kRed);
      TH1F *hdist_fit = new TH1F(hname_fit,hname_fit,10,ymc[j]-5.*eylmc[j],ymc[j]+5.*eyhmc[j]);
      hdist_fit->GetXaxis()->SetTitle("Efficiency"); hdist_fit->GetYaxis()->SetTitle("Counts");
      hdist_fit->SetLineColor(kMagenta); hdist_fit->SetMarkerColor(kMagenta);
      tr->Draw(Form("effmc[%i]>>",j) + hname);
      tr->Draw(Form("effmc_fit[%i]>>",j) + hname_fit);
      hdist->Draw(); hdist->GetYaxis()->SetRangeUser(0,100);
      hdist_fit->Draw("same");
      if (j==0) {
         TLegend *tleg = new TLegend(0.7,0.7,0.9,0.9);
         tleg->AddEntry(hdist,"Toys","lp");
         tleg->AddEntry(hdist_fit,"Fit","lp");
         tleg->Draw();
      }
      TLatex *lt1 = new TLatex();
      lt1->SetNDC();
      lt1->DrawLatex(0.43,0.95,Form("p_{T} #in [%.1f,%.1f]",xmc[j]-exlmc[j],xmc[j]+exhmc[j]));

      // represent the input for these distribitions
      TArrow *arrlow = new TArrow(ymc[j],70.,ymc[j]-eylmc[j],70.,0.01,"|->");
      arrlow->SetLineColor(kRed); arrlow->Draw();
      TArrow *arrhigh = new TArrow(ymc[j],70.,ymc[j]+eyhmc[j],70.,0.01,"|->");
      arrhigh->SetLineColor(kRed); arrhigh->Draw();
      TArrow *arrfit = new TArrow(fmc->Eval(xmc[j]),50.,fmc->Eval(xmc[j]),0.,0.01,"->");
      arrfit->SetLineColor(kMagenta); arrfit->Draw();
   }
   cmc_dists->SaveAs(TString("toysmc_dists_") + outputname + ".pdf");

   // write and close open files
   fileout.close();
   tfileout->Write();
   tfileout->Close();

   // more cleaning
   delete cdata, cmc;
   delete cdata_dists, cmc_dists;
};

#endif // #ifndef toyStudy_h
