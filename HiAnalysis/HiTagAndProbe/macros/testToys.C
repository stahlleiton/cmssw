#include "toyStudy.h"

#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TString.h"


void testToys() {
   TF1 *fdata = new TF1("fdata","[0]*TMath::Erf((x-[1])/[2])",1.5,30.);
   fdata->SetParameters(0.95,1.,2.);
   fdata->SetLineWidth(2);
   fdata->SetLineColor(kBlue);

   const int npoints = 6;
   double x[6] = {2.75,3.25,3.75,4.25,5,17.75};
   double ex[6] = {0.25,0.25,0.25,0.25,0.5,12.25};
   double y[6];
   double eydata[6]; double eymc[6];
   for (int i=0; i<6; i++) {
      y[i] = fdata->Eval(x[i]);
      eydata[i] = 0.05;
      eymc[i] = 0.005;
   }

   TGraphAsymmErrors *gdata = new TGraphAsymmErrors(npoints,x,y,ex,ex,eydata,eydata);
   TGraphAsymmErrors *gmc = new TGraphAsymmErrors(npoints,x,y,ex,ex,eymc,eymc);

   gmc->SetMarkerStyle(20);
   gmc->SetMarkerSize(1.4);
   gmc->SetMarkerColor(kRed+1);
   gmc->SetLineColor(kRed+1);
   gdata->SetMarkerStyle(25);
   gdata->SetMarkerSize(1.4);
   gdata->SetMarkerColor(kBlue+1);
   gdata->SetLineColor(kBlue+1);

   TF1 *fmc = (TF1*) fdata->Clone("fmc");
   fmc->SetLineColor(kRed);

   toyStudy(gdata, gmc, fdata, fmc, "testToys");
}
