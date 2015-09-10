#include "TTree.h"
#include "TH1.h"

void optimal_bins(TTree *tr, TString cut, int nmin)
{
   tr->Draw("pt>>hpt(60,0,30)",cut);
   TH1F *hpt = (TH1F*) gDirectory->Get("hpt");
   int nentries = hpt->Integral();
   int nbins = nentries/nmin;
   double *binedges = new double[nbins+1];
   double *bincontent = new double[nbins];
   cout << "Total events: " << nentries << ", " << nmin << " events per bin." << endl;
   int cnt=0;
   int tot=0;
   int ibin=0;
   bool isfirst=true;
   for (int i=1; i<hpt->GetNbinsX()+1; i++)
   {
      cnt += hpt->GetBinContent(i);
      tot += hpt->GetBinContent(i);
      if (isfirst&&cnt>0)
      {
         isfirst=false;
         binedges[0] = hpt->GetBinLowEdge(i);
      }
      if (cnt>= nmin&&(nentries-tot)>=nmin)
      {
         bincontent[ibin] = cnt;
         cnt=0;
         ibin++;
         binedges[ibin] = hpt->GetBinLowEdge(i+1);
         // cout << ibin << " " << binedges[ibin] << " " << bincontent[ibin] << endl;
      }
   }

   if (ibin<nbins)
   {
      bincontent[ibin] = cnt;
      cnt=0;
      ibin++;
      binedges[ibin] = 30.;
      // cout << ibin << "..." << endl;
   }

   for (ibin=0; ibin<nbins; ibin++)
   {
      cout << ibin << ": [" << binedges[ibin] << "," << binedges[ibin+1] << "]: " << bincontent[ibin] << endl;
      if (binedges[ibin+1]==30.) break;
   }
   for (int ibin=0; ibin<nbins; ibin++)
   {
      if (binedges[ibin]==30.) break;
      cout << binedges[ibin] << ", ";
   }
   cout << 30 << endl;
}
