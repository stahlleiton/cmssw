#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TKey.h"

#include <iostream>

using namespace std;

// A "simple" macro to read a ROOT file as created by the TnP fitting code, and print the canvases corresponding to the mass fits,
// as well as the final efficiency plot.
// The code is quite ugly, sorry about that... The idea was to iterate on the contents on the file to guess where are the interesting stuff.

void plot(const char *filename, const char *var)
{
   TFile *f = new TFile(filename);
   TIter next(f->GetListOfKeys()); TObject *obj;
   while ((obj = next()))
   {
      obj = ((TKey*) obj)->ReadObj();
      if (TString(obj->ClassName()) == "TDirectoryFile")
      {
         TDirectoryFile *tdir = dynamic_cast<TDirectoryFile*>(obj);
         tdir->cd();
         TIter next2(tdir->GetListOfKeys()); TObject *obj2;
         while ((obj2 = next2()))
         {
            obj2 = ((TKey*) obj2)->ReadObj();
            if (TString(obj2->ClassName()) == "TDirectoryFile")
            {
               TDirectoryFile *tdir2 = dynamic_cast<TDirectoryFile*>(obj2);
               tdir2->cd();
               TIter next3(tdir2->GetListOfKeys()); TObject *obj3;
               while ((obj3 = next3()))
               {
                  obj3 = ((TKey*) obj3)->ReadObj();
                  if (TString(obj3->ClassName()) == "TDirectoryFile" && TString(obj3->GetName()).Contains("bin"))
                  {
                     TDirectoryFile *tdir3 = dynamic_cast<TDirectoryFile*>(obj3);
                     tdir3->cd();
                     TCanvas *canv = (TCanvas*) gDirectory->Get("fit_canvas");
                     canv->SaveAs(TString(obj3->GetName()) + TString(".pdf"));
                     canv->SaveAs(TString(obj3->GetName()) + TString(".png"));
                  }
                  if (TString(obj3->GetTitle()) == "fit_eff_plots")
                  {
                     TDirectoryFile *tdir3 = dynamic_cast<TDirectoryFile*>(obj3);
                     tdir3->cd();
                     TIter next4(tdir3->GetListOfKeys()); TObject *obj4;
                     while ((obj4 = next4()))
                     {
                        obj4 = ((TKey*) obj4)->ReadObj();
                        if (TString(obj4->GetName()).Contains(TString(var)+TString("_PLOT"))
                              && !TString(obj4->GetName()).Contains(TString("_") + TString(var)) )
                        {
                           cout << obj4->GetName() << endl;
                           TCanvas *canv = (TCanvas*) obj4;
                           canv->SaveAs(TString(var) + TString(".pdf"));
                           canv->SaveAs(TString(var) + TString(".png"));
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
