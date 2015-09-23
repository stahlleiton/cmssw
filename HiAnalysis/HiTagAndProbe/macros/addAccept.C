#include "TTree.h"
#include "TFile.h"

TTree* copyTree(TTree* told) {
   TTree *tnew = told->CloneTree(0);
   int tag_accept, accept; float pt,tag_pt,abseta,tag_abseta;
   told->SetBranchAddress("pt",&pt);
   told->SetBranchAddress("tag_pt",&tag_pt);
   told->SetBranchAddress("abseta",&abseta);
   told->SetBranchAddress("tag_abseta",&tag_abseta);
   tnew->Branch("accept",&accept,"accept/I");
   tnew->Branch("tag_accept",&tag_accept,"tag_accept/I");

   int nentries = told->GetEntries();
   for (int i=0; i<nentries; i++) {
      told->GetEntry(i);
      accept = ((abseta<1.2 && pt>3.5) ||(abseta>=1.2 && abseta<2.4 && pt>=4.5-1.25*abseta));
      tag_accept = ((tag_abseta<1.2 && tag_pt>3.5) ||(tag_abseta>=1.2 && tag_abseta<2.4 && tag_pt>=4.5-1.25*tag_abseta));
      tnew->Fill();
   }

   return tnew;
}

void addAccept(const char *filein, const char *fileout) {
   TFile *fin = new TFile(filein);
   TFile *fout = new TFile(fileout,"RECREATE");

   fout->cd();
   TDirectory *tdir_trk = fout->mkdir("MuonTrk");
   tdir_trk->cd();
   TTree *tr_trk = copyTree((TTree*) fin->Get("MuonTrk/fitter_tree"));

   fout->cd();
   TDirectory *tdir_mutrg = fout->mkdir("MuonTrg");
   tdir_mutrg->cd();
   TTree *tr_muidtrg = copyTree((TTree*) fin->Get("MuonTrg/fitter_tree"));

   fout->cd();
   TDirectory *tdir_muidtrg = fout->mkdir("MuonIDTrg");
   tdir_muidtrg->cd();
   TTree *tr_muidtrg = copyTree((TTree*) fin->Get("MuonIDTrg/fitter_tree"));


   fout->Write();
   fout->Close();
}
