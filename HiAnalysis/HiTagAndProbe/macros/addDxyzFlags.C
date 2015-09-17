#include "TTree.h"
#include "TFile.h"

#define DXYCUT 3
#define DZCUT 15

TTree* copyTree(TTree* told) {
   TTree *tnew = told->CloneTree(0);
   int dxyzPVCuts; float dzPV, dxyPVdzmin;
   told->SetBranchAddress("dzPV",&dzPV);
   told->SetBranchAddress("dxyPVdzmin",&dxyPVdzmin);
   tnew->Branch("dxyzPVCuts",&dxyzPVCuts,"dxyzPVCuts/I");

   int nentries = told->GetEntries();
   for (int i=0; i<nentries; i++) {
      told->GetEntry(i);
      dxyzPVCuts = (fabs(dxyPVdzmin)<DXYCUT && fabs(dzPV)<DZCUT);
      tnew->Fill();
   }

   return tnew;
}

void addDxyzFlags(const char *filein, const char *fileout) {
   TFile *fin = new TFile(filein);
   TFile *fout = new TFile(fileout,"RECREATE");

   fout->cd();
   TDirectory *tdir_trk = fout->mkdir("MuonTrk");
   tdir_trk->cd();
   TTree *tr_trk = ((TTree*) fin->Get("MuonTrk/fitter_tree"))->CloneTree();

   fout->cd();
   TDirectory *tdir_muidtrg = fout->mkdir("MuonIDTrg");
   tdir_muidtrg->cd();
   TTree *tr_muidtrg = copyTree((TTree*) fin->Get("MuonIDTrg/fitter_tree"));

   fout->cd();
   TDirectory *tdir_trg = fout->mkdir("MuonTrg");
   tdir_trg->cd();
   TTree *tr_trg = ((TTree*) fin->Get("MuonTrg/fitter_tree"))->CloneTree();


   fout->Write();
   fout->Close();
}
