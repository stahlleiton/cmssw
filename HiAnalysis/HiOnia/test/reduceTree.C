#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"

void makeList() {
   TFile *f1 = TFile::Open(
       "/tmp/miheejo/All_Histos_cmssw445p5_RegIt_EvtPlane_noCuts.root"
//       "root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2013/pp/Prompt/TTrees/All_v2.24_Histos_Runs_211739-211831_woPileUpRej_muLessPV.root"
   );
   TTree *myTree = (TTree*) f1->Get("myTree");
   myTree->Draw(">>elist","Reco_QQ_4mom.M()>7.0&&Reco_QQ_4mom.M()<14.0");
//   myTree->Draw(">>elist","Reco_QQ_4mom.M()>2.0&&Reco_QQ_4mom.M()<5.0");
   TEventList *elist = (TEventList*)gDirectory->Get("elist");
   TFile ef("/tmp/miheejo/All_Histos_cmssw445p5_RegIt_EvtPlane_noCuts_elist.root","recreate");
   elist->Write();
   ef.Close();
}

void makeSmall() {
   TFile *f = new TFile("/tmp/miheejo/All_Histos_cmssw445p5_RegIt_EvtPlane_noCuts_elist.root");
   TEventList *elist = (TEventList*)f->Get("elist");
   TFile *f1 = TFile::Open(
       "/tmp/miheejo/All_Histos_cmssw445p5_RegIt_EvtPlane_noCuts.root"
//       "root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2013/pp/Prompt/TTrees/All_v2.24_Histos_Runs_211739-211831_woPileUpRej_muLessPV.root"
   );
   TTree *myTree = (TTree*) f1->Get("myTree");
   myTree->SetEventList(elist);
   TFile *f2 = new TFile("/tmp/miheejo/All_Histos_cmssw445p5_RegIt_EvtPlane_noCuts_M7-14GeV.root","recreate");
   TTree *mySmallTree = myTree->CopyTree("");
   mySmallTree->Write();
   mySmallTree->Print();
   f2->Close();
}

