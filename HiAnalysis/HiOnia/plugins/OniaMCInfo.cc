#include "HiAnalysis/HiOnia/interface/HiOniaAnalyzer.h"

reco::GenParticleRef HiOniaAnalyzer::findMotherRef(reco::GenParticleRef GenParticleMother, int GenParticlePDG){
  for(int i=0; i<1000; ++i) {
    if (GenParticleMother.isNonnull() && (GenParticleMother->pdgId()==GenParticlePDG) && GenParticleMother->numberOfMothers()>0) {        
      GenParticleMother = GenParticleMother->motherRef();
    } else break;
  }
  return GenParticleMother;
};

std::vector<reco::GenParticleRef> HiOniaAnalyzer::GenBrothers(reco::GenParticleRef GenParticleMother, int GenJpsiPDG){
  bool foundJpsi = false;
  std::vector<reco::GenParticleRef> res;

  if(!GenParticleMother.isNonnull()) return res;
  //if(Reco_3mu_size>0) cout<<"\nScanning daughters of Jpsi "<<GenJpsiPDG<<" mother, pdg = "<<GenParticleMother->pdgId()<<endl; 
  for(int i=0;i<(int)GenParticleMother->numberOfDaughters();i++){    
    reco::GenParticleRef dau = findDaughterRef(GenParticleMother->daughterRef(i), GenParticleMother->pdgId());
    for(int l=0;l<100;l++){ //avoid having a daughter of same pdgId
      if(!(dau.isNonnull() && dau->status()>0 && dau->status()<1000)) break;
      if(dau->pdgId()==GenParticleMother->pdgId() && dau->numberOfDaughters()==1)
	dau = findDaughterRef(dau->daughterRef(0), dau->pdgId());
      else break;
    }

    if(!(dau.isNonnull() && dau->status()>0 && dau->status()<1000) ) continue;
    //if(Reco_3mu_size>0) cout<<"Daughter #"<<i<<" pdg = "<< dau->pdgId()<<" pt,eta = "<<dau->pt()<<" "<<dau->eta()<<endl;
    if(isChargedTrack(dau->pdgId())){
      res.push_back(dau);}
    if(dau->pdgId()==GenJpsiPDG) {
      foundJpsi = true; //continue;
    }

    for(int j=0;j<(int)dau->numberOfDaughters();j++){
      reco::GenParticleRef grandDau = findDaughterRef(dau->daughterRef(j), dau->pdgId());
      for(int l=0;l<100;l++){ //avoid having a daughter of same pdgId
	if(!(grandDau.isNonnull() && grandDau->status()>0 && grandDau->status()<1000)) break;
	if(grandDau->pdgId()==dau->pdgId() && grandDau->numberOfDaughters()==1)
	  grandDau = findDaughterRef(grandDau->daughterRef(0), grandDau->pdgId());
	else break;
      }

      if(!(grandDau.isNonnull() && grandDau->status()>0 && grandDau->status()<1000)) continue;
      //if(Reco_3mu_size>0) cout<<"    grand-daughter #"<<j<<" pdg = "<< grandDau->pdgId()<<" pt,eta = "<<grandDau->pt()<<" "<<grandDau->eta()<<endl;
      if(isChargedTrack(grandDau->pdgId())){
	res.push_back(grandDau);}
      if(grandDau->pdgId()==GenJpsiPDG) {
	foundJpsi = true; //continue;
      }

      for(int k=0;k<(int)grandDau->numberOfDaughters();k++){
	reco::GenParticleRef ggrandDau = findDaughterRef(grandDau->daughterRef(k), grandDau->pdgId());
	for(int l=0;l<100;l++){ //avoid having a daughter of same pdgId
	  if(!(ggrandDau.isNonnull() && ggrandDau->status()>0 && ggrandDau->status()<1000)) break;
	  if(ggrandDau->pdgId()==grandDau->pdgId() && ggrandDau->numberOfDaughters()==1)
	    ggrandDau = findDaughterRef(ggrandDau->daughterRef(0), ggrandDau->pdgId());
	  else break;
	}

	if(!(ggrandDau.isNonnull() && ggrandDau->status()>0 && ggrandDau->status()<1000)) continue;
	//if(Reco_3mu_size>0) cout<<"        grand-grand-daughter #"<<k<<" pdg = "<< ggrandDau->pdgId()<<" pt,eta = "<<ggrandDau->pt()<<" "<<ggrandDau->eta()<<endl;
	if(isChargedTrack(ggrandDau->pdgId())){
	  res.push_back(ggrandDau);}
	if(ggrandDau->pdgId()==GenJpsiPDG) {
	  foundJpsi = true; //continue;
	}
      }
    }
  }

  if(!foundJpsi){
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Incoherence in genealogy: Jpsi not found in the daughters!\n"<<endl;
  }
  // if(!isAbHadron(GenParticleMother->pdgId())){
  //   cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!! Jpsi ancestor is not a b-hadron! pdgID(mother of this ancestor) = "<<findMotherRef(GenParticleMother->motherRef() , GenParticleMother->pdgId())->pdgId()<<endl;
  // }  

  return res;
};

std::pair< std::vector<reco::GenParticleRef>, std::pair<float, float> > HiOniaAnalyzer::findGenMCInfo(const reco::GenParticle* genJpsi){

  float trueLife = -99.;
  float trueLife3D = -99.;
  std::vector<reco::GenParticleRef> JpsiBrothers;

  if (genJpsi->numberOfMothers()>0) {
    TVector3 trueVtx(0.0,0.0,0.0);
    TVector3 trueP(0.0,0.0,0.0);
    TVector3 trueVtxMom(0.0,0.0,0.0);

    trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
    trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

    bool aBhadron = false;
    reco::GenParticleRef Jpsimom_final;
    reco::GenParticleRef Jpsimom = findMotherRef(genJpsi->motherRef(), genJpsi->pdgId());      

    if (Jpsimom.isNull()) {
      std::pair<float, float> trueLifePair = std::make_pair(trueLife, trueLife3D);
      std::pair< std::vector<reco::GenParticleRef>, std::pair<float, float> > result = std::make_pair(JpsiBrothers, trueLifePair);
      return result;
    } 
    else if (Jpsimom->numberOfMothers()<=0) {
      if (isAbHadron(Jpsimom->pdgId())) {  
	Jpsimom_final = Jpsimom;
	aBhadron = true;
      }
    } 

    else {
      reco::GenParticleRef Jpsigrandmom = findMotherRef(Jpsimom->motherRef(), Jpsimom->pdgId());   
      if (isAbHadron(Jpsimom->pdgId())) {
        if (Jpsigrandmom.isNonnull() && isAMixedbHadron(Jpsimom->pdgId(),Jpsigrandmom->pdgId())) {       
	  Jpsimom_final = Jpsigrandmom;
        } 
        else {                  
	Jpsimom_final = Jpsimom;
        }
        aBhadron = true;
      } 

      else if (Jpsigrandmom.isNonnull() && isAbHadron(Jpsigrandmom->pdgId()))  {  
        if (Jpsigrandmom->numberOfMothers()<=0) {
	  Jpsimom_final = Jpsigrandmom;
        } 
        else { 
          reco::GenParticleRef JpsiGrandgrandmom = findMotherRef(Jpsigrandmom->motherRef(), Jpsigrandmom->pdgId());
          if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(Jpsigrandmom->pdgId(),JpsiGrandgrandmom->pdgId())) {
	    Jpsimom_final = JpsiGrandgrandmom;
          } 
          else {
	    Jpsimom_final = Jpsigrandmom;
          }
        }
        aBhadron = true;
      }

      //This is to forcefully find the b-like mother of Jpsi 
      else if (Jpsigrandmom.isNonnull() && Jpsigrandmom->numberOfMothers()>0){
	reco::GenParticleRef JpsiGrandgrandmom = findMotherRef(Jpsigrandmom->motherRef(), Jpsigrandmom->pdgId());
	if(JpsiGrandgrandmom.isNonnull() && isAbHadron(JpsiGrandgrandmom->pdgId())){
          Jpsimom_final = JpsiGrandgrandmom;
	  aBhadron = true;
	}
      }

    }
    if (!aBhadron) {
      Jpsimom_final = Jpsimom;
    }

    if (Jpsimom_final.isNonnull()){
      trueVtxMom.SetXYZ(Jpsimom_final->vertex().x(),Jpsimom_final->vertex().y(),Jpsimom_final->vertex().z());
      if(_genealogyInfo && Reco_3mu_size>0
	 ){
	JpsiBrothers = GenBrothers(Jpsimom_final, genJpsi->pdgId() );
      }
      JpsiBrothers.insert(JpsiBrothers.begin(), Jpsimom_final);
    }
    
    TVector3 vdiff = trueVtx - trueVtxMom;
    trueLife = vdiff.Perp()*JpsiPDGMass/trueP.Perp();
    trueLife3D = vdiff.Mag()*JpsiPDGMass/trueP.Mag();

  }

  std::pair<float, float> trueLifePair = std::make_pair(trueLife, trueLife3D);
  std::pair<std::vector<reco::GenParticleRef>, std::pair<float, float> > result = std::make_pair(JpsiBrothers, trueLifePair);
  return result;

};