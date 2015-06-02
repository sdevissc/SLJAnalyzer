#pragma once
#ifndef GENTORECOCLASS_H_
#define GENTORECOCLASS_H_


#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

class GenToRecoFiller{
	public:
	inline GenToRecoFiller(const TString &);
	inline ~GenToRecoFiller();
	inline void initGenToRecoFillerObject();
	inline void WriteInFileAndCloseIt();
	TTree *tree_efficiency;
	TFile *f;
	int pdgId;
	int origin;
        int nPV;
        float ecalEnergy;
        float pIn;
        float pOut;
        float EemPinRatio;
        float EemPoutRatio;
        float GGE_p;
        float GGE_pt;
        float GGE_eta;
        float GGE_phi;
        float GGE_energy;
        float fbrem;
        float dRGsfTrackElectron;
        float inversedRFirstLastHit;
        float radiusFirstHit;
        float zFirstHit;
        float insidePFJet;
        float ptGen;
        float etaGen;
        float phiGen;
        float ptRel;
        float mindistJR;
        float pGen;
	float dRJetReco;
	int gen_match_;
	float gsftrkSeedPt;
	float gsftrkSeedEta;
	float gsftrkSeedPhi;
	float gedgsfElecPt;
        float gedgsfElecEta;
        float gedgsfElecPhi;
	int inRecoJet;
	int isMatchedWithASeed;
	int isEcalDrivenSeeded;
	int isDeltaRMatchedWithASeed;
        int isDeltaREcalDrivenSeeded;
	int isMatchedWithAGedGsfElec;
	int isMatchedWithAPFElec;
	float mva_e_pi;	
	float sharedHits;
        float mva_e_pi_PF;
	float PFElecPt;
        float PFElecEta;
        float PFElecPhi;
};

GenToRecoFiller::GenToRecoFiller(const TString & tag){
	f=new TFile("GenToRecoFill.root","recreate");
	tree_efficiency = new TTree("tree_efficiency","Reconst ntuple");
	tree_efficiency->Branch("pdgId",&pdgId,"pdgId/I");	
        tree_efficiency->Branch("origin",&origin,"origin/I");
        tree_efficiency->Branch("nPV",&nPV,"nPV/I");
        tree_efficiency->Branch("ecalEnergy",&ecalEnergy,"ecalEnergy/F");
        tree_efficiency->Branch("pIn",&pIn,"pIn/F");
        tree_efficiency->Branch("pOut",&pOut,"pOut/F");
        tree_efficiency->Branch("EemPinRatio",&EemPinRatio,"EemPinRatio/F");
        tree_efficiency->Branch("EemPoutRatio",&EemPoutRatio,"EemPoutRatio/F");
        tree_efficiency->Branch("GGE_p",&GGE_p,"GGE_p/F");
        tree_efficiency->Branch("GGE_pt",&GGE_pt,"GGE_pt/F");
        tree_efficiency->Branch("GGE_eta",&GGE_eta,"GGE_eta/F");
        tree_efficiency->Branch("GGE_phi",&GGE_phi,"GGE_phi/F");
        tree_efficiency->Branch("GGE_energy",&GGE_energy,"GGE_energy/F");
        tree_efficiency->Branch("fbrem",&fbrem,"fbrem/F");
        tree_efficiency->Branch("dRGsfTrackElectron",&dRGsfTrackElectron,"dRGsfTrackElectron/F");
        tree_efficiency->Branch("inversedRFirstLastHit",&inversedRFirstLastHit,"inversedRFirstLastHit/F");
        tree_efficiency->Branch("radiusFirstHit",&radiusFirstHit,"radiusFirstHit/F");
        tree_efficiency->Branch("zFirstHit",&zFirstHit,"zFirstHit/F");
        tree_efficiency->Branch("insidePFJet",&insidePFJet,"insidePFJet/F");
        tree_efficiency->Branch("ptGen",&ptGen,"ptGen/F");
        tree_efficiency->Branch("etaGen",&etaGen,"etaGen/F");
        tree_efficiency->Branch("phiGen",&phiGen,"phiGen/F");
        tree_efficiency->Branch("pGen",&pGen,"pGen/F");
        tree_efficiency->Branch("dRJetReco",&dRJetReco,"dRJetReco/F");
        tree_efficiency->Branch("ptRel",&ptRel,"ptRel/F");
        tree_efficiency->Branch("mindistJR",&mindistJR,"mindistJR/F");	
	tree_efficiency->Branch("gen_match_",&gen_match_,"gen_match_/I");
	tree_efficiency->Branch("gsftrkSeedPt",&gsftrkSeedPt,"gsftrkSeedPt/F");
        tree_efficiency->Branch("gsftrkSeedEta",&gsftrkSeedEta,"gsftrkSeedEta/F");
        tree_efficiency->Branch("gsftrkSeedPhi",&gsftrkSeedPhi,"gsftrkSeedPhi/F");
	tree_efficiency->Branch("inRecoJet",&inRecoJet,"inRecoJet/I");
	tree_efficiency->Branch("isMatchedWithASeed",&isMatchedWithASeed,"isMatchedWithASeed/I");
	tree_efficiency->Branch("isEcalDrivenSeeded",&isEcalDrivenSeeded,"isEcalDrivenSeeded/I");
        tree_efficiency->Branch("isDeltaRMatchedWithASeed",&isDeltaRMatchedWithASeed,"isDeltaRMatchedWithASeed/I");
        tree_efficiency->Branch("isDeltaREcalDrivenSeeded",&isDeltaREcalDrivenSeeded,"isDeltaREcalDrivenSeeded/I");
	tree_efficiency->Branch("gedgsfElecPt",&gedgsfElecPt,"gedgsfElecPt/F");
        tree_efficiency->Branch("gedgsfElecEta",&gedgsfElecEta,"gedgsfElecEta/F");
        tree_efficiency->Branch("gedgsfElecPhi",&gedgsfElecPhi,"gedgsfElecPhi/F");
	tree_efficiency->Branch("isMatchedWithAGedGsfElec",&isMatchedWithAGedGsfElec,"isMatchedWithAGedGsfElec/I");
	tree_efficiency->Branch("mve_e_pi",&mva_e_pi,"mva_e_pi/F");		
	tree_efficiency->Branch("sharedHits",&sharedHits,"sharedHits/F");
	tree_efficiency->Branch("isMatchedWithAPFElec",&isMatchedWithAPFElec,"isMatchedWithAPFElec/I");
        tree_efficiency->Branch("mva_e_pi_PF",&mva_e_pi_PF,"mva_e_pi_PF/F");
        tree_efficiency->Branch("PFElecPt",&PFElecPt,"PFElecPt/F");
        tree_efficiency->Branch("PFElecEta",&PFElecEta,"PFElecEta/F");
        tree_efficiency->Branch("PFElecPhi",&PFElecPhi,"PFElecPhi/F");
}

void GenToRecoFiller::initGenToRecoFillerObject(){
	pdgId= -777;	
	nPV                         =-1;
        ecalEnergy                  = -777.0;
        pIn                         = -777.0;
        pOut                        = -777.0;
        EemPinRatio                 = -777.0;
        EemPoutRatio                = -777.0;
        GGE_p                      = -777.0;
        GGE_pt                      = -777.0;
        GGE_eta                     = -777.0;
        GGE_phi                     = -777.0;
        GGE_energy                  = -777.0;
        fbrem                       = -777.0;
        dRGsfTrackElectron          = -777.0;
        insidePFJet                 = 0;
        ptGen                       =-777.0;
        etaGen                      =-777.0;
        phiGen                      =-777.0;
        pGen                        =-777.0;

	gen_match_ = -1;
	gsftrkSeedPt = -777.0;
	gsftrkSeedEta = -777.0;
	gsftrkSeedPhi = -777.0;
	gedgsfElecPt = -777.0;
        gedgsfElecEta = -777.0;
        gedgsfElecPhi = -777.0;
	inRecoJet = 0;
	isMatchedWithASeed=0;
	isEcalDrivenSeeded=0;
        isDeltaRMatchedWithASeed=0;
        isDeltaREcalDrivenSeeded=0;
	isMatchedWithAGedGsfElec=0;	
	origin=-777;
	mva_e_pi=-777;
	sharedHits=-777.0;
	isMatchedWithAPFElec=-0;
        mva_e_pi_PF=-777.0;
        PFElecPt=-777.0;
        PFElecEta=-777.0;
        PFElecPhi=-777.0;
	
}


void GenToRecoFiller::WriteInFileAndCloseIt(){
	f->cd();
	tree_efficiency->Write();
	f->Close();
}
GenToRecoFiller::~GenToRecoFiller(){};

#endif
