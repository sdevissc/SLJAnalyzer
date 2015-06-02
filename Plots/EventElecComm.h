//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan 27 12:14:28 2013 by ROOT version 5.32/00
// from TTree tree_efficiency/Reconst ntuple
// found on file: ElectronTree.root
//////////////////////////////////////////////////////////

#ifndef EventElecComm_h
#define EventElecComm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class EventElecComm {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree *tree;
   // Declaration of leaf types
Int_t           origin;
Int_t		pdgId;	
   Int_t           nPV;
   Float_t         ecalEnergy;
   Float_t         pIn;
   Float_t         pOut;
   Float_t         EemPinRatio;
   Float_t         EemPoutRatio;
   Float_t         GGE_p;
   Float_t         GGE_pt;
   Float_t         GGE_eta;
   Float_t         GGE_phi;
   Float_t         GGE_energy;
   Float_t         fbrem;
   Float_t         dRGsfTrackElectron;
   Float_t         inversedRFirstLastHit;
   Float_t         radiusFirstHit;
   Float_t         zFirstHit;
   Float_t         insidePFJet;
   Float_t         ptGen;
   Float_t         etaGen;
   Float_t         phiGen;
   Float_t         pGen;
   Float_t         dRJetReco;
   Float_t         ptRel;
   Float_t         mindistJR;
   Int_t           gen_match_;
   Float_t         gsftrkSeedPt;
   Float_t         gsftrkSeedEta;
   Float_t         gsftrkSeedPhi;
   Int_t           inRecoJet;
   Int_t           isMatchedWithASeed;
   Int_t 	   isEcalDrivenSeeded;
   Float_t         gedgsfElecPt;
   Float_t         gedgsfElecEta;
   Float_t         gedgsfElecPhi;
   Int_t           isMatchedWithAGedGsfElec;	
   Int_t           isMatchedWithAPFElec;
   Float_t 	   mva_e_pi;
   Float_t         mva_e_pi_PF;
   // List of branches
TBranch        *b_origin;   //!
   TBranch        *b_pdgId;   //!	
   TBranch        *b_nPV;   //!
   TBranch        *b_ecalEnergy;   //!
   TBranch        *b_pIn;   //!
   TBranch        *b_pOut;   //!
   TBranch        *b_EemPinRatio;   //!
   TBranch        *b_EemPoutRatio;   //!
   TBranch        *b_GGE_p;   //!
   TBranch        *b_GGE_pt;   //!
   TBranch        *b_GGE_eta;   //!
   TBranch        *b_GGE_phi;   //!
   TBranch        *b_GGE_energy;   //!
   TBranch        *b_fbrem;   //!
   TBranch        *b_dRGsfTrackElectron;   //!
   TBranch        *b_inversedRFirstLastHit;   //!
   TBranch        *b_radiusFirstHit;   //!
   TBranch        *b_zFirstHit;   //!
   TBranch        *b_insidePFJet;   //!
   TBranch        *b_ptGen;   //!
   TBranch        *b_etaGen;   //!
   TBranch        *b_phiGen;   //!
   TBranch        *b_pGen;   //!
   TBranch        *b_dRJetReco;   //!
   TBranch        *b_ptRel;   //!
   TBranch        *b_mindistJR;   //!
   TBranch        *b_gen_match_;   //!
   TBranch        *b_gsftrkSeedPt;   //!
   TBranch        *b_gsftrkSeedEta;   //!
   TBranch        *b_gsftrkSeedPhi;   //!
   TBranch        *b_inRecoJet;   //!
   TBranch        *b_isMatchedWithASeed;   //!
   TBranch 	  *b_isEcalDrivenSeeded;
   TBranch        *b_gedgsfElecPt;   //!
   TBranch        *b_gedgsfElecEta;   //!
   TBranch        *b_gedgsfElecPhi;   //!
   TBranch        *b_isMatchedWithAGedGsfElec;   //!	
   TBranch        *b_isMatchedWithAPFElec;   //! 
   TBranch        *b_mva_e_pi;   //!
   TBranch        *b_mva_e_pi_PF;

   EventElecComm( TString filename );
   virtual ~EventElecComm();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
//   virtual Float_t StringToFloat(TString);	
};

#endif

#ifdef EventElecComm_cxx

EventElecComm::EventElecComm( TString filename )
{

  tree = 0;

  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f) {
         f = new TFile(filename);
      }
      tree = (TTree*)gDirectory->Get("tree_efficiency");

   }
   Init(tree);
}


/*EventElecComm::EventElecComm(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ElectronTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ElectronTree.root");
      }
      f->GetObject("tree_efficiency",tree);

   }
   Init(tree);
}*/

EventElecComm::~EventElecComm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventElecComm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventElecComm::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventElecComm::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("origin", &origin, &b_origin);
   fChain->SetBranchAddress("pdgId", &pdgId, &b_pdgId);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("ecalEnergy", &ecalEnergy, &b_ecalEnergy);
   fChain->SetBranchAddress("pIn", &pIn, &b_pIn);
   fChain->SetBranchAddress("pOut", &pOut, &b_pOut);
   fChain->SetBranchAddress("EemPinRatio", &EemPinRatio, &b_EemPinRatio);
   fChain->SetBranchAddress("EemPoutRatio", &EemPoutRatio, &b_EemPoutRatio);
   fChain->SetBranchAddress("GGE_p", &GGE_p, &b_GGE_p);
   fChain->SetBranchAddress("GGE_pt", &GGE_pt, &b_GGE_pt);
   fChain->SetBranchAddress("GGE_eta", &GGE_eta, &b_GGE_eta);
   fChain->SetBranchAddress("GGE_phi", &GGE_phi, &b_GGE_phi);
   fChain->SetBranchAddress("GGE_energy", &GGE_energy, &b_GGE_energy);
   fChain->SetBranchAddress("fbrem", &fbrem, &b_fbrem);
   fChain->SetBranchAddress("dRGsfTrackElectron", &dRGsfTrackElectron, &b_dRGsfTrackElectron);
   fChain->SetBranchAddress("inversedRFirstLastHit", &inversedRFirstLastHit, &b_inversedRFirstLastHit);
   fChain->SetBranchAddress("radiusFirstHit", &radiusFirstHit, &b_radiusFirstHit);
   fChain->SetBranchAddress("zFirstHit", &zFirstHit, &b_zFirstHit);
   fChain->SetBranchAddress("insidePFJet", &insidePFJet, &b_insidePFJet);
   fChain->SetBranchAddress("ptGen", &ptGen, &b_ptGen);
   fChain->SetBranchAddress("etaGen", &etaGen, &b_etaGen);
   fChain->SetBranchAddress("phiGen", &phiGen, &b_phiGen);
   fChain->SetBranchAddress("pGen", &pGen, &b_pGen);
   fChain->SetBranchAddress("dRJetReco", &dRJetReco, &b_dRJetReco);
   fChain->SetBranchAddress("ptRel", &ptRel, &b_ptRel);
   fChain->SetBranchAddress("mindistJR", &mindistJR, &b_mindistJR);
   fChain->SetBranchAddress("gen_match_", &gen_match_, &b_gen_match_);
   fChain->SetBranchAddress("gsftrkSeedPt", &gsftrkSeedPt, &b_gsftrkSeedPt);
   fChain->SetBranchAddress("gsftrkSeedEta", &gsftrkSeedEta, &b_gsftrkSeedEta);
   fChain->SetBranchAddress("gsftrkSeedPhi", &gsftrkSeedPhi, &b_gsftrkSeedPhi);
   fChain->SetBranchAddress("inRecoJet", &inRecoJet, &b_inRecoJet);
   fChain->SetBranchAddress("isMatchedWithASeed", &isMatchedWithASeed, &b_isMatchedWithASeed);
   fChain->SetBranchAddress("isEcalDrivenSeeded",&isEcalDrivenSeeded,&b_isEcalDrivenSeeded);
   fChain->SetBranchAddress("gedgsfElecPt", &gedgsfElecPt, &b_gedgsfElecPt);
   fChain->SetBranchAddress("gedgsfElecEta", &gedgsfElecEta, &b_gedgsfElecEta);
   fChain->SetBranchAddress("gedgsfElecPhi", &gedgsfElecPhi, &b_gedgsfElecPhi);
   fChain->SetBranchAddress("isMatchedWithAGedGsfElec", &isMatchedWithAGedGsfElec, &b_isMatchedWithAGedGsfElec);
   fChain->SetBranchAddress("isMatchedWithAPFElec", &isMatchedWithAPFElec, &b_isMatchedWithAPFElec);
   fChain->SetBranchAddress("mve_e_pi", &mva_e_pi, &b_mva_e_pi);
   fChain->SetBranchAddress("mva_e_pi_PF", &mva_e_pi_PF, &b_mva_e_pi_PF);
   Notify();
}

Bool_t EventElecComm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventElecComm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventElecComm::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif 

/*
Float_t EventElecComm::StringToFloat(TString VarString)
{

  Float_t VarFloat = 0;

  if( VarString == "origin")VarFloat =origin;
  else if( VarString == "ecalEnergy")VarFloat =ecalEnergy;
  else if( VarString == "pIn")VarFloat =pIn;
  else if( VarString == "pOut")	VarFloat =pOut;
  else if( VarString == "GSF_GSFSuperClusterEnergy")VarFloat =GSF_GSFSuperClusterEnergy;
  else if( VarString == "GSF_PFSuperClusterEnergy")VarFloat =GSF_PFSuperClusterEnergy;
  else if( VarString == "GSF_momentum")VarFloat =GSF_momentum;
  else if( VarString == "PF_pt")VarFloat =PF_pt;
  else if( VarString == "PF_eta")VarFloat =PF_eta;
  else if( VarString == "PF_energy")VarFloat =PF_energy;
  else if( VarString == "GSF_pt")VarFloat =GSF_pt;
  else if( VarString == "GSF_eta")VarFloat =GSF_eta;
  else if( VarString == "GSF_phi")VarFloat =GSF_phi;
  else if( VarString == "GSF_Track_pt")VarFloat =GSF_Track_pt;
  else if( VarString == "GSF_Track_eta")VarFloat =GSF_Track_eta;
  else if( VarString == "GSF_Track_phi")VarFloat =GSF_Track_phi;	
  else if( VarString == "GSF_energy")VarFloat =GSF_energy;
  else if( VarString == "fbrem")VarFloat =fbrem;
  else if( VarString == "dRGsfTrackElectron")VarFloat =dRGsfTrackElectron;
  else if( VarString == "mva")VarFloat =mva;
  else if( VarString == "newmva")VarFloat =newmva;
  else if( VarString == "inversedRFirstLastHit")VarFloat =inversedRFirstLastHit;
  else if( VarString == "radiusFirstHit")VarFloat =radiusFirstHit;
  else if( VarString == "zFirstHit")VarFloat =zFirstHit;
  else if( VarString == "insidePFJet")VarFloat =insidePFJet;
  else if( VarString == "ptGen")VarFloat =ptGen;
  else if( VarString == "etaGen")VarFloat =etaGen;
  else if( VarString == "phiGen")VarFloat =phiGen;	
  else if( VarString == "pGen")VarFloat =pGen;
  else if( VarString == "dRJetReco")VarFloat =dRJetReco;
  else if( VarString == "ptRel")VarFloat =ptRel;
  else if( VarString == "pid")VarFloat =pid;
  else if( VarString == "mindistJR")VarFloat =mindistJR;
  else if( VarString == "GSF_fbrem")VarFloat =GSF_fbrem;
  else if( VarString == "GSF_dRGsfTrackElectron")VarFloat =GSF_dRGsfTrackElectron;
  else if( VarString == "PF_EemPinRatio")VarFloat =PF_EemPinRatio;
  else if( VarString == "PF_fbrem")VarFloat =PF_fbrem;
  else if( VarString == "PF_dRGsfTrackElectron")VarFloat =PF_dRGsfTrackElectron;
  else if( VarString == "PF_hitCondition")VarFloat =PF_hitCondition;
  else if( VarString == "GSF_hitCondition")VarFloat =GSF_hitCondition;
  else if( VarString == "GSF_mva")VarFloat =GSF_mva;
// DeltaR GSF vs True
  else if( VarString == "GSF_GSFElec_Gen_deltaR")VarFloat =sqrt(pow(GSF_eta-etaGen,2)+pow(GSF_phi-phiGen,2));
  else if( VarString == "GSF_GSFTrack_Gen_deltaR")VarFloat =sqrt(pow(GSF_Track_eta-etaGen,2)+pow(GSF_Track_phi-phiGen,2));

// recoed momentum vs pGen
  else if( VarString == "TrackPOverpGen")VarFloat =GSF_Track_pMode/pGen;
  else if( VarString == "TrackPtOverptGen")VarFloat =GSF_Track_pt/ptGen;

// Energy in the SC vs True E
  else if( VarString == "PFEscOverPgen")VarFloat =GSF_PFSuperClusterEnergy/pGen;
  else if( VarString == "GSFEscOverPgen")VarFloat =GSF_GSFSuperClusterEnergy/pGen;

// Energy in the SC vs Track P/Electron P
  else if( VarString == "PFEscOverP")VarFloat =GSF_PFSuperClusterEnergy/GSF_Track_pMode;
  else if( VarString == "GSFEscOverP")VarFloat =GSF_GSFSuperClusterEnergy/GSF_Track_pMode;
  else if( VarString == "EemPinRatio")VarFloat =(GSF_PFSuperClusterEnergy/GSF_Track_pMode-1)/(GSF_PFSuperClusterEnergy/GSF_Track_pMode+1);
return VarFloat;

}

*/
