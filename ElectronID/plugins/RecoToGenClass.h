#ifndef jprtag
#define jprtag

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>


class JPsiReco{
	public:
	TFile *g;
	TTree *tree_JpsiUpsilon;
	inline void initJPsiRecoFillerObject();
	inline void WriteInFileAndCloseIt();
	float Mass;
	inline JPsiReco();

	inline ~JPsiReco();
};

#endif

JPsiReco::JPsiReco()
{
	g=new TFile("JpsiUpsilonFill.root","recreate");
	tree_JpsiUpsilon= new TTree("tree_JpsiUpsilon","JU ntuple");
	tree_JpsiUpsilon->Branch("Mass",&Mass,"Mass/F");
}

void JPsiReco::initJPsiRecoFillerObject(){
	Mass=-777.0;	
}

void JPsiReco::WriteInFileAndCloseIt(){ 
	g->cd();
        tree_JpsiUpsilon->Write();
        g->Close();
}

JPsiReco::~JPsiReco(){}	

#ifndef rtgtag
#define rtgtag

class RecoToGenFiller{
	public:
	inline RecoToGenFiller(const TString &);
	inline ~RecoToGenFiller();
	inline void initRecoToGenFillerObject();
	inline void WriteInFileAndCloseIt();
	TTree *tree_purity;
	TTree *tree_purity_GenLepB;
	TTree *tree_purity_GenPi;
	TTree *tree_purity_GenX;
	TTree *tree_purity_NonGenElec;
        TTree *tree_purity_NonGenPi;
        TTree *tree_purity_NonGenX;
	TFile *f;
	int origin;
        int nPV;
	int type;
        float ecalEnergy;
        float pIn;
        float pOut;
        float fbrem;
        float insidePFJet;
        float ptGen;
        float etaGen;
        float phiGen;
        float ptRel;
        float mindistJR;
        float pGen;
	float dRJetReco;
	int gen_match_;
	float RecoPt;
        float RecoEta;
        float RecoPhi;
        float RecoE;
	float TotalEnergy;
	int inRecoJet;
	int isMatchedWithASeed;
	int isMatchedWithAGedGsfElec;
	int Vtx;
        int pdgId;
	int tppdgId;
	float EtotOvePin;
	float EClusOverPout;
        float EBremOverDeltaP;
        float logSigmaEtaEta;
        float DeltaEtaTrackEcalSeed;
        float HOverE;
        float Chi2GSF;
        float CHi2KF;
        float nHits;
        float SigmaPtOverPt;
        float lnPt;
	float deta;
        float dphi;
        float detacalo;
        float see;
	float spp;
        float PreShowerOverRaw;
        float etawidth;
        float phiwidth;
        float e1x5e5x5;
        float R9;
        float HoE;
        float EoP;
        float IoEmIoP;
        float eleEoPout;
        float d0;
        float ip3d;
        float ip3dSig;
	float mva_e_pi;
	float weight;
	int fromConversion;
	float sip2d;
        float sip3d;
        float deltaR;
        float etaRel;
        float ratio;
        float ratioRel;
	int fl;
};

RecoToGenFiller::RecoToGenFiller(const TString & tag){
	f=new TFile("RecoToGenFill_"+tag+".root","recreate");
	tree_purity_GenLepB = new TTree("tree_purity_GenLepB","Reconst ntuple");
        tree_purity_GenLepB->Branch("origin",&origin,"origin/I");
        tree_purity_GenLepB->Branch("nPV",&nPV,"nPV/I");
	tree_purity_GenLepB->Branch("type",&type,"type/I");
        tree_purity_GenLepB->Branch("ecalEnergy",&ecalEnergy,"ecalEnergy/F");
        tree_purity_GenLepB->Branch("pIn",&pIn,"pIn/F");
        tree_purity_GenLepB->Branch("pOut",&pOut,"pOut/F");
        tree_purity_GenLepB->Branch("insidePFJet",&insidePFJet,"insidePFJet/F");
        tree_purity_GenLepB->Branch("ptGen",&ptGen,"ptGen/F");
        tree_purity_GenLepB->Branch("etaGen",&etaGen,"etaGen/F");
        tree_purity_GenLepB->Branch("phiGen",&phiGen,"phiGen/F");
        tree_purity_GenLepB->Branch("pGen",&pGen,"pGen/F");
        tree_purity_GenLepB->Branch("dRJetReco",&dRJetReco,"dRJetReco/F");
        tree_purity_GenLepB->Branch("ptRel",&ptRel,"ptRel/F");
        tree_purity_GenLepB->Branch("mindistJR",&mindistJR,"mindistJR/F");	
	tree_purity_GenLepB->Branch("gen_match_",&gen_match_,"gen_match_/I");
	tree_purity_GenLepB->Branch("inRecoJet",&inRecoJet,"inRecoJet/I");
	tree_purity_GenLepB->Branch("isMatchedWithASeed",&isMatchedWithASeed,"isMatchedWithASeed/I");
	tree_purity_GenLepB->Branch("pt",&RecoPt,"pt/F");
        tree_purity_GenLepB->Branch("eta",&RecoEta,"eta/F");
        tree_purity_GenLepB->Branch("phi",&RecoPhi,"phi/F");
        tree_purity_GenLepB->Branch("RecoE",&RecoE,"RecoE/F");
	tree_purity_GenLepB->Branch("TotalEnergy",&TotalEnergy,"TotalEnergy/F");
	tree_purity_GenLepB->Branch("isMatchedWithAGedGsfElec",&isMatchedWithAGedGsfElec,"isMatchedWithAGedGsfElec/I");
	tree_purity_GenLepB->Branch("Vtx",&Vtx,"Vtx/I");
        tree_purity_GenLepB->Branch("pdgId",&pdgId,"pdgId/I");
	tree_purity_GenLepB->Branch("tppdgId",&tppdgId,"tppdgId/I");

//ElecTraining
	tree_purity_GenLepB->Branch("EtotOvePin",&EtotOvePin,"EtotOvePin/F");
        tree_purity_GenLepB->Branch("EClusOverPout",&EClusOverPout,"EClusOverPout/F");
        tree_purity_GenLepB->Branch("fbrem",&fbrem,"fbrem/F");
        tree_purity_GenLepB->Branch("EBremOverDeltaP",&EBremOverDeltaP,"EBremOverDeltaP/F");
        tree_purity_GenLepB->Branch("logSigmaEtaEta",&logSigmaEtaEta,"logSigmaEtaEta/F");
        tree_purity_GenLepB->Branch("DeltaEtaTrackEcalSeed",&DeltaEtaTrackEcalSeed,"DeltaEtaTrackEcalSeed/F");
        tree_purity_GenLepB->Branch("gsfchi2",&Chi2GSF,"gsfchi2/F");
        tree_purity_GenLepB->Branch("kfchi2",&CHi2KF,"kfchi2/F");
        tree_purity_GenLepB->Branch("kfhits",&nHits,"kfhits/F");
        tree_purity_GenLepB->Branch("SigmaPtOverPt",&SigmaPtOverPt,"SigmaPtOverPt/F");
	tree_purity_GenLepB->Branch("deta",&deta,"deta/F");
        tree_purity_GenLepB->Branch("dphi",&dphi,"dphi/F");
        tree_purity_GenLepB->Branch("detacalo",&detacalo,"detacalo/F");
        tree_purity_GenLepB->Branch("see",&see,"see/F");
	tree_purity_GenLepB->Branch("spp",&spp,"spp/F");
        tree_purity_GenLepB->Branch("PreShowerOverRaw",&PreShowerOverRaw,"PreShowerOverRaw/F");
        tree_purity_GenLepB->Branch("etawidth",&etawidth,"etawidth/F");
        tree_purity_GenLepB->Branch("phiwidth",&phiwidth,"phiwidth/F");
        tree_purity_GenLepB->Branch("e1x5e5x5",&e1x5e5x5,"e1x5e5x5/F");
        tree_purity_GenLepB->Branch("R9",&R9,"R9/F");
        tree_purity_GenLepB->Branch("HoE",&HoE,"HoE/F");
        tree_purity_GenLepB->Branch("EoP",&EoP,"EoP/F");
        tree_purity_GenLepB->Branch("IoEmIoP",&IoEmIoP,"IoEmIoP/F");
        tree_purity_GenLepB->Branch("ip3d",&ip3d,"ip3d/F");
        tree_purity_GenLepB->Branch("ip3dSig",&ip3dSig,"ip3dSig/F");
	tree_purity_GenLepB->Branch("mva_e_pi",&mva_e_pi,"mva_e_pi/F");
	tree_purity_GenLepB->Branch("weight",&weight,"weight/F");
	tree_purity_GenLepB->Branch("fromConversion",&fromConversion,"fromConversion/I");
// BTagb training
        tree_purity_GenLepB->Branch("sip2d",&sip2d,"sip2d/F");
        tree_purity_GenLepB->Branch("sip3d",&sip3d,"sip3d/F");
        tree_purity_GenLepB->Branch("deltaR",&deltaR,"deltaR/F");
        tree_purity_GenLepB->Branch("etaRel",&etaRel,"etaRel/F");
        tree_purity_GenLepB->Branch("ratio",&ratio,"ratio/F");
        tree_purity_GenLepB->Branch("ratioRel",&ratioRel,"ratioRel/F");
	tree_purity_GenLepB->Branch("fl",&fl,"fl/I");
	tree_purity_GenPi=(TTree*)tree_purity_GenLepB->CloneTree();
	tree_purity_GenX=(TTree*)tree_purity_GenLepB->CloneTree();
	tree_purity_NonGenElec=(TTree*)tree_purity_GenLepB->CloneTree();
	tree_purity_NonGenPi=(TTree*)tree_purity_GenLepB->CloneTree();
	tree_purity_NonGenX=(TTree*)tree_purity_GenLepB->CloneTree();
        tree_purity=(TTree*)tree_purity_GenLepB->CloneTree();

	tree_purity_GenPi->SetName("tree_purity_GenPi");
        tree_purity_GenX->SetName("tree_purity_GenX");
        tree_purity_NonGenElec->SetName("tree_purity_NonGenElec");
        tree_purity_NonGenPi->SetName("tree_purity_NonGenPi");
        tree_purity_NonGenX->SetName("tree_purity_NonGenX");	
        tree_purity->SetName("tree_purity");
	
}

void RecoToGenFiller::initRecoToGenFillerObject(){
	tppdgId=-777;	
	origin=				-777;
	nPV                         =-1;
	type				=-1;
        ecalEnergy                  = -777.0;
        pIn                         = -777.0;
        pOut                        = -777.0;
        fbrem                       = -777.0;
        insidePFJet                 = 0;
        ptGen                       =-777.0;
        etaGen                      =-777.0;
        phiGen                      =-777.0;
        pGen                        =-777.0;

	gen_match_ = -1;
	RecoPt = -777.0;
        RecoEta = -777.0;
        RecoPhi = -777.0;
	RecoE = -777.0;
	TotalEnergy = -777.0;
	inRecoJet = 0;
	isMatchedWithASeed=0;
	isMatchedWithAGedGsfElec=0;	
	Vtx=-777;
        pdgId=-777;
	EtotOvePin=-777;
        EClusOverPout=-777;
        EBremOverDeltaP=-777;;
        logSigmaEtaEta=-777;;
        DeltaEtaTrackEcalSeed=-777;;
        HOverE=-777;;
        Chi2GSF=-777;;
        CHi2KF=-777;;
        nHits=-777;;
        SigmaPtOverPt=-777;;


	deta=-777;
        dphi=-777;
        detacalo=-777;
        see=-777;
	spp=-777.0;
   	PreShowerOverRaw=-777.0;
	etawidth=-777;
        phiwidth=-777;
        e1x5e5x5=-777;
        R9=-777;
        HoE=-777;
        EoP=-777;
        IoEmIoP=-777;
	ip3d=-777;
	ip3dSig=-777;
	mva_e_pi = -777.0;
	weight=777.0;
	fromConversion=-777.0;
        sip2d = -777.0;
        sip3d = -777.0;
        deltaR = -777.0;
        ptRel = -777.0;
        etaRel = -777.0;
        ratio = -777.0;
        ratioRel = -777.0;
	fl=-777.0;
}


void RecoToGenFiller::WriteInFileAndCloseIt(){
	f->cd();
	tree_purity->Write();
	tree_purity_GenLepB->Write();
	tree_purity_GenPi->Write();
        tree_purity_GenX->Write();
        tree_purity_NonGenElec->Write();
        tree_purity_NonGenPi->Write();
        tree_purity_NonGenX->Write();
	f->Close();
}
RecoToGenFiller::~RecoToGenFiller(){};

#endif
