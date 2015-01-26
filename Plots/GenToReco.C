#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1F.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <./EventElecComm.C>
#include <TPaveStats.h>
using namespace std;

const int Nbin2D=300; 

class HistoMaker
{
  public:
	HistoMaker(const TString &,const TString &);
	~HistoMaker();
	void SetStyle(TH1F *,int, int, const TString&,const TString&, float, float);
        void SetStyle(TProfile *,int, int, const TString&,const TString&, float, float);
	void Combine();

	TH1F *Pt_gen_all;
	TH1F *Pt_gen_MatchedWithASeed;
	TH1F *Pt_gen_MatchedWithAGedGsfElec;
        TH1F *Pt_gen_MatchedWithAPFElec;
        
        TH1F *Eta_gen_all;
	TH1F *Eta_gen_MatchedWithASeed;
	TH1F *Eta_gen_MatchedWithAGedGsfElec;
        TH1F *Eta_gen_MatchedWithAPFElec;

	TH1F *EfficiencyVsPt_Seed;
        TH1F *EfficiencyVsEta_Seed;
        TH1F *EfficiencyVsPt_GedGsfElec;
        TH1F *EfficiencyVsEta_GedGsfElec;
        TH1F *EfficiencyVsPt_PFElec;
        TH1F *EfficiencyVsEta_PFElec;
	
	TH1F *mva_e_pi,*mva_e_pi_PF;	

	TH1F* ForIntegral;
	TGraph *eEffvspiEff;
	float vecXmva[Nbin2D],vecYmva[Nbin2D];

};

HistoMaker::HistoMaker(const TString& inputname,const TString & tag)
{
	const int nbinpt=6,nbineta=13;
	 TH1::SetDefaultSumw2(kTRUE);
	cout<<"define the binning"<<endl;
//	float ptx[17]= {0.0 , 2.0 , 3.0 , 4.0 , 6.0 , 8.0 , 10.0 , 13.0 , 16.0  , 20.0 , 24.0 , 30.0 , 40.0 , 55.0 , 70.0 , 90.0 , 120.0};
	float etax[14]={-2.5 , -2.1  , -1.7 , -1.3 , -0.9 , -0.5 , -0.1 , 0.1 , 0.5 , 0.9 , 1.3 , 1.7 , 2.1 , 2.5 };
        float ptx[nbinpt+1]= { 2.0 , 5.0 , 10.0 , 20.0 , 30.0 , 50.0 , 120.0};

	 cout<<"define pt and eta basic histos"<<endl;
	ForIntegral=new TH1F(inputname+tag+"forIntegral","",Nbin2D,-100000,100000);	

        Pt_gen_all=new TH1F(inputname+tag+"Pt_gen_all","",nbinpt,ptx);
	cout<<"ok"<<endl;
        Eta_gen_all=new TH1F(inputname+tag+"Eta_gen_all","",nbineta,etax);
	
	Pt_gen_MatchedWithASeed=new TH1F(inputname+tag+"Pt_gen_MatchedWithSeed","",nbinpt,ptx);
        Eta_gen_MatchedWithASeed=new TH1F(inputname+tag+"Eta_gen_MatchedWithSeed","",nbineta,etax);
	
	Pt_gen_MatchedWithAGedGsfElec=new TH1F(inputname+tag+"Pt_gen_MatchedWithGedGsfElec","",nbinpt,ptx);
        Eta_gen_MatchedWithAGedGsfElec=new TH1F(inputname+tag+"Eta_gen_MatchedWithGedGsfElec","",nbineta,etax);
	
        Pt_gen_MatchedWithAPFElec=new TH1F(inputname+tag+"Pt_gen_MatchedWithPFElec","",nbinpt,ptx);
        Eta_gen_MatchedWithAPFElec=new TH1F(inputname+tag+"Eta_gen_MatchedWithPFElec","",nbineta,etax);

	 cout<<"define the efficiency histos"<<endl;
	EfficiencyVsPt_Seed=new TH1F(inputname+tag+"effptseed","",nbinpt,ptx);
        EfficiencyVsEta_Seed=new TH1F(inputname+tag+"effetaseed","",nbineta,etax);
        EfficiencyVsPt_GedGsfElec=new TH1F(inputname+tag+"effptged","",nbinpt,ptx);
        EfficiencyVsEta_GedGsfElec=new TH1F(inputname+tag+"effgedeta","",nbineta,etax);
	EfficiencyVsPt_PFElec=new TH1F(inputname+tag+"effptpf","",nbinpt,ptx);
        EfficiencyVsEta_PFElec=new TH1F(inputname+tag+"effetapf","",nbineta,etax);

	mva_e_pi=new TH1F(inputname+tag+"mva","",Nbin2D,-1,1.0);
        mva_e_pi_PF=new TH1F(inputname+tag+"mvapf","",Nbin2D,-1,1.0);

	SetStyle(Pt_gen_all,1,22,"p_{T} of the generated electron","",0.001,1.2);
	SetStyle(Eta_gen_all,1,22,"eta_{T} of the generated electron","",0.001,1.2);	

        SetStyle(Pt_gen_MatchedWithASeed,1,22,"p_{T} of the generated electron","",0.001,1.2);
        SetStyle(Pt_gen_MatchedWithAGedGsfElec,2,22,"p_{T} of the generated electron","",0.001,1.2);
        SetStyle(Pt_gen_MatchedWithAPFElec,2,22,"p_{T} of the generated electron","",0.001,1.2);

        SetStyle(Eta_gen_MatchedWithASeed,1,22,"#eta of the generated electron","",0.001,1.2);
        SetStyle(Eta_gen_MatchedWithAGedGsfElec,2,22,"#eta of the generated electron","",0.001,1.2);
        SetStyle(Eta_gen_MatchedWithAPFElec,2,22,"#eta of the generated electron","",0.001,1.2);

	SetStyle(EfficiencyVsPt_Seed,1,22,"p_{T} of the generated electron","",0.001,1.2);
        SetStyle(EfficiencyVsPt_GedGsfElec,2,22,"p_{T} of the generated electron","",0.001,1.2);
        SetStyle(EfficiencyVsPt_PFElec,2,22,"p_{T} of the generated electron","",0.001,1.2);
	
	SetStyle(EfficiencyVsEta_Seed,1,22,"#eta of the generated electron","",0.001,1.2);
	SetStyle(EfficiencyVsEta_GedGsfElec,2,22,"#eta of the generated electron","",0.001,1.2);
	SetStyle(EfficiencyVsEta_PFElec,2,22,"#eta of the generated electron","",0.001,1.2);
	

}


void HistoMaker::SetStyle(TH1F *histo,int color, int style, const TString& xtitle,const TString& ytitle, float min, float max){
	histo->SetLineColor(color);
        histo->SetMarkerColor(color);
        histo->SetMarkerStyle(style);
        histo->GetXaxis()->SetTitle(xtitle);
	histo->GetYaxis()->SetTitle(ytitle);
		histo->SetMinimum(min);
      	  	histo->SetMaximum(max);
}

void HistoMaker::SetStyle(TProfile *histo,int color, int style, const TString& xtitle,const TString& ytitle, float min, float max){
        histo->SetLineColor(color);
        histo->SetMarkerColor(color);
        histo->SetMarkerStyle(style);
        histo->GetXaxis()->SetTitle(xtitle);
        histo->GetYaxis()->SetTitle(ytitle);
        if(min!=-777.0 && max!=-777.0 ){
                histo->SetMinimum(min);
                histo->SetMaximum(max);
        }
}

HistoMaker::~HistoMaker()                 // destructor, just an example
{
}





void HistoMaker::Combine(){
	TH1::SetDefaultSumw2(kTRUE);
	EfficiencyVsPt_Seed->Add(Pt_gen_MatchedWithASeed);
        EfficiencyVsPt_Seed->Divide(Pt_gen_all);
        EfficiencyVsEta_Seed->Add(Eta_gen_MatchedWithASeed);
        EfficiencyVsEta_Seed->Divide(Eta_gen_all);

	EfficiencyVsPt_GedGsfElec->Add(Pt_gen_MatchedWithAGedGsfElec);
        EfficiencyVsPt_GedGsfElec->Divide(Pt_gen_all);
        EfficiencyVsEta_GedGsfElec->Add(Eta_gen_MatchedWithAGedGsfElec);
        EfficiencyVsEta_GedGsfElec->Divide(Eta_gen_all);

        EfficiencyVsPt_PFElec->Add(Pt_gen_MatchedWithAPFElec);
        EfficiencyVsPt_PFElec->Divide(Pt_gen_all);
        EfficiencyVsEta_PFElec->Add(Eta_gen_MatchedWithAPFElec);
        EfficiencyVsEta_PFElec->Divide(Eta_gen_all);

	
	cout<<EfficiencyVsPt_Seed->Integral()<<" "<<EfficiencyVsEta_Seed->Integral()<<" "<<EfficiencyVsPt_GedGsfElec->Integral()<<" "<<EfficiencyVsEta_GedGsfElec->Integral()<<endl;
}


class SLPlotter                   // begin declaration of the class
{
  public:                    // begin public section
    	SLPlotter(const TString &,const TString &);     // constructor
    	~SLPlotter();                  // destructor
	void initialize();
	void setTDRStyle();
	EventElecComm *Evt;
	int nentries;
	TFile *file;
        TTree *tp,*te;
	void getPlot(int,TString&); 
	void GetEfficiencies(int,TString &);
//	void GetEscOverPGen(TString &);
	void FillHistos();
	void TheFill(HistoMaker *,EventElecComm *);
	void MakeEffvsEffPlot(TString &);
	int nmvasteps;
//	TH1F* GetEffPvsPt();
//	TH1F* GetEffPvsEta();

	TString DefElecFromB;
	TString DefElecFromV;
	TString DefPiKaInJet;
	TString LimitBarrelEndcap;
	TString ProcessFile;
        TString ProcName;
	TString TheTag;
	HistoMaker *electron;
	HistoMaker *pion;
	HistoMaker *kaon;
	HistoMaker *electron_lowPU;
        HistoMaker *electron_midPU;
        HistoMaker *electron_highPU;
        HistoMaker *pion_lowPU;
        HistoMaker *pion_midPU;
        HistoMaker *pion_highPU;
	HistoMaker *kaon_lowPU;
        HistoMaker *kaon_midPU;
        HistoMaker *kaon_highPU;	
	HistoMaker *background;	
	

		
	TString ProcessName;
	
}
;

SLPlotter::SLPlotter(const TString &input,const TString &tag)
{
	ProcessFile = input;
	TheTag=tag;
	electron=new HistoMaker("elec",TheTag);
	pion=new HistoMaker("pion",TheTag);
	kaon=new HistoMaker("kaon",TheTag);
	electron_lowPU=new HistoMaker("elec_lpu",TheTag);
	electron_midPU=new HistoMaker("elec_mpu",TheTag);
	electron_highPU=new HistoMaker("elec_hpu",TheTag);
        pion_lowPU=new HistoMaker("pion_lpu",TheTag);
	pion_midPU=new HistoMaker("pion_mpu",TheTag);
	pion_highPU=new HistoMaker("pion_hpu",TheTag);
	kaon_lowPU=new HistoMaker("kaon_lpu",TheTag);
        kaon_midPU=new HistoMaker("kaon_mpu",TheTag);
        kaon_highPU=new HistoMaker("kaon_hpu",TheTag);
	background=new HistoMaker("bckg",TheTag);
	
}

void SLPlotter::initialize(){
        ifstream Processes;
        Processes.open(ProcessFile);
        string processline;
        Processes >>   ProcName;
        Evt = new EventElecComm(ProcName);
        nentries =/* 1500000;*/Evt->fChain->GetEntriesFast();
}


SLPlotter::~SLPlotter()                 // destructor, just an example
{
}

void SLPlotter::getPlot(int log,TString& process){
	initialize();
	ProcessName = process;
//	GetEscOverPGen(ProcessName);
	FillHistos();
	GetEfficiencies(log,ProcessName);	
	MakeEffvsEffPlot(ProcessName);
}


void SLPlotter::GetEfficiencies(int log,TString &process){
	TH1::SetDefaultSumw2(kTRUE);
	ProcessName = process;
	setTDRStyle();
        TH1::SetDefaultSumw2(kTRUE);
        TCanvas *can = new TCanvas("h", "bla",600,600);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.05);
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);
	TLegend *legend=new TLegend(0.65,0.8,0.95,0.96);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.020);
	

//-----------------------------------

	gPad->SetGridx();
        gPad->SetGridy();
	if(log==1)gPad->SetLogy(1);
	//Plotting Pt plot


	// GSFTrack_all vs GSFTracl_ecalSeeded
	legend->AddEntry(electron->EfficiencyVsPt_Seed,"e_{B}#rightarrow Seed (all PU)");
//	legend->AddEntry(electron->EfficiencyVsPt_GedGsfElec,"e_{B}#rightarrow GED Elec");
  //      legend->AddEntry(electron->EfficiencyVsPt_PFElec,"e_{B}#rightarrow PF Elec");	
        legend->AddEntry(electron_lowPU->EfficiencyVsPt_Seed,"e_{B}#rightarrow Seed (<30 PU)","p");
//        legend->AddEntry(electron_lowPU->EfficiencyVsPt_GedGsfElec,"GED Gsf elec low PU","p");
//        legend->AddEntry(electron_lowPU->EfficiencyVsPt_PFElec,"PF elec low PU","p");
        legend->AddEntry(electron_highPU->EfficiencyVsPt_Seed,"e_{B}#rightarrow Seed (>30 PU)","p");
//        legend->AddEntry(electron_highPU->EfficiencyVsPt_GedGsfElec,"GED Gsf elec high PU","p");
//        legend->AddEntry(electron_highPU->EfficiencyVsPt_PFElec,"PF elec high PU","p");

	electron->EfficiencyVsPt_Seed->SetMaximum(1.2);
	electron->EfficiencyVsPt_Seed->Draw("pe2");
//        electron->EfficiencyVsPt_GedGsfElec->Draw("pe1same");
//        electron->EfficiencyVsPt_PFElec->Draw("pe1same");
	
        electron->EfficiencyVsPt_Seed->SetMarkerStyle(22);
	electron->EfficiencyVsPt_Seed->SetFillColor(1);
	electron->EfficiencyVsPt_Seed->SetFillStyle(3004);
	electron->EfficiencyVsPt_GedGsfElec->SetMarkerStyle(24);
        electron->EfficiencyVsPt_PFElec->SetMarkerStyle(25);
	electron->EfficiencyVsPt_Seed->SetMarkerColor(1);
        electron->EfficiencyVsPt_GedGsfElec->SetMarkerColor(1);	
        electron->EfficiencyVsPt_PFElec->SetMarkerColor(1);

	electron_lowPU->EfficiencyVsPt_Seed->SetMarkerStyle(22);
	electron_lowPU->EfficiencyVsPt_Seed->SetFillColor(2);
        electron_lowPU->EfficiencyVsPt_Seed->SetFillStyle(3005);
	electron_lowPU->EfficiencyVsPt_GedGsfElec->SetMarkerStyle(24);	
        electron_lowPU->EfficiencyVsPt_PFElec->SetMarkerStyle(25);
	electron_lowPU->EfficiencyVsPt_Seed->SetMarkerColor(2);
        electron_lowPU->EfficiencyVsPt_GedGsfElec->SetMarkerColor(2);
        electron_lowPU->EfficiencyVsPt_PFElec->SetMarkerColor(2);

	electron_highPU->EfficiencyVsPt_Seed->SetMarkerStyle(22);
	electron_highPU->EfficiencyVsPt_Seed->SetFillColor(4);
        electron_highPU->EfficiencyVsPt_Seed->SetFillStyle(3003);

        electron_highPU->EfficiencyVsPt_GedGsfElec->SetMarkerStyle(24);
        electron_highPU->EfficiencyVsPt_PFElec->SetMarkerStyle(25);
	electron_highPU->EfficiencyVsPt_Seed->SetMarkerColor(4);
        electron_highPU->EfficiencyVsPt_GedGsfElec->SetMarkerColor(3);
        electron_highPU->EfficiencyVsPt_PFElec->SetMarkerColor(3);

	electron_lowPU->EfficiencyVsPt_Seed->Draw("pe2same");
//        electron_lowPU->EfficiencyVsPt_GedGsfElec->Draw("pe1same");
//        electron_lowPU->EfficiencyVsPt_PFElec->Draw("pe1same");
	electron_highPU->EfficiencyVsPt_Seed->Draw("pe2same");
//        electron_highPU->EfficiencyVsPt_GedGsfElec->Draw("pe1same");	
//        electron_highPU->EfficiencyVsPt_PFElec->Draw("pe1same");

	legend->Draw();	

	can->Print(ProcessName+"_Electron_RecoEff_vs_Pt.pdf");	

	electron->EfficiencyVsEta_Seed->SetMaximum(1.2);
	electron->EfficiencyVsEta_Seed->SetMinimum(0.0001);	
	electron->EfficiencyVsEta_Seed->Draw("pe1");
        electron->EfficiencyVsEta_GedGsfElec->Draw("pe1same");
        electron->EfficiencyVsEta_PFElec->Draw("pe1same");


        electron->EfficiencyVsEta_Seed->SetMarkerColor(1);
        electron->EfficiencyVsEta_GedGsfElec->SetMarkerColor(1);
        electron->EfficiencyVsEta_PFElec->SetMarkerColor(1);
        electron->EfficiencyVsEta_Seed->SetMarkerStyle(22);
        electron->EfficiencyVsEta_GedGsfElec->SetMarkerStyle(24);
        electron->EfficiencyVsEta_PFElec->SetMarkerStyle(25);

        electron_lowPU->EfficiencyVsEta_Seed->SetMarkerColor(2);
        electron_lowPU->EfficiencyVsEta_GedGsfElec->SetMarkerColor(2);
        electron_lowPU->EfficiencyVsEta_PFElec->SetMarkerColor(2);
        electron_highPU->EfficiencyVsEta_Seed->SetMarkerColor(3);
        electron_highPU->EfficiencyVsEta_GedGsfElec->SetMarkerColor(3);
        electron_highPU->EfficiencyVsEta_PFElec->SetMarkerColor(3);

	electron_lowPU->EfficiencyVsEta_Seed->SetMarkerStyle(22);
        electron_lowPU->EfficiencyVsEta_GedGsfElec->SetMarkerStyle(24);
        electron_lowPU->EfficiencyVsEta_PFElec->SetMarkerStyle(25);
        electron_highPU->EfficiencyVsEta_Seed->SetMarkerStyle(22);
        electron_highPU->EfficiencyVsEta_GedGsfElec->SetMarkerStyle(24);
        electron_highPU->EfficiencyVsEta_PFElec->SetMarkerStyle(25);
/*
        electron_lowPU->EfficiencyVsEta_Seed->Draw("pe1same");
        electron_lowPU->EfficiencyVsEta_GedGsfElec->Draw("pe1same");
        electron_lowPU->EfficiencyVsEta_PFElec->Draw("pe1same");
        electron_highPU->EfficiencyVsEta_Seed->Draw("pe1same");
        electron_highPU->EfficiencyVsEta_GedGsfElec->Draw("pe1same");
        electron_highPU->EfficiencyVsEta_PFElec->Draw("pe1same");
*/
	legend->Draw();
	
        can->Print(ProcessName+"_Electron_RecoEff_vs_Eta.pdf");
	//-------------------------------------------------------------------

	gPad->SetLogy(1);	
	legend->Clear();
	
	pion->EfficiencyVsPt_Seed->SetMaximum(1.2);
	pion->EfficiencyVsPt_Seed->SetMinimum(0.001);
        pion->EfficiencyVsPt_Seed->Draw("pe2");
	pion_lowPU->EfficiencyVsEta_Seed->SetMarkerColor(2);
        pion_highPU->EfficiencyVsEta_Seed->SetMarkerColor(3);
	pion_lowPU->EfficiencyVsEta_Seed->SetMarkerStyle(22);
	pion_highPU->EfficiencyVsEta_Seed->SetMarkerStyle(22);
	pion_lowPU->EfficiencyVsPt_Seed->SetFillColor(2);
        pion_lowPU->EfficiencyVsPt_Seed->SetFillStyle(3005);
	pion_highPU->EfficiencyVsPt_Seed->SetFillColor(4);
        pion_highPU->EfficiencyVsPt_Seed->SetFillStyle(3003);
	pion_lowPU->EfficiencyVsPt_Seed->Draw("pe2same");
	pion_highPU->EfficiencyVsPt_Seed->Draw("pe2same");
	


        pion->EfficiencyVsPt_GedGsfElec->Draw("pe2same");
  

        pion->EfficiencyVsPt_PFElec->Draw("pe2same");
	pion->EfficiencyVsPt_Seed->SetMarkerStyle(22);
        pion->EfficiencyVsPt_GedGsfElec->SetMarkerStyle(22);
        pion->EfficiencyVsPt_PFElec->SetMarkerStyle(24);
	pion->EfficiencyVsPt_Seed->SetMarkerColor(1);
        pion->EfficiencyVsPt_GedGsfElec->SetMarkerColor(1);
        pion->EfficiencyVsPt_PFElec->SetMarkerColor(1);
	legend->AddEntry(pion->EfficiencyVsPt_Seed,"#pi#rightarrow Seed (all PU)");
	legend->AddEntry(pion->EfficiencyVsPt_GedGsfElec,"#pi#rightarrow Seed (<30)");
	legend->AddEntry(pion->EfficiencyVsPt_PFElec,"#pi#rightarrow Seed (>30)");

        legend->Draw();

        can->Print(ProcessName+"_Pion_RecoEff_vs_Pt.pdf");

	pion->EfficiencyVsEta_Seed->SetMarkerStyle(22);
        pion->EfficiencyVsEta_GedGsfElec->SetMarkerStyle(23);
        pion->EfficiencyVsEta_PFElec->SetMarkerStyle(24);
        pion->EfficiencyVsEta_Seed->SetMarkerColor(1);
        pion->EfficiencyVsEta_GedGsfElec->SetMarkerColor(1);
        pion->EfficiencyVsEta_PFElec->SetMarkerColor(1);

        pion->EfficiencyVsEta_Seed->SetMaximum(1.2);
        pion->EfficiencyVsEta_Seed->SetMinimum(0.001);
        pion->EfficiencyVsEta_Seed->Draw("pe1");
        pion->EfficiencyVsEta_GedGsfElec->Draw("pe1same");	
        pion->EfficiencyVsEta_PFElec->Draw("pe1same");
	legend->Draw();
	can->Print(ProcessName+"_Pion_RecoEff_vs_Eta.pdf");

	legend->Clear();

	delete legend;
	delete can;
}

void SLPlotter::MakeEffvsEffPlot(TString &process){
	setTDRStyle();
	TCanvas *can = new TCanvas("h", "bla",600,600);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.05);
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);
	ProcessName=process;
	float x_ged_lpu[Nbin2D],y_ged_lpu[Nbin2D];
	float x_ged_mpu[Nbin2D],y_ged_mpu[Nbin2D];
	float x_ged_hpu[Nbin2D],y_ged_hpu[Nbin2D];
        float x_ged_all[Nbin2D],y_ged_all[Nbin2D];
        float x_ged_all_effvsmistag[Nbin2D],y_ged_all_effvsmistag[Nbin2D];
	float x_pf_lpu[Nbin2D],y_pf_lpu[Nbin2D];
        float x_pf_mpu[Nbin2D],y_pf_mpu[Nbin2D];
        float x_pf_hpu[Nbin2D],y_pf_hpu[Nbin2D];
        float x_pf_all[Nbin2D],y_pf_all[Nbin2D];
	cout<<"elpu ForIntegral: "<<electron_lowPU->ForIntegral->Integral()<<endl;
        cout<<"plpu ForIntegral: "<<pion_lowPU->ForIntegral->Integral()<<endl;

        for(int k=0;k<Nbin2D;k++){
                float fr1=electron->mva_e_pi->Integral(k+1,Nbin2D)/electron->Pt_gen_all->Integral();
                float fr2=pion->mva_e_pi->Integral(k+1,Nbin2D)/pion->Pt_gen_all->Integral();
                x_ged_all[k]=fr1;
                y_ged_all[k]=fr2;
        }
        TGraph *gr_ged_all=new TGraph(Nbin2D,x_ged_all,y_ged_all);

        for(int k=0;k<Nbin2D;k++){
                float fr1=electron->mva_e_pi->Integral(k+1,Nbin2D)/electron->Pt_gen_all->Integral();
                float fr2=background->mva_e_pi->Integral(k+1,Nbin2D)/background->mva_e_pi->Integral();
                x_ged_all_effvsmistag[k]=fr1;
                y_ged_all_effvsmistag[k]=fr2;
        }
        TGraph *gr_ged_all_effvsmistag=new TGraph(Nbin2D,x_ged_all_effvsmistag,y_ged_all_effvsmistag);


        for(int k=0;k<Nbin2D;k++){
		float fr1=electron_lowPU->mva_e_pi->Integral(k+1,Nbin2D)/electron_lowPU->Pt_gen_all->Integral();
		float fr2=pion_lowPU->mva_e_pi->Integral(k+1,Nbin2D)/pion_lowPU->Pt_gen_all->Integral();
                x_ged_lpu[k]=fr1;
                y_ged_lpu[k]=fr2;
		cout<<"lowpu: bin"<<k<<" center of bin "<<electron_lowPU->mva_e_pi->GetBinCenter(k+1)<<" "<<electron_lowPU->mva_e_pi->Integral(k+1,Nbin2D)<<" "<<electron_lowPU->Pt_gen_all->Integral()<<" "<<fr1<<" "<<fr2<<endl;
        }
	TGraph *gr_ged_lpu=new TGraph(Nbin2D,x_ged_lpu,y_ged_lpu);
	for(int k=0;k<Nbin2D;k++){
		float fr1=electron_midPU->mva_e_pi->Integral(k+1,Nbin2D)/electron_midPU->Pt_gen_all->Integral();
                float fr2=pion_midPU->mva_e_pi->Integral(k+1,Nbin2D)/pion_midPU->Pt_gen_all->Integral();
                x_ged_mpu[k]=fr1;
                y_ged_mpu[k]=fr2;
		cout<<"midpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_ged_mpu=new TGraph(Nbin2D,x_ged_mpu,y_ged_mpu);
	for(int k=0;k<Nbin2D;k++){
		float fr1=electron_highPU->mva_e_pi->Integral(k+1,Nbin2D)/electron_highPU->Pt_gen_all->Integral();
                float fr2=pion_highPU->mva_e_pi->Integral(k+1,Nbin2D)/pion_highPU->Pt_gen_all->Integral();
                x_ged_hpu[k]=fr1;
                y_ged_hpu[k]=fr2;
		cout<<"highpu:"<<fr1<<" "<<fr2<<endl;
        }
	TGraph *gr_ged_hpu=new TGraph(Nbin2D,x_ged_hpu,y_ged_hpu);

/////////////////
        for(int k=0;k<Nbin2D;k++){
                float fr1=electron_lowPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron_lowPU->Pt_gen_all->Integral();
                float fr2=pion_lowPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion_lowPU->Pt_gen_all->Integral();
                x_pf_lpu[k]=fr1;
                y_pf_lpu[k]=fr2;
                cout<<"lowpu: bin"<<k<<" center of bin "<<electron_lowPU->mva_e_pi->GetBinCenter(k+1)<<" "<<electron_lowPU->mva_e_pi->Integral(k+1,Nbin2D)<<" "<<electron_lowPU->Pt_gen_all->Integral()<<" "<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_lpu=new TGraph(Nbin2D,x_pf_lpu,y_pf_lpu);
        for(int k=0;k<Nbin2D;k++){
                float fr1=electron_midPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron_midPU->Pt_gen_all->Integral();
                float fr2=pion_midPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion_midPU->Pt_gen_all->Integral();
                x_pf_mpu[k]=fr1;
                y_pf_mpu[k]=fr2;
                cout<<"midpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_mpu=new TGraph(Nbin2D,x_pf_mpu,y_pf_mpu);
        for(int k=0;k<Nbin2D;k++){
                float fr1=electron_highPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron_highPU->Pt_gen_all->Integral();
                float fr2=pion_highPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion_highPU->Pt_gen_all->Integral();
                x_pf_hpu[k]=fr1;
                y_pf_hpu[k]=fr2;
                cout<<"highpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_hpu=new TGraph(Nbin2D,x_pf_hpu,y_pf_hpu);

        for(int k=0;k<Nbin2D;k++){
                float fr1=electron->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron->Pt_gen_all->Integral();
                float fr2=pion->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion->Pt_gen_all->Integral();
                x_pf_all[k]=fr1;
                y_pf_all[k]=fr2;
                cout<<"highpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_all=new TGraph(Nbin2D,x_pf_all,y_pf_all);



	gr_ged_lpu->SetMarkerStyle(21);
	gr_ged_mpu->SetMarkerStyle(22);
	gr_ged_hpu->SetMarkerStyle(23);
	gr_ged_lpu->SetMarkerColor(kRed);
        gr_ged_mpu->SetMarkerColor(kBlue);
	gr_ged_lpu->GetXaxis()->SetLimits(0.0,1.0);
	gr_ged_lpu->SetMaximum(1.0);
	gr_ged_lpu->SetMinimum(0.00001);
	gr_ged_lpu->GetXaxis()->SetLimits(0.05,0.7);
	gr_ged_mpu->GetXaxis()->SetLimits(0.05,0.7);
	gr_ged_hpu->GetXaxis()->SetLimits(0.05,0.7);
	gr_ged_lpu->GetXaxis()->SetTitle("e#rightarrow gedgsf el efficiency");
	gr_ged_lpu->GetYaxis()->SetTitle("#pi#rightarrow gedgsf el efficiency");

        gr_pf_lpu->SetMarkerStyle(25);
        gr_pf_mpu->SetMarkerStyle(26);
        gr_pf_hpu->SetMarkerStyle(32);
        gr_pf_lpu->SetMarkerColor(kRed);
        gr_pf_mpu->SetMarkerColor(kBlue);
        gr_pf_lpu->GetXaxis()->SetLimits(0.0,1.0);
        gr_pf_lpu->SetMaximum(1.0);
        gr_pf_lpu->SetMinimum(0.00001);
        gr_pf_lpu->GetXaxis()->SetLimits(0.05,0.7);
        gr_pf_mpu->GetXaxis()->SetLimits(0.05,0.7);
        gr_pf_hpu->GetXaxis()->SetLimits(0.05,0.7);
        gr_pf_lpu->GetXaxis()->SetTitle("e#rightarrow gedgsf el efficiency");
        gr_pf_lpu->GetYaxis()->SetTitle("#pi#rightarrow gedgsf el efficiency");

	
	gr_ged_lpu->Draw("Ap");
	gr_ged_mpu->Draw("psame");
	gr_ged_hpu->Draw("psame");

//        gr_pf_lpu->Draw("psame");
//        gr_pf_mpu->Draw("psame");
//        gr_pf_hpu->Draw("psame");

	gPad->SetLogy(1);
	gPad->SetGridy(1);
        gPad->SetGridx(1);
	TLegend *legend=new TLegend(0.65,0.8,0.95,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
	legend->AddEntry(gr_ged_lpu,"GED PU<10","p");
	legend->AddEntry(gr_ged_mpu,"GED 10<=PU<25","p");
	legend->AddEntry(gr_ged_hpu,"GED PU>25","p");
//        legend->AddEntry(gr_ged_lpu,"PF PU<10","p");
//        legend->AddEntry(gr_ged_mpu,"PF 10<=PU<25","p");
//        legend->AddEntry(gr_ged_hpu,"PF PU>25","p");
	legend->Draw();
	can->Print(ProcessName+"effvseff.pdf");	

	gr_ged_all->SetMarkerStyle(22);
	gr_pf_all->SetMarkerStyle(25);
	gr_ged_all->Draw("Ap");
        gr_pf_all->Draw("psame");

	gr_ged_all->SetMaximum(1.0);
        gr_pf_all->SetMinimum(0.00001);
	gr_ged_all->GetXaxis()->SetLimits(0.05,0.7);
        gr_pf_all->GetXaxis()->SetLimits(0.05,0.7);
	gr_ged_all->GetXaxis()->SetTitle("e#rightarrow gedgsf el efficiency");
        gr_ged_all->GetYaxis()->SetTitle("#pi#rightarrow gedgsf el efficiency");	

	legend->Clear();
	legend->AddEntry(gr_ged_all,"GED new mva_e_pi","p");
        legend->AddEntry(gr_pf_all,"PF  mva_e_pi","p");
	legend->Draw();
        can->Print(ProcessName+"effvseff_GEDvsPF.pdf");
		

	gr_ged_all_effvsmistag->SetMarkerStyle(22);
	gr_ged_all_effvsmistag->Draw("Ap");
	gr_ged_all_effvsmistag->SetMaximum(1.0);
        gr_ged_all_effvsmistag->SetMinimum(0.00001);
        gr_ged_all_effvsmistag->GetXaxis()->SetLimits(0.05,0.7);
        gr_ged_all_effvsmistag->GetXaxis()->SetLimits(0.05,0.7);
	legend->Clear();
        legend->AddEntry(gr_ged_all,"GED new mva_e_pi","p");
        legend->Draw();
	gr_ged_all_effvsmistag->GetXaxis()->SetTitle("e#rightarrow gedgsf el efficiency");
	gr_ged_all_effvsmistag->GetYaxis()->SetTitle("mistag rate");
        can->Print(ProcessName+"effvsmistagerate_GED.pdf");


	
}


void SLPlotter::TheFill(HistoMaker *obj,EventElecComm *Evt){
        int flagEta=fabs(Evt->etaGen)<1.4?0:1;

        obj->Pt_gen_all->Fill(Evt->ptGen);
        obj->Eta_gen_all->Fill(Evt->etaGen);
        if(Evt->isMatchedWithASeed==1){
	      	obj->Pt_gen_MatchedWithASeed->Fill(Evt->ptGen);
                obj->Eta_gen_MatchedWithASeed->Fill(Evt->etaGen);
        }
	if(Evt->isMatchedWithAGedGsfElec==1 && Evt->mva_e_pi>-0.1){
        	obj->Pt_gen_MatchedWithAGedGsfElec->Fill(Evt->ptGen);
                obj->Eta_gen_MatchedWithAGedGsfElec->Fill(Evt->etaGen);
		obj->mva_e_pi->Fill(Evt->mva_e_pi);
	}
	if(Evt->isMatchedWithAPFElec==1 && Evt->mva_e_pi_PF>-0.1){
                obj->Pt_gen_MatchedWithAPFElec->Fill(Evt->ptGen);
                obj->Eta_gen_MatchedWithAPFElec->Fill(Evt->etaGen);
                obj->mva_e_pi_PF->Fill(Evt->mva_e_pi_PF);
        }		
	//cout<<"histo integral "<<obj->Pt_gen_MatchedWithASeed->Integral()<<" "<<obj->Eta_gen_MatchedWithASeed->Integral()<<endl;	
}


void SLPlotter::FillHistos(){
	for(int ientry=0; ientry<nentries; ientry++) {
		if(ientry % 20000==0)cout<<"Event number "<<ientry<<endl;
                Evt->GetEntry(ientry);
		bool basicselection= Evt->ptGen>2 && fabs(Evt->etaGen)<2.4 ;
        	bool isfromB=(Evt->origin==4 || Evt->origin==6);
        	bool isfromD=(Evt->origin==2);
        	bool isfromX=(Evt->origin==0);
	//	if(Evt->ptGen>2 && fabs(Evt->etaGen)<2.4)cout<<Evt->ptGen<<" "<<Evt->etaGen<<" "<<isfromB<<" "<<abs(Evt->pdgId)<<endl;
		if(basicselection && isfromB && abs(Evt->pdgId)==11){
	//		cout<<"an electron spotted"<<endl;
			TheFill(electron,Evt);
			if(Evt->nPV<30){
				TheFill(electron_lowPU,Evt);
			}
			if(Evt->nPV>=25 && Evt->nPV<35){
				TheFill(electron_midPU,Evt);
			}
			if(Evt->nPV>=30){
				TheFill(electron_highPU,Evt);
			}
                }
		if(basicselection && abs(Evt->pdgId)==211){
			TheFill(pion,Evt);
			if(Evt->nPV<30){
				TheFill(pion_lowPU,Evt);
			}
                        if(Evt->nPV>=25 && Evt->nPV<35){
				TheFill(pion_midPU,Evt);
			}
                        if(Evt->nPV>=30){
				TheFill(pion_highPU,Evt);
			}
                }
		if(basicselection && abs(Evt->pdgId)==311){
                        TheFill(kaon,Evt);
                        if(Evt->nPV<30){
                                TheFill(kaon_lowPU,Evt);
                        }
                        if(Evt->nPV>=25 && Evt->nPV<35){
                                TheFill(kaon_midPU,Evt);
                        }
                        if(Evt->nPV>=30){
                                TheFill(kaon_highPU,Evt);
                        }
                }
		if(basicselection && abs(Evt->pdgId)!=11){
			TheFill(background,Evt);
		}

        }
	electron->Combine();
	electron_lowPU->Combine();
	electron_highPU->Combine();
	pion->Combine();
	kaon->Combine();


}


void SLPlotter::setTDRStyle() {
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetCanvasColor(kWhite);
        gStyle->SetCanvasDefH(1500); //Height of canvas
        gStyle->SetCanvasDefW(1500); //Width of canvas
        gStyle->SetCanvasDefX(0);   //POsition on screen
        gStyle->SetCanvasDefY(0);

        gStyle->SetPadBorderMode(0);
        gStyle->SetPadColor(kWhite);
        gStyle->SetPadGridX(false);
        gStyle->SetPadGridY(false);
        gStyle->SetGridColor(0);
        gStyle->SetGridStyle(3);
        gStyle->SetGridWidth(1);

        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameBorderSize(0.1);
        gStyle->SetFrameFillColor(0);
        gStyle->SetFrameFillStyle(0);
        gStyle->SetFrameLineColor(1);
        gStyle->SetFrameLineStyle(1);
        gStyle->SetFrameLineWidth(0.1);

        gStyle->SetHistLineColor(1);
        gStyle->SetHistLineStyle(0);
        gStyle->SetHistLineWidth(0.1);

        gStyle->SetEndErrorSize(2);
//        gStyle->SetErrorX(0.);
        gStyle->SetMarkerStyle(20);

        gStyle->SetOptFit(1);
        gStyle->SetFitFormat("5.4g");
        gStyle->SetFuncColor(2);
        gStyle->SetFuncStyle(1);
        gStyle->SetFuncWidth(1);

        gStyle->SetOptDate(0);
        gStyle->SetOptStat(0);

        // Margins:
        gStyle->SetPadTopMargin(0.05);
        gStyle->SetPadBottomMargin(0.13);
        gStyle->SetPadLeftMargin(0.16);
        gStyle->SetPadRightMargin(0.02);

        // For the Global title:

        gStyle->SetOptTitle(0);
        gStyle->SetTitleFont(42);
        gStyle->SetTitleColor(1);
        gStyle->SetTitleTextColor(1);
        gStyle->SetTitleFillColor(10);
        gStyle->SetTitleFontSize(0.05);

        // For the axis titles:

        gStyle->SetTitleColor(1, "XYZ");
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetTitleSize(0.06, "XYZ");
        // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
        // gStyle->SetTitleYSize(Float_t size = 0.02);
        gStyle->SetTitleXOffset(0.9);
        gStyle->SetTitleYOffset(1.25);
        // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

        // For the axis labels:

        gStyle->SetLabelColor(1, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetLabelOffset(0.007, "XYZ");
        gStyle->SetLabelSize(0.05, "XYZ");

        // For the axis:

        gStyle->SetAxisColor(1, "XYZ");
        gStyle->SetStripDecimals(kTRUE);
        gStyle->SetTickLength(0.03, "XYZ");
        gStyle->SetNdivisions(510, "XYZ");
        gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
        gStyle->SetPadTickY(1);

        // Change for log plots:
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        gStyle->SetOptLogz(0);

        gROOT->ForceStyle();

}

class PlotTogether
{
        public:
                PlotTogether();
                ~PlotTogether();
		void EffPlotCombination();
};

PlotTogether::PlotTogether()
{

}

PlotTogether::~PlotTogether()
{
}

void PlotTogether::EffPlotCombination(){
	SLPlotter *plot[3];
	
	plot[0]	=	new SLPlotter ("input_efficiency_TTbar.txt","TTbar");
	plot[1]	=	new SLPlotter("input_efficiency_30_80.txt","QCD30_80");
	plot[2]	=	new SLPlotter("input_efficiency_80_170.txt","QCD80_170");
	
	TH1F *SeedEffPt[3];
	TH1F *SeedEffEta[3];
	TH1F *GEDEffPt[3];
        TH1F *GEDEffEta[3];

	TLegend *legend=new TLegend(0.65,0.8,0.95,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
	for(int u=0;u<3;u++){
		plot[u]->setTDRStyle();
		plot[u]->initialize();
		plot[u]->FillHistos();
		SeedEffPt[u]=plot[u]->electron->EfficiencyVsPt_Seed;
       		SeedEffEta[u]=plot[u]->electron->EfficiencyVsEta_Seed;
		GEDEffPt[u]=plot[u]->electron->EfficiencyVsPt_GedGsfElec;
                GEDEffEta[u]=plot[u]->electron->EfficiencyVsEta_GedGsfElec;	
		SeedEffPt[u]->SetMarkerColor(1+u);
                SeedEffEta[u]->SetMarkerColor(1+u);
                GEDEffPt[u]->SetMarkerColor(1+u);
                GEDEffEta[u]->SetMarkerColor(1+u);
		SeedEffPt[u]->SetMarkerStyle(22+u);
                SeedEffEta[u]->SetMarkerStyle(24+u);
                GEDEffPt[u]->SetMarkerStyle(24+u);
                GEDEffEta[u]->SetMarkerStyle(24+u);
		SeedEffPt[u]->SetMaximum(1.2);
        	SeedEffPt[u]->SetMinimum(0.0);
	}

        TCanvas *can = new TCanvas("h", "bla",600,600);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.05);
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);

	gPad->SetGridx();
        gPad->SetGridy();
	cout<<"getting the plots"<<endl;
        legend->AddEntry(SeedEffPt[0],"Seeding TTbar","p");
	legend->AddEntry(GEDEffPt[0],"GED TTbar","p");
	legend->AddEntry(SeedEffPt[1],"Seeding QCD-30-80","p");
        legend->AddEntry(GEDEffPt[1],"GED QCD-30-80","p");
        legend->AddEntry(SeedEffPt[2],"Seeding QCD-80-170","p");
        legend->AddEntry(GEDEffPt[2],"GED QCD-80-170","p");

	SeedEffPt[0]->Draw("p");
	SeedEffPt[1]->Draw("psame");
	SeedEffPt[2]->Draw("psame");
	GEDEffPt[0]->Draw("psame");
        GEDEffPt[1]->Draw("psame");
        GEDEffPt[2]->Draw("psame");
	legend->Draw();
	can->Print("SeedGED_Efficiencies_pt_comp.pdf");


	
	SeedEffEta[0]->Draw("p");
      	SeedEffEta[1]->Draw("psame");
        SeedEffEta[2]->Draw("psame");
        GEDEffEta[0]->Draw("psame");
        GEDEffEta[1]->Draw("psame");
        GEDEffEta[2]->Draw("psame");
	legend->Draw();
        can->Print("SeedGED_Efficiencies_eta_comp.pdf");

}

