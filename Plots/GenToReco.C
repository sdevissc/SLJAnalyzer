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

template <class Type>
void SetStyle(Type *a,int color, int style, int fillstyle, const TString& xtitle,const TString& ytitle, float min, float max){
        a->SetLineColor(color);
        a->SetMarkerColor(color);
        a->SetMarkerStyle(style);
        a->SetFillColor(color);
        a->SetFillStyle(fillstyle);
        a->GetXaxis()->SetTitle(xtitle);
        a->GetYaxis()->SetTitle(ytitle);
        if(min!=-777.0 && max!=-777.0 ){
                a->SetMinimum(min);
                a->SetMaximum(max);
        }
}


class HistoMaker
{
  public:
	HistoMaker(const TString& ,const TString &);
	~HistoMaker();
	void Combine();
	
	TH1F* TemplatePt;
	TH1F* TemplateEta;
	TH1F *Pt[6],*EffPt[6],*EffEta[6],*Eta[6];
/*
	TH1F *Pt[0];
	TH1F *Pt_gen_MatchedWithATrackSeed;
	TH1F *Pt_gen_MatchedWithAECALSeed;
	TH1F *Pt_gen_MatchedWithAECAL_or_TrackerSeed;
	TH1F *Pt_gen_MatchedWithAGedGsfElec;
        TH1F *Pt_gen_MatchedWithAPFElec;
        
        TH1F *Eta_gen_all;
	TH1F *Eta_gen_MatchedWithATrackSeed;
       	TH1F *Eta_gen_MatchedWithAECALSeed;
        TH1F *Eta_gen_MatchedWithAECAL_or_TrackerSeed;
	TH1F *Eta_gen_MatchedWithAGedGsfElec;
        TH1F *Eta_gen_MatchedWithAPFElec;

	TH1F *EffPt[1];
        TH1F *EffEta[1];
	TH1F *EfficiencyVsPt_ECALSeed;
        TH1F *EfficiencyVsEta_ECALSeed;
	TH1F *EfficiencyVsPt_ECAL_or_TrackerSeed;
        TH1F *EfficiencyVsEta_ECAL_or_TrackerSeed;
        TH1F *EffPt[4];
        TH1F *EffEta[4];
        TH1F *EfficiencyVsPt_PFElec;
        TH1F *EfficiencyVsEta_PFElec;
*/	
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
	float etax[14]={-2.5 , -2.1  , -1.7 , -1.3 , -0.9 , -0.5 , -0.1 , 0.1 , 0.5 , 0.9 , 1.3 , 1.7 , 2.1 , 2.5 };
        float ptx[nbinpt+1]= { 2.0 , 5.0 , 10.0 , 20.0 , 30.0 , 50.0 , 120.0};

	cout<<"define pt and eta basic histos"<<endl;
	ForIntegral=new TH1F(inputname+tag+"forIntegral","",Nbin2D,-100000,100000);	
	TemplatePt=new TH1F(inputname+tag+"Pt[0]","",nbinpt,ptx);
	TemplateEta=new TH1F(inputname+tag+"Eta[0]","",nbineta,etax);
	
	for(int i=0;i<6;i++){
		Pt[i]=(TH1F*)TemplatePt->Clone(Form("Pt_%d",i));
		Eta[i]=(TH1F*)TemplateEta->Clone(Form("Eta_%d",i));
		EffPt[i]=(TH1F*)TemplatePt->Clone(Form("EffPt_%d",i));
                EffEta[i]=(TH1F*)TemplateEta->Clone(Form("EffEta_%d",i));
	}

}

HistoMaker::~HistoMaker()                 // destructor, just an example
{
}





void HistoMaker::Combine(){
	TH1::SetDefaultSumw2(kTRUE);

	for(int i=1;i<6;i++){
		EffPt[i]->Add(Pt[i]);
		EffPt[i]->Divide(Pt[0]);
		EffEta[i]->Add(Eta[i]);
		EffEta[i]->Divide(Eta[0]);
	}
	cout<<EffPt[1]->Integral()<<" "<<EffEta[1]->Integral()<<" "<<EffPt[4]->Integral()<<" "<<EffEta[4]->Integral()<<endl;
	
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
	HistoMaker *electron_isolated;
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
	electron_isolated=new HistoMaker("elecisol_hpu",TheTag);
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

//---------------------------------------------------------	
        legend->SetTextSize(0.020);
	legend->AddEntry(electron->EffPt[1],"e_{B}#rightarrow TDS (all PU)");
        legend->AddEntry(electron_lowPU->EffPt[1],"e_{B}#rightarrow TDS (<30 PU)");
        legend->AddEntry(electron_highPU->EffPt[1],"e_{B}#rightarrow TDS (>30 PU)");
	legend->AddEntry(pion->EffPt[1],"#pi#rightarrow TDS (all PU)");
        legend->AddEntry(pion_lowPU->EffPt[1],"#pi#rightarrow TDS (<30 PU)");
        legend->AddEntry(pion_highPU->EffPt[1],"#pi#rightarrow TDS (>30 PU)");
	

	gPad->SetGridx();
        gPad->SetGridy();
	gPad->SetLogy(log);

	SetStyle(electron->EffPt[1],1,20,3003,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
	SetStyle(electron_lowPU->EffPt[1],2,21,3004,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
	SetStyle(electron_highPU->EffPt[1],4,22,3005,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
	SetStyle(pion->EffPt[1],1,23,3006,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(pion_lowPU->EffPt[1],2,33,3002,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(pion_highPU->EffPt[1],4,34,3017,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);

	electron->EffPt[1]->Draw("pe2");
	electron_lowPU->EffPt[1]->Draw("pe2same");
	electron_highPU->EffPt[1]->Draw("pe2same");
	pion->EffPt[1]->Draw("pe2same");
        pion_lowPU->EffPt[1]->Draw("pe2same");
        pion_highPU->EffPt[1]->Draw("pe2same");


	legend->Draw();	

	can->Print(Form(ProcessName+"_Electron_and_pion_TDS_SeedEff_vs_Pt_log%d.pdf",log));	
	legend->Clear();
//----------------------------------------------------------------


        legend->SetTextSize(0.020);
        legend->AddEntry(electron->EffPt[1],"e_{B}#rightarrow TDS");
        legend->AddEntry(electron_lowPU->EffPt[1],"e_{B}#rightarrow EDS ");
        legend->AddEntry(electron_highPU->EffPt[1],"e_{B}#rightarrow TDS #wedge EDS");
        legend->AddEntry(pion->EffPt[1],"#pi#rightarrow TDS");
        legend->AddEntry(pion_lowPU->EffPt[1],"#pi#rightarrow EDS");
        legend->AddEntry(pion_highPU->EffPt[1],"#pi#rightarrow TDS #wedge EDS");


        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogy(log);

        SetStyle(electron->EffPt[1],1,20,3003,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(electron->EffPt[2],2,21,3004,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(electron->EffPt[3],4,22,3005,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(pion->EffPt[1],1,23,3006,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(pion->EffPt[2],2,33,3002,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);
        SetStyle(pion->EffPt[3],4,34,3017,"p_{T} of the generated particle","Efficiency",0.02,log==0?1.2:10);

        electron->EffPt[1]->Draw("pe2");
        electron->EffPt[2]->Draw("pe2same");
        electron->EffPt[3]->Draw("pe2same");
        pion->EffPt[1]->Draw("pe2same");
        pion->EffPt[2]->Draw("pe2same");
        pion->EffPt[3]->Draw("pe2same");


        legend->Draw();

        can->Print(Form(ProcessName+"_Electron_and_pion_TDSvsEDSvsTDSOREDS_seedeff_vs_Pt_log%d.pdf",log));
        legend->Clear();

//------------------------------------------------------------------------	


        legend->SetTextSize(0.020);
        //ilegend->AddEntry(electron_isolated->EffPt[1],"e_{V}#rightarrow TDS");
        //legend->AddEntry(electron_isolated->EffPt[2],"e_{V}#rightarrow EDS");
        //legend->AddEntry(electron_isolated->EffPt[3],"e_{V}#rightarrow  TDS #wedge EDS");

	
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogy(log);

        SetStyle(electron_isolated->EffPt[1],1,20,3003,"p_{T} of the generated particle","Efficiency",0.99,550);
        SetStyle(electron_isolated->EffPt[2],2,21,3004,"p_{T} of the generated particle","Efficiency",0.99,550);
        SetStyle(electron_isolated->EffPt[3],4,22,3005,"p_{T} of the generated particle","Efficiency(#%)",0.99,550);
	
	electron_isolated->EffPt[3]->GetYaxis()->SetNdivisions(520);
	electron_isolated->EffPt[3]->Scale(100);
//        electron_isolated->EffPt[1]->Draw("pe2");
//        electron_isolated->EffPt[2]->Draw("pe2same");
//        electron_isolated->EffPt[3]->Draw("pe2same");
	TH1F*h =(TH1F*)electron_isolated->EffPt[3]->Clone("improvement");
	h->Divide(electron_isolated->EffPt[2]);
	SetStyle(h,6,23,3005,"p_{T} of the generated electron","Efficiency improvement (%)",50,1000);
	legend->AddEntry(h,"#varepsilon_{TDS #wedge EDS}/#varepsilon_{EDS}");
	h->Draw("p");


        legend->Draw();

        can->Print(Form(ProcessName+"_IsolatedElectron_TDSvsEDSvsTDSOREDS_seedeff_vs_Pt_log%d.pdf",log));
        legend->Clear();

//----------------------------------------------------------

	SetStyle(electron->EffEta[1],1,22,3004,"#eta_{T} of the generated particle","",0.001,1.2);
        SetStyle(electron_lowPU->EffEta[1],2,23,3005,"#eta_{T} of the generated particle","",0.001,1.2);
        SetStyle(electron_highPU->EffEta[1],3,24,3006,"#eta_{T} of the generated particle","",0.001,1.2);

	electron->EffEta[1]->Draw("hpe1");
        electron->EffEta[4]->Draw("hpe1same");
        electron->EffEta[5]->Draw("hpe1same");

	legend->Draw();
	
        can->Print(ProcessName+"_Electron_RecoEff_vs_Eta.pdf");
	//-------------------------------------------------------------------

	gPad->SetLogy(1);	
	legend->Clear();
	
	pion->EffPt[1]->SetMaximum(1.2);
	pion->EffPt[1]->SetMinimum(0.001);
        pion->EffPt[1]->Draw("pe2");
	pion_lowPU->EffPt[1]->Draw("pe2same");
	pion_highPU->EffPt[1]->Draw("pe2same");
	


        pion->EffPt[4]->Draw("pe2same");
  

        pion->EffPt[5]->Draw("pe2same");
	legend->AddEntry(pion->EffPt[1],"#pi#rightarrow Seed (all PU)");
	legend->AddEntry(pion->EffPt[4],"#pi#rightarrow Seed (<30)");
	legend->AddEntry(pion->EffPt[5],"#pi#rightarrow Seed (>30)");

        legend->Draw();

        can->Print(ProcessName+"_Pion_RecoEff_vs_Pt.pdf");

        pion->EffEta[1]->SetMaximum(1.2);
        pion->EffEta[1]->SetMinimum(0.001);
        pion->EffEta[1]->Draw("pe1");
        pion->EffEta[4]->Draw("pe1same");	
        pion->EffEta[5]->Draw("pe1same");
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
                float fr1=electron->mva_e_pi->Integral(k+1,Nbin2D)/electron->Pt[0]->Integral();
                float fr2=pion->mva_e_pi->Integral(k+1,Nbin2D)/pion->Pt[0]->Integral();
                x_ged_all[k]=fr1;
                y_ged_all[k]=fr2;
        }
        TGraph *gr_ged_all=new TGraph(Nbin2D,x_ged_all,y_ged_all);

        for(int k=0;k<Nbin2D;k++){
                float fr1=electron->mva_e_pi->Integral(k+1,Nbin2D)/electron->Pt[0]->Integral();
                float fr2=background->mva_e_pi->Integral(k+1,Nbin2D)/background->mva_e_pi->Integral();
                x_ged_all_effvsmistag[k]=fr1;
                y_ged_all_effvsmistag[k]=fr2;
        }
        TGraph *gr_ged_all_effvsmistag=new TGraph(Nbin2D,x_ged_all_effvsmistag,y_ged_all_effvsmistag);


        for(int k=0;k<Nbin2D;k++){
		float fr1=electron_lowPU->mva_e_pi->Integral(k+1,Nbin2D)/electron_lowPU->Pt[0]->Integral();
		float fr2=pion_lowPU->mva_e_pi->Integral(k+1,Nbin2D)/pion_lowPU->Pt[0]->Integral();
                x_ged_lpu[k]=fr1;
                y_ged_lpu[k]=fr2;
		cout<<"lowpu: bin"<<k<<" center of bin "<<electron_lowPU->mva_e_pi->GetBinCenter(k+1)<<" "<<electron_lowPU->mva_e_pi->Integral(k+1,Nbin2D)<<" "<<electron_lowPU->Pt[0]->Integral()<<" "<<fr1<<" "<<fr2<<endl;
        }
	TGraph *gr_ged_lpu=new TGraph(Nbin2D,x_ged_lpu,y_ged_lpu);
	for(int k=0;k<Nbin2D;k++){
		float fr1=electron_midPU->mva_e_pi->Integral(k+1,Nbin2D)/electron_midPU->Pt[0]->Integral();
                float fr2=pion_midPU->mva_e_pi->Integral(k+1,Nbin2D)/pion_midPU->Pt[0]->Integral();
                x_ged_mpu[k]=fr1;
                y_ged_mpu[k]=fr2;
		cout<<"midpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_ged_mpu=new TGraph(Nbin2D,x_ged_mpu,y_ged_mpu);
	for(int k=0;k<Nbin2D;k++){
		float fr1=electron_highPU->mva_e_pi->Integral(k+1,Nbin2D)/electron_highPU->Pt[0]->Integral();
                float fr2=pion_highPU->mva_e_pi->Integral(k+1,Nbin2D)/pion_highPU->Pt[0]->Integral();
                x_ged_hpu[k]=fr1;
                y_ged_hpu[k]=fr2;
		cout<<"highpu:"<<fr1<<" "<<fr2<<endl;
        }
	TGraph *gr_ged_hpu=new TGraph(Nbin2D,x_ged_hpu,y_ged_hpu);

/////////////////
        for(int k=0;k<Nbin2D;k++){
                float fr1=electron_lowPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron_lowPU->Pt[0]->Integral();
                float fr2=pion_lowPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion_lowPU->Pt[0]->Integral();
                x_pf_lpu[k]=fr1;
                y_pf_lpu[k]=fr2;
                cout<<"lowpu: bin"<<k<<" center of bin "<<electron_lowPU->mva_e_pi->GetBinCenter(k+1)<<" "<<electron_lowPU->mva_e_pi->Integral(k+1,Nbin2D)<<" "<<electron_lowPU->Pt[0]->Integral()<<" "<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_lpu=new TGraph(Nbin2D,x_pf_lpu,y_pf_lpu);
        for(int k=0;k<Nbin2D;k++){
                float fr1=electron_midPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron_midPU->Pt[0]->Integral();
                float fr2=pion_midPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion_midPU->Pt[0]->Integral();
                x_pf_mpu[k]=fr1;
                y_pf_mpu[k]=fr2;
                cout<<"midpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_mpu=new TGraph(Nbin2D,x_pf_mpu,y_pf_mpu);
        for(int k=0;k<Nbin2D;k++){
                float fr1=electron_highPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron_highPU->Pt[0]->Integral();
                float fr2=pion_highPU->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion_highPU->Pt[0]->Integral();
                x_pf_hpu[k]=fr1;
                y_pf_hpu[k]=fr2;
                cout<<"highpu:"<<fr1<<" "<<fr2<<endl;
        }
        TGraph *gr_pf_hpu=new TGraph(Nbin2D,x_pf_hpu,y_pf_hpu);

        for(int k=0;k<Nbin2D;k++){
                float fr1=electron->mva_e_pi_PF->Integral(k+1,Nbin2D)/electron->Pt[0]->Integral();
                float fr2=pion->mva_e_pi_PF->Integral(k+1,Nbin2D)/pion->Pt[0]->Integral();
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

        obj->Pt[0]->Fill(Evt->ptGen);
        obj->Eta[0]->Fill(Evt->etaGen);
        if(Evt->isMatchedWithASeed==1){
	      	obj->Pt[1]->Fill(Evt->ptGen);
                obj->Eta[1]->Fill(Evt->etaGen);
        }
	if(Evt->isEcalDrivenSeeded){
		obj->Pt[2]->Fill(Evt->ptGen);
                obj->Eta[2]->Fill(Evt->etaGen);
	}
	if(Evt->isEcalDrivenSeeded || Evt->isMatchedWithASeed){
                obj->Pt[3]->Fill(Evt->ptGen);
                obj->Eta[3]->Fill(Evt->etaGen);
        }
	if(Evt->isMatchedWithAGedGsfElec==1 && Evt->mva_e_pi>-0.1){
        	obj->Pt[4]->Fill(Evt->ptGen);
                obj->Eta[4]->Fill(Evt->etaGen);
		obj->mva_e_pi->Fill(Evt->mva_e_pi);
	}
	if(Evt->isMatchedWithAPFElec==1 && Evt->mva_e_pi_PF>-0.1){
                obj->Pt[5]->Fill(Evt->ptGen);
                obj->Eta[5]->Fill(Evt->etaGen);
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
		bool isfromV=(Evt->origin==1);
        	bool isfromX=(Evt->origin==0);
	//	if(Evt->ptGen>2 && fabs(Evt->etaGen)<2.4)cout<<Evt->ptGen<<" "<<Evt->etaGen<<" "<<isfromB<<" "<<abs(Evt->pdgId)<<endl;
		if(basicselection && isfromB && abs(Evt->pdgId)==11){
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
		if(basicselection && isfromV && abs(Evt->pdgId)==11){
                        TheFill(electron_isolated,Evt);
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
	electron_isolated->Combine();
	pion->Combine();
	pion_lowPU->Combine();
	pion_highPU->Combine();
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
		SeedEffPt[u]=plot[u]->electron->EffPt[1];
       		SeedEffEta[u]=plot[u]->electron->EffEta[1];
		GEDEffPt[u]=plot[u]->electron->EffPt[4];
                GEDEffEta[u]=plot[u]->electron->EffEta[4];	
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

