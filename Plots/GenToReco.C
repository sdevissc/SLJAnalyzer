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

const int Nbin2D=200; 

template <class Type>
void SetStyle(Type *a,int color, int style, int fillstyle, const TString& xtitle,const TString& ytitle, float min, float max){
        a->SetLineColor(color);
        a->SetMarkerColor(color);
        a->SetMarkerStyle(style);
        a->SetFillColor(color);
        a->SetFillStyle(fillstyle);
	a->SetMarkerSize(0.8);
        a->GetXaxis()->SetTitle(xtitle);
        a->GetYaxis()->SetTitle(ytitle);
        if(min!=-777.0 && max!=-777.0 ){
                a->SetMinimum(min);
                a->SetMaximum(max);
        }
}

template <class Type>
void SetStyle(Type *a,int color, int style, const TString& xtitle,const TString& ytitle, float min, float max){
        a->SetLineColor(color);
        a->SetMarkerColor(color);
        a->SetMarkerStyle(style);
        a->SetMarkerSize(0.8);
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
	void EffSelDiscr();
	TH1F* TemplatePt;
	TH1F* TemplateEta;
	TH1F *Pt[20],*EffPt[20],*EffEta[20],*Eta[20];

	TH1F* mvaDiscr[3][4];
	float Eff[3][2][2][Nbin2D];

	TH1F *mva_e_pi,*themva;	

	TH1F* ForIntegral;
	TGraph *eEffvspiEff;
	float vecXmva[Nbin2D],vecYmva[Nbin2D];

};

HistoMaker::HistoMaker(const TString& inputname,const TString & tag)
{
	const int nbinpt=17,nbineta=13;
	TH1::SetDefaultSumw2(kTRUE);
	cout<<"define the binning"<<endl;
	float etax[14]={-2.5 , -2.1  , -1.7 , -1.3 , -0.9 , -0.5 , -0.1 , 0.1 , 0.5 , 0.9 , 1.3 , 1.7 , 2.1 , 2.5 };
        float ptx[nbinpt+1]= { 2.0 , 2.5, 3.0 , 3.5, 4.0 ,  4.5, 5.0 ,5.5, 6,6.5, 7,7.5,8,10,13,17,25.0 ,50.0};

	cout<<"define pt and eta basic histos"<<endl;
	ForIntegral=new TH1F(inputname+tag+"forIntegral","",Nbin2D,-100000,100000);	
	TemplatePt=new TH1F(inputname+tag+"Pt[0]","",nbinpt,ptx);
	TemplateEta=new TH1F(inputname+tag+"Eta[0]","",nbineta,etax);
	themva=new TH1F(inputname+tag+"mva","",Nbin2D,-1.1,1.1);	
	for(int i=0;i<20;i++){
		Pt[i]=(TH1F*)TemplatePt->Clone(Form("Pt_%d",i));
		Eta[i]=(TH1F*)TemplateEta->Clone(Form("Eta_%d",i));
		EffPt[i]=(TH1F*)TemplatePt->Clone(Form("EffPt_%d",i));
                EffEta[i]=(TH1F*)TemplateEta->Clone(Form("EffEta_%d",i));
	}
	for(int i=0;i<3;i++)for(int j=0;j<4;j++){
                mvaDiscr[i][j]=(TH1F*)themva->Clone(Form("mva_%d",i));
        }

}

HistoMaker::~HistoMaker()                 // destructor, just an example
{
}



void HistoMaker::EffSelDiscr(){
	for(int k=0;k<Nbin2D;k++)for(int kin=0;kin<3;kin++)for(int seed=0;seed<2;seed++)for(int type=0;type<2;type++){
		int denom=(type==0)?2:3;
		Eff[kin][seed][type][k]         = mvaDiscr[kin][seed]->Integral(k+1,Nbin2D)/mvaDiscr[kin][denom]->Integral(0,Nbin2D+1);
	}

}

void HistoMaker::Combine(){
	TH1::SetDefaultSumw2(kTRUE);

	for(int i=1;i<20;i++){
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
	void getPlot(int,const TString&); 
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
	string ProcessFile;
        TString ProcName;
	TString TheTag;
	HistoMaker *electron;
	HistoMaker *pion;
	HistoMaker *kaon;
	HistoMaker *electron_lowPU;
        HistoMaker *electron_midPU;
        HistoMaker *electron_highPU;
	HistoMaker *electron_isolated[5];
        HistoMaker *pion_lowPU;
        HistoMaker *pion_midPU;
        HistoMaker *pion_highPU;
	HistoMaker *kaon_lowPU;
        HistoMaker *kaon_midPU;
        HistoMaker *kaon_highPU;	
	HistoMaker *background[5];	
	

		
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
	for(int i=0;i<5;i++){
		electron_isolated[i]=new HistoMaker(Form("elecisol_%d",i),TheTag);
		background[i]=new HistoMaker(Form("bckg_%d",i),TheTag);
	}
        pion_lowPU=new HistoMaker("pion_lpu",TheTag);
	pion_midPU=new HistoMaker("pion_mpu",TheTag);
	pion_highPU=new HistoMaker("pion_hpu",TheTag);
	kaon_lowPU=new HistoMaker("kaon_lpu",TheTag);
        kaon_midPU=new HistoMaker("kaon_mpu",TheTag);
        kaon_highPU=new HistoMaker("kaon_hpu",TheTag);
	
}

void SLPlotter::initialize(){
        ifstream Processes;
        Processes.open(ProcessFile.c_str());
        string processline;
        Processes >>   ProcName;
        Evt = new EventElecComm(ProcName);
//        nentries =Evt->fChain->GetEntriesFast();
	nentries =1500000;
}


SLPlotter::~SLPlotter()                 // destructor, just an example
{
}

void SLPlotter::getPlot(int log,const TString& process){
	initialize();
	ProcessName = process;
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
        can->SetLeftMargin(0.10);
        can->SetTopMargin(0.07);
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.1);

	TLegend *legend=new TLegend(0.12,0.8,0.45,0.92);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.020);

	TPaveText *text=new TPaveText(0.18,0.94,0.3,0.99,"NDC");
	text->AddText("CMS simulation");
	text->SetTextSize(0.035);
	text->SetBorderSize(0);
        text->SetFillColor(0);
	text->Draw();
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
	gPad->SetLogx(log);

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
	text->Draw();
	can->Print(Form(ProcessName+"_Electron_and_pion_TDS_SeedEff_vs_Pt_log%d_%s.pdf",log,ProcessFile.c_str()));	
	legend->Clear();
//----------------------------------------------------------------


        legend->SetTextSize(0.020);
//        legend->AddEntry(electron->EffPt[1],"e_{B}#rightarrow TDS");
        legend->AddEntry(electron->EffPt[2],"e: ECAL-driven ");
        legend->AddEntry(electron->EffPt[3],"e: ECAL-driven or Tracker-driven");
  //      legend->AddEntry(pion->EffPt[1],"#pi#rightarrow TDS");
        legend->AddEntry(pion->EffPt[2],"#pi: ECAL-driven");
        legend->AddEntry(pion->EffPt[3],"#pi:ECAL-driven or Tracker-driven");


        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx(1);

        SetStyle(electron->EffPt[1],1,20,3003,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(electron->EffPt[2],2,26,3004,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(electron->EffPt[3],4,22,3005,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(pion->EffPt[1],1,23,3006,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(pion->EffPt[2],2,24,3002,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(pion->EffPt[3],4,20,3017,"p_{T} of the generated particle","Efficiency",0.0,1.05);

    //    electron->EffPt[1]->Draw("pe2");
        electron->EffPt[2]->Draw("pe2");
        electron->EffPt[3]->Draw("pe2same");
   //     pion->EffPt[1]->Draw("pe2same");
        pion->EffPt[2]->Draw("pe2same");
        pion->EffPt[3]->Draw("pe2same");


        legend->Draw();
	text->Draw();
        can->Print(Form(ProcessName+"_Electron_and_pion_TDSvsEDSvsTDSOREDS_seedeff_vs_Pt_log%d_%s.pdf",log,ProcessFile.c_str()));
        legend->Clear();


//------------------------------------------------------------------------


        legend->SetTextSize(0.020);
        legend->AddEntry(electron->EffPt[1],"e_{B}#rightarrow TDS");
        legend->AddEntry(electron->EffPt[6],"e_{B}#rightarrow EDS (dR)");
        legend->AddEntry(electron->EffPt[2],"e_{B}#rightarrow EDS");
        legend->AddEntry(electron->EffPt[7],"e_{B}#rightarrow EDS (dR)");
	legend->AddEntry(electron_isolated[0]->EffPt[1],"e_{V}#rightarrow TDS");
        legend->AddEntry(electron_isolated[0]->EffPt[6],"e_{V}#rightarrow TDS (dR)");


        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx(log);

        SetStyle(electron->EffPt[1],1,20,3003,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(electron->EffPt[6],1,21,3004,"p_{T} of the generated particle","Efficiency",0.0,1.05);
        SetStyle(electron->EffPt[2],4,22,3005,"p_{T} of the generated particle","Efficiency",0.0,1.05);
	SetStyle(electron->EffPt[7],4,23,3006,"p_{T} of the generated particle","Efficiency",0.0,1.05);
	SetStyle(electron_isolated[0]->EffPt[1],2,20,3003,"p_{T} of the generated particle","Efficiency",0.99,550);
	SetStyle(electron_isolated[0]->EffPt[6],2,20,3004,"p_{T} of the generated particle","Efficiency",0.99,550);

        electron->EffPt[1]->Draw("pe2");
        electron->EffPt[6]->Draw("pe2same");
        electron->EffPt[2]->Draw("pe2same");
	electron->EffPt[7]->Draw("pe2same");
	electron_isolated[0]->EffPt[1]->Draw("pe2same");
	electron_isolated[0]->EffPt[6]->Draw("pe2same");


        legend->Draw();

        can->Print(Form(ProcessName+"_Electron_and_pion_TDSvsEDSvsTDSOREDS_CheckDeltaRseedeff_vs_Pt_log%d_%s.pdf",log,ProcessFile.c_str()));
        legend->Clear();

	text->Draw();

//------------------------------------------------------------------------	


        legend->SetTextSize(0.020);
        //ilegend->AddEntry(electron_isolated->EffPt[1],"e_{V}#rightarrow TDS");
        //legend->AddEntry(electron_isolated->EffPt[2],"e_{V}#rightarrow EDS");
        //legend->AddEntry(electron_isolated->EffPt[3],"e_{V}#rightarrow  TDS #wedge EDS");

	
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx(1);

        SetStyle(electron_isolated[0]->EffPt[1],1,20,3003,"p_{T} of the generated particle","Efficiency",0.99,550);
        SetStyle(electron_isolated[0]->EffPt[2],2,21,3004,"p_{T} of the generated particle","Efficiency",0.99,550);
        SetStyle(electron_isolated[0]->EffPt[3],2,22,3005,"p_{T} of the generated particle","Efficiency(#%)",0.99,550);
	
//	electron_isolated->EffPt[3]->GetYaxis()->SetNdivisions(520);
//	electron_isolated->EffPt[3]->Scale(100);
//        electron_isolated->EffPt[1]->Draw("pe2");
//        electron_isolated->EffPt[2]->Draw("pe2same");
//        electron_isolated->EffPt[3]->Draw("pe2same");
	TH1F*h =(TH1F*)electron_isolated[0]->EffPt[3]->Clone("improvement");
	h->Add(electron_isolated[0]->EffPt[2],-1);
	h->Scale(100);
	SetStyle(h,1,23,3005,"p_{T} of the generated electron","Efficiency Gain (%)",0,55);
//	legend->AddEntry(h,"#varepsilon_{TDS #wedge EDS}/#varepsilon_{EDS}");
	h->Draw("p");


  //      legend->Draw();
	text->Draw();
        can->Print(Form(ProcessName+"_IsolatedElectron_TDSvsEDSvsTDSOREDS_seedeff_vs_Pt_log%d_%s.pdf",log,ProcessFile.c_str()));
        legend->Clear();

//----------------------------------------------------------

	SetStyle(electron->EffEta[1],1,22,3004,"#eta_{T} of the generated particle","",0.001,1.2);
        SetStyle(electron_lowPU->EffEta[1],2,23,3005,"#eta_{T} of the generated particle","",0.001,1.2);
        SetStyle(electron_highPU->EffEta[1],3,24,3006,"#eta_{T} of the generated particle","",0.001,1.2);

	electron->EffEta[1]->Draw("hpe1");
        electron->EffEta[4]->Draw("hpe1same");
        electron->EffEta[5]->Draw("hpe1same");

	legend->Draw();
	text->Draw();
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
	text->Draw();
        legend->Draw();

        can->Print(ProcessName+"_Pion_RecoEff_vs_Pt.pdf");

        pion->EffEta[1]->SetMaximum(1.2);
        pion->EffEta[1]->SetMinimum(0.001);
        pion->EffEta[1]->Draw("pe1");
        pion->EffEta[4]->Draw("pe1same");	
        pion->EffEta[5]->Draw("pe1same");
	legend->Draw();
	text->Draw();
	can->Print(ProcessName+"_Pion_RecoEff_vs_Eta.pdf");

	legend->Clear();
//--------------------------------------------GSF id Pt Eff
	 legend->SetTextSize(0.020);
        legend->AddEntry(electron_isolated[0]->EffPt[8],"e: ecal-driven");
        legend->AddEntry(background[0]->EffPt[8],"#pi: ecal driven");
	legend->AddEntry(electron_isolated[0]->EffPt[9],"e: ecal-driven or tracker-driven");
        legend->AddEntry(background[0]->EffPt[9],"#pi: ecal driven or tracker driven");


        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx(1);
	gPad->SetLogy(1);
        SetStyle(electron_isolated[0]->EffPt[8],kBlack,20,"p_{T} of the generated particle","Efficiency",0.002,log==0?1.2:1.2);
        SetStyle(background[0]->EffPt[8],1,21,"p_{T} of the generated particle","Efficiency",0.002,log==0?1.2:1.2);
	SetStyle(electron_isolated[0]->EffPt[9],kRed,24,"p_{T} of the generated particle","Efficiency",0.002,log==0?1.2:1.2);
        SetStyle(background[0]->EffPt[9],2,25,"p_{T} of the generated particle","Efficiency",0.002,log==0?1.2:1.20);
	SetStyle(electron_isolated[0]->EffPt[19],kRed,24,"p_{T} of the generated particle","Efficiency",0.002,log==0?1.2:1.2);
        SetStyle(background[0]->EffPt[19],2,25,"p_{T} of the generated particle","Efficiency",0.002,log==0?1.2:1.20);

        electron_isolated[0]->EffPt[8]->Draw("p");
        background[0]->EffPt[8]->Draw("psame");
	electron_isolated[0]->EffPt[9]->Draw("psame");
        background[0]->EffPt[9]->Draw("psame");


        legend->Draw();
        text->Draw();
        can->Print(Form(ProcessName+"_Electron_and_background_EDS_TDS_SeedEff_vs_Pt_log%d_%s.pdf",log,ProcessFile.c_str()));

	legend->Clear();


	TH1F *div[5];
	div[0]=(TH1F*)electron_isolated[0]->EffPt[9]->Clone();
	div[0]->Add(electron_isolated[0]->EffPt[8],-1);
	div[0]->Scale(100);
	
	div[1]=(TH1F*)electron_isolated[0]->EffPt[19]->Clone();
        div[1]->Add(electron_isolated[0]->EffPt[18],-1);
        div[1]->Scale(100);

	div[4]=(TH1F*)electron_isolated[0]->EffPt[17]->Clone();
        div[4]->Add(electron_isolated[0]->EffPt[16],-1);
        div[4]->Scale(100);
		
	gPad->SetLogy(0);
//	div[0]->Draw("p");
	div[1]->Draw("p");
	div[4]->Draw("psame");

	TLegend *legend2=new TLegend(0.32,0.8,0.65,0.92);
        legend2->SetBorderSize(0);
        legend2->SetFillColor(0);
        legend2->SetTextSize(0.035);

	 legend2->SetTextSize(0.020);
   //     legend2->AddEntry(div[0],"mva>0");
	legend2->AddEntry(div[1],"H#rightarrow(ZZ)#rightarrow 4e WP.");
	legend2->AddEntry(div[4],"mva>0.8 WP");

	legend2->Draw();
	
	TPaveText *text2=new TPaveText(0.18,0.94,0.3,0.99,"NDC");
        text2->AddText("CMS simulation");
        text2->SetTextSize(0.035);
        text2->SetBorderSize(0);
        text2->SetFillColor(0);
	text2->Draw();		
	SetStyle(div[0],1,25,"p_{T} of the generated particle","Efficiency Gain (%)",0.002,log==0?30:30);
	SetStyle(div[1],2,26,"p_{T} of the generated particle","Efficiency Gain (%)",0.002,log==0?30:30);
	SetStyle(div[4],4,29,"p_{T} of the generated particle","Efficiency Gain (%)",0.002,log==0?30:30);
	can->Print(Form(ProcessName+"_DIV_Electron_and_background_EDS_TDS_SeedEff_vs_Pt_log%d_%s.pdf",log,ProcessFile.c_str()));
	gPad->SetLogy(1);
//-----------------------------------------Gsf Id Eta Eff
//	legend->Clear();
	legend->SetTextSize(0.020);
 //       legend->AddEntry(electron_isolated[0]->EffEta[8],"e: ecal-driven");
 //       legend->AddEntry(background[0]->EffEta[8],"#pi: ecal driven");
 //       legend->AddEntry(electron_isolated[0]->EffEta[9],"e: ecal-driven or tracker-driven");
 //       legend->AddEntry(background[0]->EffEta[9],"#pi: ecal driven or tracker driven");
        
        
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogy(0);
	gPad->SetLogx(0);
        SetStyle(electron_isolated[0]->EffEta[8],kBlack,20,"#eta of the generated particle","Efficiency",0.002,log==0?1.2:2);
        SetStyle(background[0]->EffEta[8],1,21,3006,"#eta of the generated particle","Efficiency",0.002,log==0?1.2:2);
        SetStyle(electron_isolated[0]->EffEta[9],kRed,24,"#eta of the generated particle","Efficiency",0.002,log==0?1.2:2);
        SetStyle(background[0]->EffEta[9],2,25,"#eta of the generated particle","Efficiency",0.002,log==0?1.2:2);

        electron_isolated[0]->EffEta[8]->Draw("p");
        background[0]->EffEta[8]->Draw("psame");
        electron_isolated[0]->EffEta[9]->Draw("psame");
        background[0]->EffEta[9]->Draw("psame");


        legend->Draw();
        text->Draw();
        can->Print(Form(ProcessName+"_Electron_and_background_EDS_TDS_SeedEff_vs_Eta_log%d_%s.pdf",log,ProcessFile.c_str()));


	delete legend;
	delete can;
}

void SLPlotter::MakeEffvsEffPlot(TString &process){
	setTDRStyle();
	TCanvas *can = new TCanvas("h", "bla",600,600);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.07);
        can->SetRightMargin(0.07);
        can->SetBottomMargin(0.15);
	ProcessName=process;
	 TGraph *gr_pf_all[3][2][2][2][2];

	int col[24];
        col[0]=1;
        col[1]=2;
        col[2]=1;
        col[3]=2;
        col[4]=1;
        col[5]=2;
	col[6]=1;
        col[7]=2;
        col[8]=1;
        col[9]=2;
        col[10]=1;
        col[11]=2;
	col[12]=1;
        col[13]=2;
        col[14]=1;
        col[15]=2;
        col[16]=1;
        col[17]=2;
	col[18]=1;
        col[19]=2;
        col[20]=1;
        col[21]=2;
        col[22]=1;
        col[23]=2;

int mark[24];
        mark[0]=21;
        mark[1]=25;
        mark[2]=21;
        mark[3]=25;
        mark[4]=21;
        mark[5]=25;
        mark[6]=22;
        mark[7]=26;
        mark[8]=22;
        mark[9]=26;
        mark[10]=22;
        mark[11]=26;
        mark[12]=23;
        mark[13]=32;
        mark[14]=23;
        mark[15]=32;
        mark[16]=23;
        mark[17]=32;
        mark[18]=34;
        mark[19]=28;
        mark[20]=34;
        mark[21]=28;
        mark[22]=34;
        mark[23]=28;
	
int mark2[12];
	mark2[0]=20;
	mark2[1]=24;
	mark2[2]=20;
	mark2[3]=24;
	mark2[4]=20;
	mark2[5]=24;
	mark2[6]=20;
        mark2[7]=24;
        mark2[8]=20;
        mark2[9]=24;
        mark2[10]=20;
        mark2[11]=24;
	TPaveText *text=new TPaveText(0.18,0.94,0.3,0.99,"NDC");
        text->AddText("CMS simulation");
        text->SetTextSize(0.035);
        text->SetBorderSize(0);
        text->SetFillColor(0);
	TPaveText *textkin[3];
	for(int kin=0;kin<3;kin++){
		textkin[kin]=new TPaveText(0.3,0.88,0.8,0.88,"NDC");
		textkin[kin]->SetTextSize(0.03);
                textkin[kin]->SetBorderSize(0);
                textkin[kin]->SetFillColor(0);
		textkin[kin]->SetTextAlign(10);
	}
        textkin[0]->AddText("P_{T}< 10 GeV");
//	textkin[1]->AddText("10 <P_{T}< 20 GeV");
	textkin[1]->AddText("P_{T}> 10 GeV");
	textkin[2]->AddText("P_{T} > 20 GeV");
	cout<<"Now filling the graphs"<<endl;
	TPaveText *textDef[2][2];
	
	for(int sig=0;sig<2;sig++)for(int bck=0;bck<2;bck++){
		textDef[sig][bck]=new TPaveText(0.3,0.81,0.8,0.86,"NDC");
		textDef[sig][bck]->SetTextSize(0.03);
		textDef[sig][bck]->SetTextAlign(10);
        	textDef[sig][bck]->SetBorderSize(0);
        	textDef[sig][bck]->SetFillColor(0);
	}	
	textDef[0][0]->AddText("Signal:     electrons from Z decay");
	textDef[0][0]->AddText("Background: any meson/baryon");
	textDef[1][0]->AddText("Signal:     electrons from B hadron decay");
        textDef[1][0]->AddText("Background: any meson/baryon");
	textDef[0][1]->AddText("Signal:     electrons from Z decay");
        textDef[0][1]->AddText("Background: pions");
	textDef[1][1]->AddText("Signal:     electrons from B hadron decay");
        textDef[1][1]->AddText("Background: pions");
	for(int kin=0;kin<3;kin++){
		textkin[kin]->SetTextSize(0.025);
		textkin[kin]->SetBorderSize(0);
		textkin[kin]->SetFillColor(0);
		for(int seed=0;seed<2;seed++)for(int type=0;type<2;type++) {
			gr_pf_all[kin][seed][type][0][0]=new TGraph(Nbin2D,background[0]->Eff[kin][seed][type],electron_isolated[0]->Eff[kin][seed][type]);
			gr_pf_all[kin][seed][type][1][0]=new TGraph(Nbin2D,background[0]->Eff[kin][seed][type],electron->Eff[kin][seed][type]);
			gr_pf_all[kin][seed][type][0][1]=new TGraph(Nbin2D,pion->Eff[kin][seed][type],electron_isolated[0]->Eff[kin][seed][type]);
			gr_pf_all[kin][seed][type][1][1]=new TGraph(Nbin2D,pion->Eff[kin][seed][type],electron->Eff[kin][seed][type]);
			for(int sig=0;sig<2;sig++)for(int bck=0;bck<2;bck++){
				gr_pf_all[kin][seed][type][sig][bck]->SetMarkerStyle(20+4*seed);
				gr_pf_all[kin][seed][type][sig][bck]->SetMarkerColor(1+seed);
				gr_pf_all[kin][seed][type][sig][bck]->SetMinimum(0.5);
				gr_pf_all[kin][seed][type][sig][bck]->GetXaxis()->SetLimits(0.00,(type==0)?0.4:(kin==0)?0.008:0.02);
				gr_pf_all[kin][seed][type][sig][bck]->GetXaxis()->SetNdivisions(404,0);
				gr_pf_all[kin][seed][type][sig][bck]->GetYaxis()->SetRangeUser(0.0055,1.4);
				gr_pf_all[kin][seed][type][sig][bck]->GetXaxis()->SetTitle("Background efficiency");
				gr_pf_all[kin][seed][type][sig][bck]->GetYaxis()->SetTitle("Signal efficiency");
				gr_pf_all[kin][seed][type][sig][bck]->GetXaxis()->SetTitleSize(0.04);
				gr_pf_all[kin][seed][type][sig][bck]->GetYaxis()->SetTitleSize(0.04);
				gr_pf_all[kin][seed][type][sig][bck]->GetYaxis()->SetTitleOffset(1.5);
				gr_pf_all[kin][seed][type][sig][bck]->GetXaxis()->SetTitleOffset(1.5);
			}
		}
	}
	cout<<"Now drawing"<<endl;
	for(int kin=0;kin<3;kin++)for(int type=0;type<2;type++)for(int sig=0;sig<2;sig++)for(int bck=0;bck<2;bck++){
		cout<<kin<<" "<<type<<" "<<sig<<" "<<bck<<endl;
		gr_pf_all[kin][0][type][sig][bck]->Draw("ap");
		gr_pf_all[kin][1][type][sig][bck]->Draw("psame");
		TLegend *legend=new TLegend(0.2,0.7,0.92,0.80);
        	legend->SetBorderSize(0);
        	legend->SetFillColor(0);
        	legend->SetTextSize(0.03);
		TString efftype = (type==0)?"eff ID":"eff RECO+ID";
		legend->AddEntry(gr_pf_all[kin][1][type][sig][bck],      "ecal-or tracker-driven "+efftype,"p");
		legend->AddEntry(gr_pf_all[kin][0][type][sig][bck],      "ecal-driven "+efftype,"p");
		legend->Draw();
		textkin[kin]->Draw();
		text->Draw();
		textDef[sig][bck]->Draw();
		can->Print(Form(ProcessName+"effvseff_%d_%d_%d_%d.pdf",kin,type,sig,bck));	
		delete legend;
	}
	// Getting the ID and RECO+ID plots merged
        for(int kin=0;kin<3;kin++)for(int type=0;type<2;type++)for(int sig=0;sig<2;sig++)for(int bck=0;bck<2;bck++){
                cout<<kin<<" "<<type<<" "<<sig<<" "<<bck<<endl;
                gr_pf_all[kin][0][0][sig][bck]->Draw("ap");
                gr_pf_all[kin][1][0][sig][bck]->Draw("psame");
		gr_pf_all[kin][0][1][sig][bck]->Draw("psame");
                gr_pf_all[kin][1][1][sig][bck]->Draw("psame");
                TLegend *legend=new TLegend(0.2,0.7,0.92,0.8);
                legend->SetBorderSize(1);
                legend->SetFillColor(0);
                legend->SetTextSize(0.03);
                TString efftype = (type==0)?"eff ID":"eff RECO+ID";
                legend->AddEntry(gr_pf_all[kin][0][0][sig][bck], "eff_{ID}: ecal-driven "+efftype,"p");
                legend->AddEntry(gr_pf_all[kin][1][0][sig][bck], "eff_{ID}: ecal-or tracker-driven "+efftype,"p");
		legend->AddEntry(gr_pf_all[kin][0][1][sig][bck], "eff_{RECO+ID}: ecal-driven "+efftype,"p");
                legend->AddEntry(gr_pf_all[kin][1][1][sig][bck], "eff_{RECO+ID}: ecal-or tracker-driven "+efftype,"p");
		gr_pf_all[kin][0][1][sig][bck]->SetMarkerStyle(21);
		gr_pf_all[kin][1][1][sig][bck]->SetMarkerStyle(25);
		
                legend->Draw();
                textkin[kin]->Draw();
                text->Draw();
                textDef[sig][bck]->Draw();
                can->Print(Form(ProcessName+"merged_effvseff_%d_%d_%d.pdf",kin,sig,bck));
                delete legend;
        }
	
}


void SLPlotter::TheFill(HistoMaker *obj,EventElecComm *Evt){

// CAUTION: THIS SECTION IS NOT PROTECTED ANYMORE BY A MVA MINIMAL CONDITION!!!!

        int flagEta=fabs(Evt->etaGen)<1.4?0:1;
	obj->ForIntegral->Fill(Evt->ptGen);
        obj->Pt[0]->Fill(Evt->ptGen);
        obj->Eta[0]->Fill(Evt->etaGen);
	obj->themva->Fill(Evt->mva);
        if(Evt->isDeltaRMatchedWithASeed==1){
	      	obj->Pt[1]->Fill(Evt->ptGen);
                obj->Eta[1]->Fill(Evt->etaGen);
        }
	if(Evt->isDeltaREcalDrivenSeeded){
		obj->Pt[2]->Fill(Evt->ptGen);
                obj->Eta[2]->Fill(Evt->etaGen);
	}
	if(Evt->isMatchedWithASeed==1){
                obj->Pt[6]->Fill(Evt->ptGen);
                obj->Eta[6]->Fill(Evt->etaGen);
        }
        if(Evt->isEcalDrivenSeeded==1){
                obj->Pt[7]->Fill(Evt->ptGen);
                obj->Eta[7]->Fill(Evt->etaGen);
        }

	if(Evt->isDeltaREcalDrivenSeeded || Evt->isDeltaRMatchedWithASeed){
                obj->Pt[3]->Fill(Evt->ptGen);
                obj->Eta[3]->Fill(Evt->etaGen);
        }
	if(Evt->isMatchedWithAGedGsfElec==1 && Evt->mva_e_pi>-0.1){
        	obj->Pt[4]->Fill(Evt->ptGen);
                obj->Eta[4]->Fill(Evt->etaGen);
		obj->mva_e_pi->Fill(Evt->mva_e_pi);
	}
	if(Evt->isMatchedWithAPFElec==1 && Evt->mva>-0.1){
                obj->Pt[5]->Fill(Evt->ptGen);
                obj->Eta[5]->Fill(Evt->etaGen);
                obj->themva->Fill(Evt->mva);
        }		
	
	int kincond=-1;
	if(Evt->ptGen<10)kincond=0;
//        if(Evt->ptGen>10 && Evt->ptGen<20)kincond=1;
	if(Evt->ptGen>10)kincond=1;
//        if(Evt->ptGen>20)kincond=2;
	obj->mvaDiscr[kincond][3]->Fill(Evt->mva);	
	if(Evt->ecalseed> -10 || Evt->trkseed>-10){
		obj->mvaDiscr[kincond][2]->Fill(Evt->mva);
		if(Evt->ecalseed==1){
			obj->mvaDiscr[kincond][0]->Fill(Evt->mva);
			if(Evt->mva>0.0){
				obj->Pt[8]->Fill(Evt->ptGen);
				obj->Eta[8]->Fill(Evt->etaGen);
			}
			if(Evt->mva>0.2){
                                obj->Pt[10]->Fill(Evt->ptGen);
                                obj->Eta[10]->Fill(Evt->etaGen);
                        }
			if(Evt->mva>0.4){
                                obj->Pt[12]->Fill(Evt->ptGen);
                                obj->Eta[12]->Fill(Evt->etaGen);
                        }
			if(Evt->mva>0.6){
                                obj->Pt[14]->Fill(Evt->ptGen);
                                obj->Eta[14]->Fill(Evt->etaGen);
                        }
			if(Evt->mva>0.8){
                                obj->Pt[16]->Fill(Evt->ptGen);
                                obj->Eta[16]->Fill(Evt->etaGen);
                        }
			if(Evt->ptGen<10){
				float feta=fabs(Evt->etaGen);
				float bdt=Evt->mva;
				bool cond1=feta<0.8 && bdt>0.47;
				bool cond2=feta>0.8 && feta<1.479 && bdt>0.004;
				bool cond3=feta>1.479 && bdt>0.295;
				if(cond1 || cond2 || cond3){
					obj->Pt[18]->Fill(Evt->ptGen);
                                	obj->Eta[18]->Fill(Evt->etaGen);
                        	}	
			}else{
				float feta=fabs(Evt->etaGen);
                                float bdt=Evt->mva;
                                bool cond1=feta<0.8 && bdt>-0.34;
                                bool cond2=feta>0.8 && feta<1.479 && bdt>-0.65;
                                bool cond3=feta>1.479 && bdt>0.6;
                                if(cond1 || cond2 || cond3){
                                        obj->Pt[18]->Fill(Evt->ptGen);
                                        obj->Eta[18]->Fill(Evt->etaGen);
                                }
			}
		}
		if(Evt->ecalseed==1 || Evt->trkseed==1){
			obj->mvaDiscr[kincond][1]->Fill(Evt->mva);
			if(Evt->mva>0.0){
				obj->Pt[9]->Fill(Evt->ptGen);
                        	obj->Eta[9]->Fill(Evt->etaGen);
			}
			if(Evt->mva>0.2){
                                obj->Pt[11]->Fill(Evt->ptGen);
                                obj->Eta[11]->Fill(Evt->etaGen);
                        }
			if(Evt->mva>0.4){
                                obj->Pt[13]->Fill(Evt->ptGen);
                                obj->Eta[13]->Fill(Evt->etaGen);
                        }
			if(Evt->mva>0.6){
                                obj->Pt[15]->Fill(Evt->ptGen);
                                obj->Eta[15]->Fill(Evt->etaGen);
                        }
			if(Evt->mva>0.8){
                                obj->Pt[17]->Fill(Evt->ptGen);
                                obj->Eta[17]->Fill(Evt->etaGen);
                        }
			if(Evt->ptGen<10){
                                float feta=fabs(Evt->etaGen);
                                float bdt=Evt->mva;
                                bool cond1=feta<0.8 && bdt>0.47;
                                bool cond2=feta>0.8 && feta<1.479 && bdt>0.004;
                                bool cond3=feta>1.479 && bdt>0.295;
                                if(cond1 || cond2 || cond3){
                                        obj->Pt[19]->Fill(Evt->ptGen);
                                        obj->Eta[19]->Fill(Evt->etaGen);
                                }
                        }else{
                                float feta=fabs(Evt->etaGen);
                                float bdt=Evt->mva;
                                bool cond1=feta<0.8 && bdt>-0.34;
                                bool cond2=feta>0.8 && feta<1.479 && bdt>-0.65;
                                bool cond3=feta>1.479 && bdt>0.6;
                                if(cond1 || cond2 || cond3){
                                        obj->Pt[19]->Fill(Evt->ptGen);
                                        obj->Eta[19]->Fill(Evt->etaGen);
                                }
                        }

			
		}
	}
}


void SLPlotter::FillHistos(){
	for(int ientry=0; ientry<nentries; ientry++) {
		if(ientry % 100000==0)cout<<"Event number "<<ientry<<endl;
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
			TheFill(electron_isolated[0],Evt);
			if(Evt->id_veto==1)TheFill(electron_isolated[1],Evt);
			if(Evt->id_loose==1)TheFill(electron_isolated[2],Evt);
			if(Evt->id_medium==1)TheFill(electron_isolated[3],Evt);
                        if(Evt->id_tight==1)TheFill(electron_isolated[4],Evt);
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
			TheFill(background[0],Evt);
                        if(Evt->id_veto==1)TheFill(background[1],Evt);
                        if(Evt->id_loose==1)TheFill(background[2],Evt);
                        if(Evt->id_medium==1)TheFill(background[3],Evt);
                        if(Evt->id_tight==1)TheFill(background[4],Evt);
		}
        }
	electron->Combine();
	electron_lowPU->Combine();
	electron_highPU->Combine();
	electron_isolated[0]->Combine();
	pion->Combine();
	pion_lowPU->Combine();
	pion_highPU->Combine();
	kaon->Combine();
	background[0]->Combine();
	electron_isolated[0]->EffSelDiscr();
	background[0]->EffSelDiscr();
	electron->EffSelDiscr();
	pion->EffSelDiscr();

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
        gStyle->SetTitleXOffset(1.25);
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
	
	plot[0]	=	new SLPlotter("inputSeedPlot_30_80_only","test1");
	plot[1]	=	new SLPlotter("inputSeedPlot_80_170_only","test2");
	plot[2]	=	new SLPlotter("inputSeedPlot_170_250_only","test3");
	
	TH1F *electron_EDS[3];
	TH1F *electron_any[3];
	TH1F *pion_EDS[3];
        TH1F *pion_any[3];
	TLegend *legend=new TLegend(0.18,0.8,0.5,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.012);
	for(int i=0;i<3;i++){
		plot[i]->setTDRStyle();
		TString test="test";
		plot[i]->getPlot(0,test);
	}
	for(int u=0;u<3;u++){
		electron_EDS[u]=(TH1F*)plot[u]->electron->EffPt[2];
       		electron_any[u]=(TH1F*)plot[u]->electron->EffPt[3];
		pion_EDS[u]=(TH1F*)plot[u]->pion->EffPt[2];
                pion_any[u]=(TH1F*)plot[u]->pion->EffPt[3];
		SetStyle(electron_EDS[u],1+u,26,3003,"p_{T} of the generated particle","Efficiency",0.0,1.25);
		SetStyle(electron_any[u],1+u,22,3003,"p_{T} of the generated particle","Efficiency",0.0,1.05);
		SetStyle(pion_EDS[u],1+u,24,3003,"p_{T} of the generated particle","Efficiency",0.0,1.05);
		SetStyle(pion_any[u],1+u,20,3003,"p_{T} of the generated particle","Efficiency",0.0,1.05);
	}

        TCanvas *can2 = new TCanvas("h", "bla",600,600);
        can2->SetLeftMargin(0.17);
        can2->SetTopMargin(0.05);
        can2->SetRightMargin(0.05);
        can2->SetBottomMargin(0.15);

//	gPad->SetGridx();
        gPad->SetGridy();
	gPad->SetLogx(1);	
	cout<<"getting the plots"<<endl;
        legend->AddEntry(electron_EDS[0],"e: ECAL-driven, 30-80","p");
	legend->AddEntry(electron_EDS[1],"e: ECAL-driven, 80-170","p");
	legend->AddEntry(electron_EDS[2],"e: ECAL-driven,170-250","p");
	legend->AddEntry(electron_any[0],"e: any, 30-80","p");
        legend->AddEntry(electron_any[1],"e: any, 80-170","p");
        legend->AddEntry(electron_any[2],"e: any, 170-250","p");
	legend->AddEntry(pion_EDS[0],"#pi: ECAL-driven, 30-80","p");
        legend->AddEntry(pion_EDS[1],"#pi: ECAL-driven, 80-170","p");
        legend->AddEntry(pion_EDS[2],"#pi: ECAL-driven,170-250","p");
        legend->AddEntry(pion_any[0],"#pi: any, 30-80","p");
        legend->AddEntry(pion_any[1],"#pi: any, 80-170","p");
        legend->AddEntry(pion_any[2],"#pi: any, 170-250","p");

	electron_EDS[0]->Draw("p");
	for(int u=0;u<3;u++){
                electron_EDS[u]->Draw("psame");
                electron_any[u]->Draw("psame");
                pion_EDS[u]->Draw("psame");
                pion_any[u]->Draw("psame");
        }
	legend->Draw();
	can2->Print("SeedG_Efficiencies_pt_comp.pdf");
}

