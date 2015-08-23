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

const int Nbin2D=10000; 

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


class HistoMaker
{
  public:
	HistoMaker(const TString& ,const TString &);
	~HistoMaker();
	void Combine();
	
	TH1F* TemplatePt;
	TH1F* TemplateEta;
	TH1F *Pt[8],*EffPt[8],*EffEta[8],*Eta[8];
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
	for(int i=0;i<8;i++){
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

	for(int i=1;i<8;i++){
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
	string ProcessFile;
        TString ProcName;
	TString TheTag;
	HistoMaker *electron;
	HistoMaker *pion;
	HistoMaker *kaon;
	HistoMaker *electron_lowPU;
        HistoMaker *electron_midPU;
        HistoMaker *electron_highPU;
	HistoMaker *electron_isolated[9];
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
	for(int i=0;i<9;i++){
		electron_isolated[i]=new HistoMaker(Form("elecisol_hpu_%d",i),TheTag);
	}
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
        Processes.open(ProcessFile.c_str());
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
        float x_pf_all[6][Nbin2D],y_pf_all[Nbin2D];
	 TGraph *gr_pf_all[6];
        for(int k=0;k<Nbin2D;k++){
                y_pf_all[k]=background->themva->Integral(k+1,Nbin2D)/background->themva->Integral();
                for(int i=0;i<6;i++) {
                        x_pf_all[i][k]=electron_isolated[i+1]->themva->Integral(k+1,Nbin2D)/electron_isolated[i+1]->themva->Integral();
                }

        }
	for(int i=0;i<6;i++) {
		gr_pf_all[i]=new TGraph(Nbin2D,y_pf_all, x_pf_all[i]);
		gr_pf_all[i]->SetMarkerStyle(25+i);
		gr_pf_all[i]->SetMarkerColor(1+i);
		gr_pf_all[i]->SetMarkerSize(0.1);
		gr_pf_all[i]->SetMinimum(0.5);
		gr_pf_all[i]->GetXaxis()->SetLimits(0.00,1.005);
        	gr_pf_all[i]->GetYaxis()->SetRangeUser(0.5,1.1);
	}
	TLegend *legend=new TLegend(0.25,0.8,0.55,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
	legend->Draw();

        gr_pf_all[0]->Draw("ap");
	gr_pf_all[1]->Draw("psame");
	gr_pf_all[2]->Draw("psame");
	gr_pf_all[3]->Draw("psame");

	legend->Clear();

        legend->AddEntry(gr_pf_all[0],"10<Pt<15 EDS","p");
	legend->AddEntry(gr_pf_all[1],"15<Pt<20 EDS","p");
	legend->AddEntry(gr_pf_all[2],"20<Pt EDS","p");
	legend->AddEntry(gr_pf_all[3],"10<Pt<15 EDS || TDS","p");
        legend->AddEntry(gr_pf_all[4],"15<Pt<20 EDS || TDS","p");
        legend->AddEntry(gr_pf_all[5],"20<Pt EDS || TDS","p");
	legend->Draw();
	gPad->SetGridx(1);
	gPad->SetGridy(1);	
        can->Print(ProcessName+"effvseff.pdf");
}


void SLPlotter::TheFill(HistoMaker *obj,EventElecComm *Evt){

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
	if(Evt->isMatchedWithAPFElec==1 && Evt->themva>-0.1){
                obj->Pt[5]->Fill(Evt->ptGen);
                obj->Eta[5]->Fill(Evt->etaGen);
                obj->themva->Fill(Evt->themva);
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
			TheFill(electron_isolated[0],Evt);
                        if(Evt->ptGen<10 			&& Evt->ecalseed==1)TheFill(electron_isolated[1],Evt);
			if(Evt->ptGen<20 && Evt->ptGen>10 	&& Evt->ecalseed==1)TheFill(electron_isolated[2],Evt);
			if(Evt->ptGen>20        	        && Evt->ecalseed==1)TheFill(electron_isolated[3],Evt);
			if(Evt->ptGen<10 			&& (Evt->trkseed==1 || Evt->ecalseed==1))TheFill(electron_isolated[4],Evt);
                        if(Evt->ptGen<20 &&  Evt->ptGen>10 	&& (Evt->trkseed==1 || Evt->ecalseed==1))TheFill(electron_isolated[5],Evt);
			if(Evt->ptGen>20                        && (Evt->trkseed==1 || Evt->ecalseed==1))TheFill(electron_isolated[6],Evt);
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
	electron_isolated[0]->Combine();
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

