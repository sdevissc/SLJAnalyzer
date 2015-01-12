{
	TFile *f=new TFile("TTbar_test.root")	;
	f->cd();
	float Xaxis[9]={2,5,10,15,25,50,75,100,200};
	float Yaxis[9]={-2.5,-1.7,-1.3,-0.5,0,0.5,1.3,1.7,2.5};
	TH2F *h1 =new TH2F("h1","",8,Xaxis,8,Yaxis);
	TH2F *h2 =new TH2F("h2","",8,Xaxis,8,Yaxis);
	TH2F *Division =new TH2F("div","",8,Xaxis,8,Yaxis);
	tree_purity->Draw("eta:pt>>h1","ptGen>2 && abs(etaGen)<2.5 && abs(pdgId)==11 && (origin==4 || origin==6)");
	tree_purity->Draw("eta:pt>>h2","ptGen>2 && abs(etaGen)<2.5 && abs(pdgId)==211 && (origin>0)");
	h1->Scale(1.0/h1->Integral());
        h2->Scale(1.0/h2->Integral());
	Division->Add(h1);
	Division->Divide(h2);
	TCanvas *can = new TCanvas("h", "bla",800,300);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.05);
        can->SetRightMargin(0.15);
        can->SetBottomMargin(0.15);
	can->Divide(3,0);
	can->cd(1);
	h1->Draw("colz");
	can->cd(2);
        h2->Draw("colz");
	can->cd(3);
        Division->Draw("colz");
	TFile *output=new TFile("weightFile.root","recreate");
	output->cd();
	cout<<"Writing in the output rootfile"<<endl;
	Division->Write();
	output->Close();
	
}
