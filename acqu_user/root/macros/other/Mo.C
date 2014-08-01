void Cuts()
{
	c1 = new TCanvas("c1", "Random Subtraction", 800, 10, 700, 500);
	c1->SetFillColor(18);
	gStyle->SetPalette(1);   
	//gStyle->SetCanvasColor(33);

//	TFile *inFile = new TFile( "histograms/realdata2.root", "READ");
	TFile *inFile = new TFile( "histograms/ARH.root.root", "READ");

//	TH3D *hP = (TH3D*)inFile->Get( "ThetaCMCut1P_v_MMissCut1P_v_TChanCut1P");
	TH3D *hP = (TH3D*)inFile->Get( "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP");
//	TH3D *hR = (TH3D*)inFile->Get( "ThetaCMCut1R_v_MMissCut1R_v_TChanCut1R");
	TH3D *hR = (TH3D*)inFile->Get( "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR");
	TH3D *hS = (TH3D*)hP->Clone();
	hS->Sumw2();
	hS->Add( hR, -0.6667);
	gStyle->SetOptStat(0);
//	hS->Draw("colz");

//	hS->GetXaxis()->SetRangeUser( 137, 154); //tagger channel
	hS->GetYaxis()->SetRangeUser( 0, 180); //angle
	TH2D *hS_xy = hS->Project3D( "zy");

	/*
	TString name = "Pi0 Scattering Angle vs. Missing Mass";
	hS_xy->SetTitle( name);

	hS_xy->GetXaxis()->SetTitleOffset( 1.2);
	hS_xy->GetXaxis()->SetTitle("M_{miss} (MeV)");
	hS_xy->GetXaxis()->SetLabelSize( 0.03);
	hS_xy->GetXaxis()->CenterTitle();

	hS_xy->GetYaxis()->SetTitleOffset( 1.2);
	hS_xy->GetYaxis()->SetTitle("#theta_{#gamma} (deg)");
	hS_xy->GetYaxis()->SetLabelSize( 0.03);
	hS_xy->GetYaxis()->CenterTitle();

	//hS_xy->Draw("colz");
	*/

	TH1D *hS_x = hS->Project3D( "z");
	hS_x->SetTitle( "Missing Mass");
	hS_x->GetXaxis()->SetTitleOffset( 1.2);
	hS_x->GetXaxis()->SetTitle("M_{miss} (MeV)");
	hS_x->GetXaxis()->SetLabelSize( 0.03);
	hS_x->GetXaxis()->CenterTitle();
	hS_x->Draw();
	hS_x->Rebin( 4);

	l = new TLine( 800, 0, 1100, 0);
	l->SetLineStyle( 2);
	l->Draw();

//	hS_x->Draw("colz");

//	c1->Print( "plots/eps/Mmiss_Theta_Cut.eps");
	c1->Print( "plots/Mmiss_Theta_Cut.pdf");

}
void TChan()
{
	TString name;

	c1 = new TCanvas( "c1", "Random Subtraction", 20, 0, 700, 500);
	c1->SetFillColor( 18);
	gStyle->SetOptStat( 0);
	gStyle->SetPalette( 1);

	name = "histograms/JanCarbon.root";
	TFile *inFile = new TFile( name, "READ");

	name = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP";
	TH3D *hP = (TH3D*)inFile->Get( name);
	name = "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR";
	TH3D *hR = (TH3D*)inFile->Get( name);
	TH3D *hS = (TH3D*)hP->Clone();
	hS->Sumw2();
	hS->Add( hR, -0.625);

	hS->GetXaxis()->SetRangeUser( 10, 229); //tagger channel
//	hS->GetXaxis()->SetRangeUser( 137, 154); //tagger channel
	hS->GetYaxis()->SetRangeUser( 0, 180); //angle

	TH2D *hS_proj = hS->Project3D( "x");
	hS_proj->Draw();
	hS_proj->Rebin( 4);

	Double_t left = hS_proj->GetXaxis()->GetXmin();
	Double_t right = hS_proj->GetXaxis()->GetXmax();
	cout << left << " " << right << endl;
	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);
	l->Draw();

//	c1->Print( "plots/Mmiss_Theta_Cut.pdf");

}
