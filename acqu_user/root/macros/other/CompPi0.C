void Pi03D()
{
	c1 = new TCanvas("c1", "Random Subtraction", 800, 10, 700, 500);
	c1->SetFillColor(18);

	TH3D *hP = (TH3D*)inFile.Get( "MMissPi0_Theta_TChan_P");
	TH3D *hR = (TH3D*)inFile.Get( "MMissPi0_Theta_TChan_R");
	TH3D *hS = (TH3D*)hP->Clone();
	hS->Sumw2();
	hS->Add( hR, -0.5);

	gStyle->SetOptStat(0);

	hS->GetZaxis()->SetRangeUser( 900, 960);
//	hS->GetXaxis()->SetRangeUser( 200, 250);
//	hS->GetXaxis()->SetRangeUser( 301, 308);		// This is 145-163 MeV
	hS->GetXaxis()->SetRangeUser( 305, 308);		// This is 145-155 MeV
	TH2D *hS_yz = hS->Project3D( "yz");

	TString name = "Pi0 Scattering Angle vs. Missing Mass";
	hS_yz->SetTitle( name);

	hS_yz->GetXaxis()->SetTitleOffset( 1.2);
	hS_yz->GetXaxis()->SetTitle("M_{miss} (MeV)");
	hS_yz->GetXaxis()->SetLabelSize( 0.03);
	hS_yz->GetXaxis()->CenterTitle();

	hS_yz->GetYaxis()->SetTitleOffset( 1.2);
	hS_yz->GetYaxis()->SetTitle("#theta_{#gamma} (deg)");
	hS_yz->GetYaxis()->SetLabelSize( 0.03);
	hS_yz->GetYaxis()->CenterTitle();

	hS_yz->Draw("colz");

//	c1->Print( "plots/eps/Mmiss_Theta_Cut.eps");
}

void Pi0MmissCut()
{
	Double_t left, right;
	c1 = new TCanvas("c1", "Random Subtraction", 1000, 0, 700, 500);
	c1->SetFillColor(18);

	gStyle->SetOptStat(0);

	TH3D *hP = (TH3D*)inFile.Get( "MMissPi0_TGG_TChan_P");
	TH3D *hR = (TH3D*)inFile.Get( "MMissPi0_TGG_TChan_R");
	TH3D *hS = (TH3D*)hP->Clone();
	hS->Sumw2();
	hS->Add( hR, -0.5);

	left = 895;
	right = 950;

	hS->GetZaxis()->SetRangeUser( left, right);
	hS->GetYaxis()->SetRangeUser( 150, 180);
//	hS->GetXaxis()->SetRangeUser( 200, 250);
//	hS->GetXaxis()->SetRangeUser( 301, 308);		// This is 145-163 MeV
	hS->GetXaxis()->SetRangeUser( 306, 307);		// This is 148-153 MeV
	TH2D *hS_z = hS->Project3D( "z");

	TString name = "Cut Pi0 Missing Mass";
	hS_z->SetTitle( name);

	hS_z->GetXaxis()->SetTitleOffset( 1.2);
	hS_z->GetXaxis()->SetTitle("M_{miss} (MeV)");
	hS_z->GetXaxis()->SetLabelSize( 0.03);
	hS_z->GetXaxis()->CenterTitle();

	hS_z->Draw();
	l1 = new TLine( left, 0, right, 0);
	l1->SetLineStyle(1);
	l1->SetLineWidth(1);
	l1->Draw();

	Double_t max = hS_z->GetMaximum();
	max *= 1.4;
	hS_z->SetMaximum( max);

	Double_t min = hS_z->GetMinimum();
	min *= 1.4;
	hS_z->SetMinimum( min);

	l2 = new TLine( 934, min, 934, max);
	l2->SetLineStyle(2);
	l2->SetLineWidth(2);
	l2->Draw();

	c1->Print( "plots/pdf/MmissPi0_Cut.pdf");
}

void Pi0TGG2D()
{
	c1 = new TCanvas("c1", "Random Subtraction", 1000, 0, 700, 500);
	c1->SetFillColor(18);

	gStyle->SetOptStat(0);

	TH3D *hP = (TH3D*)inFile.Get( "MMissPi0_TGG_TChan_P");
	TH3D *hR = (TH3D*)inFile.Get( "MMissPi0_TGG_TChan_R");
	TH3D *hS = (TH3D*)hP->Clone();
	hS->Sumw2();
	hS->Add( hR, -0.5);

	hS->GetZaxis()->SetRangeUser( 900, 960);
	hS->GetYaxis()->SetRangeUser( 100, 180);
//	hS->GetXaxis()->SetRangeUser( 200, 250);
//	hS->GetXaxis()->SetRangeUser( 301, 308);		// This is 145-163 MeV
	hS->GetXaxis()->SetRangeUser( 305, 308);		// This is 145-155 MeV
	TH2D *hS_yz = hS->Project3D( "yz");

	TString name = "#gamma#gamma Opening Angle vs. Missing Mass";
	hS_yz->SetTitle( name);

	hS_yz->GetXaxis()->SetTitleOffset( 1.2);
	hS_yz->GetXaxis()->SetTitle("M_{miss} (MeV)");
	hS_yz->GetXaxis()->SetLabelSize( 0.03);
	hS_yz->GetXaxis()->CenterTitle();

	hS_yz->GetYaxis()->SetTitleOffset( 1.2);
	hS_yz->GetYaxis()->SetTitle("#theta_{#gamma#gamma} (deg)");
	hS_yz->GetYaxis()->SetLabelSize( 0.03);
	hS_yz->GetYaxis()->CenterTitle();

	hS_yz->Draw("colz");

//	c1->Print( "plots/eps/Mmiss_Theta_Cut.eps");
}

