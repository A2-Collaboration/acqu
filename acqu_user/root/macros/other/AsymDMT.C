//
// AsymDMT
//
// Compares asymmetry at approx. 160 MeV to ChPT, DMT2001, DMT2001 with M1-
// reduced by 15%, and the Schmidt measurement.
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

Double_t egam[] = { 155.0, 159.8, 164.6, 168.2, 170.6, 173.0, 175.4, 177.8,
	180.2, 182.6, 185.0, 187.4, 189.8, 192.3, 194.7, 197.1, 199.5, 201.9 };

void AsymDMT( Bool_t save = kFALSE)
{
	UInt_t i;
	Double_t theta[100], dtheta[100], as[100], das[100], eg[100];
	Double_t egg[1], degg[1], ass[1], dass[1];
	Double_t th, sig;
	TString name;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymDMT", 400, 10, 1000, 500);
	c1->Divide( 2, 1);

	// TAPS data at 160 MeV.
	c1->cd( 1);
	name = "xs/pi0/previous/asym_160MeV.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		inFile >> theta[i] >> as[i] >> das[i];
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	// Getting 90 deg point for later...
	egg[0] = 159.5;
	degg[0] = 0;
	ass[0] = as[3];
	dass[0] = das[3];

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
	name = "E_{#gamma} = 160 #pm 10 MeV";
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerColor( 4);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{#pi} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetLimits( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetMaximum( 0.4);
	gr->SetMinimum( -0.2);
	gr->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	// Theory DMT2001 with standard parameters
	name = "xs/pi0/theory/dmt2001_160MeV.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%lf%*s%*s%*s%*s%lf", &th, &sig);
			theta[i] = th;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr2 = new TGraph( i-1, theta, as);
	gr2->SetLineWidth( 4);
	gr2->SetLineColor( 2);
	gr2->Draw( "Lsame");

	// Theory DMT2001 with M1- adjusted by 15% parameters
	name = "xs/pi0/theory/dmt2001_160MeV_adj.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%lf%*s%*s%*s%*s%lf", &th, &sig);
			theta[i] = th;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr3 = new TGraph( i-1, theta, as);
	gr3->SetLineWidth( 4);
	gr3->SetLineColor( 2);
	gr3->SetLineStyle( 3);
	gr3->Draw( "Lsame");

	// Theory ChPT
	name = "xs/pi0/theory/chpt_160MeV.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%*s%lf%*s%lf", &th, &sig);
			theta[i] = th;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr4 = new TGraph( i-1, theta, as);
	gr4->SetLineWidth( 4);
	gr4->SetLineStyle( 2);
	gr4->SetLineColor( 6);
	gr4->Draw( "Lsame");

	pt = new TLegend( 0.15, 0.15, 0.7, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "Schmidt TAPS", "p");
	pt->AddEntry( gr2, "DMT 2001", "l");
	pt->AddEntry( gr3, "DMT 2001, M_{1-} lowered by 15%", "l");
	pt->AddEntry( gr4, "ChPT with LEC b_{p} fit to TAPS data", "l");
	pt->Draw();

	// At 90 deg as a function of Egamma.
	c1->cd( 2);

	// Theory DMT2001
	name = "xs/pi0/theory/dmt2001_90deg.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%*s%lf%*s%*s%*s%lf", &th, &sig);
			eg[i] = th;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr5 = new TGraph( i-1, eg, as);
	name = "#theta_{#pi} = 90 deg";
	gr5->SetTitle( name);
	gr5->SetLineWidth( 4);
	gr5->SetLineColor( 2);
	gr5->GetXaxis()->SetTitleOffset( 1.1);
	gr5->GetYaxis()->SetTitleOffset( 1.0);
	gr5->GetYaxis()->SetTitleSize( 0.05);
	gr5->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr5->GetYaxis()->SetTitle( "#Sigma");
	gr5->GetXaxis()->SetLabelSize( 0.03);
	gr5->GetYaxis()->SetLabelSize( 0.03);
	gr5->GetXaxis()->SetRangeUser( 145, 200);
	gr5->GetXaxis()->CenterTitle();
	gr5->GetYaxis()->CenterTitle();
	gr5->SetMaximum( 0.4);
	gr5->SetMinimum( -0.2);
	gr5->Draw( "AL");

	l1 = new TLine( 145, 0, 200, 0);
	l1->Draw();

	// Theory DMT2001 with M1- lowered by 15%
	name = "xs/pi0/theory/dmt2001_90deg_adj.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%*s%lf%*s%*s%*s%lf", &th, &sig);
			eg[i] = th;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr6 = new TGraph( i-1, eg, as);
	gr6->SetLineWidth( 4);
	gr6->SetLineColor( 2);
	gr6->SetLineStyle( 3);
	gr6->Draw( "Lsame");

	gr7 = new TGraphErrors( 1, egg, ass, degg, dass);
	gr7->SetMarkerStyle( 21);
	gr7->SetMarkerSize( 1.2);
	gr7->SetMarkerColor( 4);
	gr7->SetLineWidth( 2);
	gr7->SetLineColor( 4);
	gr7->Draw( "Psame");

	name = "plots/Pi0/AsymDMT";
//	name.Append( ".pdf");
	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void AsymEbin( UInt_t ebin = 0, Bool_t save = kFALSE)
{
	UInt_t i;
	Double_t junk, theta[100], dtheta[100], as[100], das[100], eg[100];
	TString name;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymEbin", 400, 10, 700, 500);

	// Our new data
	name = Form( "xs/pi0/asym_new/asym_%5.1fMeV.dat", egam[ebin]);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		inFile >> junk >> theta[i] >> as[i] >> das[i];
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
	name = Form( "Photon Asymmetry for E_{#gamma} = %5.1f MeV", egam[ebin]);
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerColor( 4);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{#pi} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetLimits( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
//	gr->SetMaximum( 0.3);
//	gr->SetMinimum( -0.3);
	gr->Draw( "AP");

	Double_t ymax = gr->GetYaxis()->GetXmax();
	Double_t ymin = gr->GetYaxis()->GetXmin();

	Double_t pos = 0.6*ymax;

	tex = new TLatex( 20, pos, "PRELIMINARY");              
	tex->SetTextColor( 17);
	tex->SetTextSize( 0.15);
//	tex->SetTextAngle( 15);
//	tex->SetLineWidth( 2);
	tex->Draw();

	gr->Draw( "P");
	gPad->RedrawAxis();

	l1 = new TLine( 0, 0, 180, 0);
	if ( ymin <= 0) l1->Draw();

	// Theory DMT2001 with standard parameters
	name = Form( "xs/pi0/theory/asym_dmt_%5.1fMeV.dat", egam[ebin]);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		inFile >> theta[i] >> as[i];
		i++;
	}
	inFile.close();

	gr2 = new TGraph( i-1, theta, as);
	gr2->SetLineWidth( 4);
	gr2->SetLineColor( 2);
	gr2->Draw( "Lsame");

	pt = new TLegend( 0.35, 0.15, 0.55, 0.25);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "CB-TAPS", "p");
	pt->AddEntry( gr2, "DMT 2001", "l");
	pt->Draw();

	name = Form( "plots/Pi0/AsymEbin_%d", ebin);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void Asym90deg( Bool_t save = kFALSE)
{
	UInt_t i;
	Double_t junk, eg[100], deg[100], as[100], das[100];
	TString name;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Asym90deg", 400, 10, 700, 500);

	// Our new data
	name = "xs/pi0/asym_new/asym_90deg.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		inFile >> eg[i] >> as[i] >> das[i];
		deg[i] = 0;
		i++;
	}
	inFile.close();

	gr = new TGraphErrors( i-1, eg, as, deg, das);
	name = "Photon Asymmetry for #theta_{#pi} = 90 deg";
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerColor( 4);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetLimits( 140, 205);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetMaximum( 0.3);
	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	tex = new TLatex( 146, 0.10, "PRELIMINARY");              
	tex->SetTextColor( 17);
	tex->SetTextSize( 0.15);
//	tex->SetTextAngle( 15);
//	tex->SetLineWidth( 2);
	tex->Draw();

	gr->Draw( "P");
	gPad->RedrawAxis();

	l1 = new TLine( 140, 0, 205, 0);
	l1->Draw();

	// Theory DMT2001 with standard parameters
	name = "xs/pi0/theory/asym_dmt_90deg.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while ( !inFile.eof())
	{
		inFile >> eg[i] >> as[i];
		i++;
	}
	inFile.close();

	gr2 = new TGraph( i-1, eg, as);
	gr2->SetLineWidth( 4);
	gr2->SetLineColor( 2);
	gr2->Draw( "Lsame");

	pt = new TLegend( 0.15, 0.75, 0.35, 0.85);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "CB-TAPS", "p");
	pt->AddEntry( gr2, "DMT 2001", "l");
	pt->Draw();

	name = "plots/Pi0/Asym90deg";
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}
