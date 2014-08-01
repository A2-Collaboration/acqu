//
// TheoryComp
//
// Compares asymmetries and differential cross sections to various theory
// predictions and also at some energies to TAPS data from Schmidt and Leukel.
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

//
// Asym
//
// This compares asymmetries
//
void Asym( UInt_t eg = 160, Bool_t save = kFALSE)
{
	Bool_t maid;
	UInt_t i, j, cthbins;
	Double_t sig, dsig, min, max, deg;
	Double_t theta[100], dtheta[100], as[100], das[100];
	Double_t interval, cth, th;
	TString name;

	maid = kFALSE;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Asym", 400, 10, 700, 500);

	min = -0.2;
	max = 0.4;
	if ( eg == 160) {
		cthbins = 7;
		name = "xs/pi0/asym_sergey/asym_sp_b2.dat";
		deg = 2.4;
	}
	else if ( eg == 170) {
		cthbins = 9;
		name = "xs/pi0/asym_sergey/asym_sp_b5.dat";
		deg = 1.2;
	}
	else if ( eg == 180) {
		cthbins = 9;
		name = "xs/pi0/asym_sergey/asym_sp_b9.dat";
		deg = 1.2;
	}
	else if ( eg == 280) {
		cthbins = 9;
		name = "xs/pi0/asym_sergey/asym_sp_b27.dat";
		min = -0.4;
		max = 0.6;
		deg = 1.2;
	}
	interval = 2/Double_t( cthbins);

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
			sscanf( name.Data(), "%d%lf%lf", &j, &sig, &dsig);
			cth = -1 + interval*(j - 0.5);
			theta[i] = acos( cth)/kD2R;
			dtheta[i] = 0;
			as[i] = sig;
			das[i] = dsig;
			i++;
		}
	}
	inFile.close();

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
	name = Form ( "E_{#gamma} = %d #pm %d MeV", eg, Int_t(deg));
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.5);
//	gr->SetMarkerColor( 2);
	gr->SetLineWidth( 4);
//	gr->SetLineColor( 2);
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
	gr->SetMinimum( min);
	gr->SetMaximum( max);
	gr->Draw( "AP");

	tex = new TLatex( 20, 0.20, "PRELIMINARY");              
	tex->SetTextColor( 17);
	tex->SetTextSize( 0.15);
//	tex->SetTextAngle( 15);
//	tex->SetLineWidth( 2);
	tex->Draw();
	gr->Draw( "P");
	gPad->RedrawAxis();

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	pt = new TLegend( 0.25, 0.15, 0.8, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);

	pt->AddEntry( gr, "CB-TAPS", "p");
	pt->Draw();

	// TAPS data.
	if ( ( eg == 160) || ( eg == 280)) {

		name = Form( "xs/pi0/previous/asym_%dMeV.dat", eg);
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

		gr2 = new TGraphErrors( i-1, theta, as, dtheta, das);
		gr2->SetMarkerStyle( 20);
		gr2->SetMarkerSize( 1.5);
		gr2->SetLineWidth( 4);
		gr2->SetMarkerColor( 4);
		gr2->SetLineColor( 4);
		gr2->Draw( "Psame");
		pt->AddEntry( gr2, "TAPS 2001", "p");
	}

	// Theory DMT2001 with standard parameters
	name = Form( "xs/pi0/theory/dmt2001_%dMeV.dat", eg);
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
	gr3->SetLineWidth( 5);
	gr3->SetLineColor( 2);
	gr3->Draw( "Lsame");
	pt->AddEntry( gr3, "DMT 2001", "l");

	// Theory MAID2007 with standard parameters
	if ( maid == kTRUE) {
		name = Form( "xs/pi0/theory/maid2007_%dMeV.dat", eg);
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
		gr5 = new TGraph( i-1, theta, as);
		gr5->SetLineWidth( 5);
//		gr5->SetLineStyle( 3);
		gr5->SetLineColor( 5);
		gr5->Draw( "Lsame");
		pt->AddEntry( gr5, "MAID 2007", "l");
	}

	// Theory ChPT
	if ( eg != 280) {
		name = Form( "xs/pi0/theory/chpt_%dMeV.dat", eg);
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
		gr4->SetLineWidth( 5);
		gr4->SetLineStyle( 2);
		gr4->SetLineColor( 6);
		gr4->Draw( "Lsame");
		pt->AddEntry( gr4, "ChPT with LEC b_{p} fit to TAPS data", "l");
	}

	name = Form( "plots/Pi0/AsymComp_%dMeV", eg);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void Asym90deg( Bool_t save = kFALSE)
{
	UInt_t i;
	Double_t en, sig, dsig;
	Double_t eg[100], deg[100], as[100], das[100];
	TString name;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Asym90deg", 400, 10, 700, 500);

	// Our new data
	name = "xs/pi0/asym_sergey/asym_sp_90deg.dat";
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
			sscanf( name.Data(), "%lf%lf%lf", &en, &sig, &dsig);
			eg[i] = en;
			deg[i] = 0;
			as[i] = sig;
			das[i] = dsig;
			i++;
		}
	}
	inFile.close();

	gr = new TGraphErrors( i-1, eg, as, deg, das);
	name = "#theta_{#pi} = 90 deg";
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.5);
//	gr->SetMarkerColor( 4);
	gr->SetLineWidth( 4);
//	gr->SetLineColor( 4);
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
	gr->SetMaximum( 0.4);
	gr->SetMinimum( -0.2);
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
			sscanf( name.Data(), "%*s%lf%*s%*s%*s%lf", &en, &sig);
			eg[i] = en;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr2 = new TGraph( i-1, eg, as);
	gr2->SetLineWidth( 5);
	gr2->SetLineColor( 2);
	gr2->Draw( "Lsame");

	// Axel's point at 160 MeV
	Double_t as_eg[1], as_as[1], as_deg[1], as_das[1];
	as_eg[0] = 159.5;
	as_deg[0] = 0;
	as_as[0] = 0.217;
	as_das[0] = 0.046;
	gr3 = new TGraphErrors( 1, as_eg, as_as, as_deg, as_das);
	gr3->SetMarkerStyle( 20);
	gr3->SetMarkerSize( 1.5);
	gr3->SetMarkerColor( 4);
	gr3->SetLineWidth( 4);
	gr3->SetLineColor( 4);
	gr3->Draw( "Psame");

	pt = new TLegend( 0.40, 0.15, 0.60, 0.30);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "CB-TAPS", "p");
	pt->AddEntry( gr2, "DMT 2001", "l");
	pt->AddEntry( gr3, "TAPS 2001", "p");
	pt->Draw();

	name = "plots/Pi0/Asym90deg";
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void AsymDMT( Bool_t save = kFALSE)
{
	UInt_t i, j, cthbins;
	Double_t theta[100], dtheta[100], as[100], das[100], eg[100], deg[100];
	Double_t th, sig, dsig, min, max, interval, en, cth;
	TString name;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymDMT", 400, 10, 1000, 500);
	c1->Divide( 2, 1);

	min = -0.2;
	max = 0.3;
	cthbins = 7;
	interval = 2/Double_t( cthbins);

	// Our data at 160 MeV
	c1->cd( 1);
	name = "xs/pi0/asym_sergey/asym_sp_b2.dat";
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
			sscanf( name.Data(), "%d%lf%lf", &j, &sig, &dsig);
			cth = -1 + interval*(j - 0.5);
			theta[i] = acos( cth)/kD2R;
			dtheta[i] = 0;
			as[i] = sig;
			das[i] = dsig;
			i++;
		}
	}
	inFile.close();

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
	name = "E_{#gamma} = 160 #pm 2 MeV";
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.5);
	gr->SetLineWidth( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{#pi} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetLimits( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetMinimum( min);
	gr->SetMaximum( max);
	gr->Draw( "AP");

	tex = new TLatex( 15, 0.20, "PRELIMINARY");              
	tex->SetTextColor( 17);
	tex->SetTextSize( 0.11);
//	tex->SetTextAngle( 15);
//	tex->SetLineWidth( 2);
	tex->Draw();
	gr->Draw( "P");
	gPad->RedrawAxis();

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	pt = new TLegend( 0.15, 0.15, 0.7, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);

	pt->AddEntry( gr, "CB-TAPS", "p");
	pt->Draw();

	// TAPS data at 160 MeV.
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

	gr2 = new TGraphErrors( i-1, theta, as, dtheta, das);
	gr2->SetMarkerStyle( 20);
	gr2->SetMarkerSize( 1.5);
	gr2->SetLineWidth( 4);
	gr2->SetMarkerColor( 4);
	gr2->SetLineColor( 4);
	gr2->Draw( "Psame");
	pt->AddEntry( gr2, "TAPS 2001", "p");

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

	gr3 = new TGraph( i-1, theta, as);
	gr3->SetLineWidth( 4);
	gr3->SetLineColor( 2);
	gr3->Draw( "Lsame");
	pt->AddEntry( gr3, "DMT 2001", "l");

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
	pt->AddEntry( gr4, "ChPT with LEC b_{p} fit to TAPS data", "l");

	// At 90 deg as a function of Egamma.
	c1->cd( 2);

	// Our new data
	name = "xs/pi0/asym_sergey/asym_sp_90deg.dat";
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
			sscanf( name.Data(), "%lf%lf%lf", &en, &sig, &dsig);
			eg[i] = en;
			deg[i] = 0;
			as[i] = sig;
			das[i] = dsig;
			i++;
		}
	}
	inFile.close();

	gr5 = new TGraphErrors( i-1, eg, as, deg, das);
	name = "#theta_{#pi} = 90 #pm 10 deg";
	gr5->SetTitle( name);
	gr5->SetMarkerStyle( 21);
	gr5->SetMarkerSize( 1.5);
//	gr5->SetMarkerColor( 4);
	gr5->SetLineWidth( 4);
//	gr5->SetLineColor( 4);
	gr5->GetXaxis()->SetTitleOffset( 1.1);
	gr5->GetYaxis()->SetTitleOffset( 1.0);
	gr5->GetYaxis()->SetTitleSize( 0.05);
	gr5->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr5->GetYaxis()->SetTitle( "#Sigma");
	gr5->GetXaxis()->SetLabelSize( 0.03);
	gr5->GetYaxis()->SetLabelSize( 0.03);
	gr5->GetXaxis()->SetLimits( 140, 205);
	gr5->GetXaxis()->CenterTitle();
	gr5->GetYaxis()->CenterTitle();
	gr5->SetMinimum( min);
	gr5->SetMaximum( max);
	gr5->Draw( "AP");

	l1 = new TLine( 140, 0, 205, 0);
	l1->Draw();

	// Theory DMT2001 with standard parameters
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
			sscanf( name.Data(), "%*s%lf%*s%*s%*s%lf", &en, &sig);
			eg[i] = en;
			as[i] = sig;
			i++;
		}
	}
	inFile.close();

	gr6 = new TGraph( i-1, eg, as);
	gr6->SetLineWidth( 5);
	gr6->SetLineColor( 2);
	gr6->Draw( "Lsame");

	// Axel's point at 160 MeV
	Double_t as_eg[1], as_as[1], as_deg[1], as_das[1];
	as_eg[0] = 159.5;
	as_deg[0] = 0;
	as_as[0] = 0.217;
	as_das[0] = 0.046;
	gr7 = new TGraphErrors( 1, as_eg, as_as, as_deg, as_das);
	gr7->SetMarkerStyle( 20);
	gr7->SetMarkerSize( 1.5);
	gr7->SetMarkerColor( 4);
	gr7->SetLineWidth( 4);
	gr7->SetLineColor( 4);
	gr7->Draw( "Psame");

	name = "plots/Pi0/AsymWith";
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}
