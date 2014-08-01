// Macro Theory
//
//  Bin		Eg			Chan
// 	2		147		250
// 	3		153		245
// 	4		157		242
// 	5		163		237
// 	6		168		233
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void TheoryDiff( Int_t bin)
{
	Int_t i, ct, eg;
	Double_t eg_th;
	Double_t eng, dsig, asy, th, x, dx;
	Double_t theta[100], dtheta[100], xs[100], dxs[100], asym[100];
	Double_t max;
	Char_t lnstr[256];

	if ( bin == 2) {
		eg_th = 147.5;
		eg = 147;
	}
	else if ( bin == 3) {
		eg_th = 152.5;
		eg = 153;
	}
	else if ( bin == 4) {
		eg_th = 157.5;
		eg = 157;
	}
	else if ( bin == 5) {
		eg_th = 162.5;
		eg = 163;
	}
	else if ( bin == 6) {
		eg_th = 167.5;
		eg = 168;
	}
	else {
		cout << "<E> bin number ";
		cout << bin;
		cout << " is not valid.  Must be 2-6.";
		cout << endl;
		break;
	}

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// My results
	i = 0;
	TString filename( Form( "xs/diff_xs/diff_xs_%d.out", eg));
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	ifstream inFile( filename);
	max = 0;
	while( !inFile.eof()) {
		inFile >> th >> x >> dx;
		theta[i] = th;
		dtheta[i] = 0;
		xs[i] = x;
		if ( max < xs[i]) max = xs[i];
		dxs[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr = new TGraphErrors( ct, theta, xs, dtheta, dxs);
	TString title( Form( "Preliminary #gammap#rightarrowp#pi^{0} Differential"
				" Cross Section for E_{#gamma} = %d MeV", eg));
	gr->SetTitle( title);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 20);
	gr->SetLineWidth( 2);
	gr->GetXaxis()->SetTitleOffset( 1.15);
	gr->GetYaxis()->SetTitleOffset( 0.95);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("#theta* (deg)");
	gr->GetYaxis()->SetTitle("d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	if ( bin == 2) max *= 1.8;
	else max *= 1.4;
	gr->GetYaxis()->SetRangeUser(0,max);
	gr->GetXaxis()->SetRangeUser(0,180);
	gr->Draw( "AP");

	// Axel's results
	if ( bin == 6) filename = "xs/previous/dxs_167.dat";
	else filename = TString( Form( "xs/previous/dxs_%d.dat", eg));
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> th >> x >> dx;
		theta[i] = th;
		dtheta[i] = 0;
		xs[i] = x;
		dxs[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	// Plot Axel's Results
	gr1 = new TGraphErrors( ct, theta, xs, dtheta, dxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 21);
	gr1->Draw( "Psame");

	// Theory
	filename = "xs/theory/hbchpt_1996.dat";
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile.getline( lnstr, 256);
		if ( lnstr[0] != '#') {
			sscanf( lnstr, "%lf%lf%lf%lf", &eng, &th, &dsig, &asy);
			if ( eng == eg_th) {
				theta[i] = th;
				xs[i] = dsig*19.99/1000;
				asym[i++] = asy;
			}
		}
	}
	inFile.close();
	ct = i;

	// Plot Theory Results
	gr2 = new TGraph( ct, theta, xs);
	gr2->SetLineWidth( 4);
	gr2->SetLineColor( 4);
	gr2->Draw( "Csame");

	// Legend
	Double_t x1, x2, y1, y2;
	if ( bin == 2) {
		x1 = 0.25;
		x2 = 0.45;
		y1 = 0.7;
		y2 = 0.85;
	}
	else {
		x1 = 0.6;
		x2 = 0.8;
		y1 = 0.25;
		y2 = 0.40;
	}
	pt = new TLegend( x1, y1, x2, y2);
	pt->SetTextSize(0.04);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr1, "Schmidt - TAPS", "p");
	pt->AddEntry( gr2, "ChPT - Meissner", "l");
	pt->Draw();

	TString psfile( Form( "plots/eps/dxs_comp_%d.eps", eg));
	c1->Print( psfile);
}

void TheoryAsym( Int_t bin)
{
	Int_t i, ct, eg;
	Double_t eg_th;
	Double_t eng, th, x, dx, jnk;
	Double_t theta[100], dtheta[100], xs[100], dxs[100], asym[100], dasym[100];
	Double_t min, max;
	Char_t lnstr[256];

	if ( bin == 2) {
		eg_th = 147.5;
		eg = 147;
	}
	else if ( bin == 3) {
		eg_th = 152.5;
		eg = 153;
	}
	else if ( bin == 4) {
		eg_th = 157.5;
		eg = 158;
	}
	else if ( bin == 5) {
		eg_th = 162.5;
		eg = 163;
	}
	else if ( bin == 6) {
		eg_th = 167.5;
		eg = 168;
	}

	c1 = new TCanvas( "c1", "Asymmetries", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// My results
	i = 0;
	TString filename( Form( "xs/asym/asym_corr_%d.out", eg));
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	ifstream inFile( filename);
	max = 0;
	min = 1;
	while( !inFile.eof()) {
		inFile >> th >> x >> dx >> jnk;
		theta[i] = th;
		dtheta[i] = 0;
		asym[i] = x;
		if ( max < asym[i]) max = asym[i];
		if ( min > asym[i]) min = asym[i];
		dasym[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	if ( bin == 2) max *= 1.4;
	else if ( bin == 3) max *= 1.5;
	else if ( bin == 4) max *= 2.0;
	else if ( bin == 5) max *= 2.5;
	else if ( bin == 6) max *= 3.0;

	min *= 1.4;

	// Plot Results
	gr = new TGraphErrors( ct, theta, asym, dtheta, dasym);
	TString title( Form( "Preliminary #gammap#rightarrowp#pi^{0} Asymmetries"
				" for E_{#gamma} = %d MeV", eg));
	gr->SetTitle( title);
	gr->SetMarkerColor( 4);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.3);
	gr->GetYaxis()->SetTitleOffset( 0.7);
	gr->GetYaxis()->SetTitleSize( 0.06);
	gr->GetXaxis()->SetTitle("#theta* (deg)");
	gr->GetYaxis()->SetTitle("#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetRangeUser(min,max);
	gr->GetXaxis()->SetRangeUser(0,180);
	gr->Draw( "AP");

	// Theory
	filename = "xs/theory/hbchpt_1996.dat";
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile.getline( lnstr, 256);
		if ( lnstr[0] != '#') {
			sscanf( lnstr, "%lf%lf%lf%lf", &eng, &th, &x, &dx);
			if ( eng == eg_th) {
				theta[i] = th;
				xs[i] = x*19.99/1000;
				asym[i++] = dx;
			}
		}
	}
	inFile.close();
	ct = i;

	// Plot Results
	gr1 = new TGraph( ct, theta, asym);
	gr1->SetLineWidth( 4);
	gr1->SetLineColor( 2);
	gr1->Draw( "Csame");

	// Legend
	pt = new TLegend(0.6,0.70,0.8,0.85);
	pt->SetTextSize(0.04);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr1, "ChPT - Meissner", "l");
	pt->Draw();

	TString psfile( Form( "plots/eps/asym_comp_%d.eps", eg));
	c1->Print( psfile);
}
