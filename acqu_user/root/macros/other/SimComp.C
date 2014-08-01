// Macro SimComp.C
//
// Compares simulated EmPi0 and Tgg with actual data.  Also needed to calculate
// the carbon contamination for the asymmetry extraction.
// The contamination runs are:
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

TFile data( "histograms/ARH_full.root");

// For EmPi0
TH2D *d2p = (TH2D*)data.Get( "THits_EmPi0_P");
TH2D *d2r = (TH2D*)data.Get( "THits_EmPi0_R");

// For Tgg
TH3D *h3p = (TH3D*)data.Get( "THits_TGG_ThetaCM_P");
TH3D *h3r = (TH3D*)data.Get( "THits_TGG_ThetaCM_R");

void SimComp( Int_t chan)
{
	Int_t i, bin, eg;
	Int_t ch1, ch2, ch3, ch4;
	Int_t cdat_c, cdat_h, csim_c, csim_h;
	Double_t egamma[284];
	Double_t fact, eng, junk;
	Double_t x1, x2, x3, a, b, max, min;
	Double_t pa, pa_ratio[284];
	Double_t y1, y2;
	Double_t c, d, e, f, g, h;
	Double_t xx[1], yy[1], dxx[1], dyy[1];
	Char_t name[64], cond[] = "masscorr";

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 500, 0, 700, 500);
//   c1->SetFillColor( 38);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	strcpy( name, "xs/chan_subt/chan_subt_full.out");
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> pa ;
		pa_ratio[i] = pa;
	}
	inFile.close();

	strcpy( name, "xs/tageng450.out");
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> junk >> junk >> junk >> eng ;
		egamma[i] = eng;
	}
	inFile.close();

	eg = (int)(egamma[chan] + 0.5);
	pa = pa_ratio[chan];

// Data
	// Prompt EmPi0
	bin = d2p->GetYaxis()->FindBin( chan);
	TH1D *d1p = d2p->ProjectionX( "Data", bin, bin);

	// Random EmPi0
	bin = d2r->GetYaxis()->FindBin( chan);
	TH1D *d1r = d2r->ProjectionX( "Random", bin, bin);

/*	x1 = Linear( 145, 936, 400, 880, eg);
	x2 = Linear( 145, 944, 400, 1000, eg);
	x3 = Linear( 145, 910, 400, 880, eg);
	a = Linear( 145, 900, 400, 600, eg);
	b = Linear( 145, 950, 400, 1000, eg);
*/

	x1 = Linear( 145, 933, 400, 880, eg);
	x2 = Linear( 145, 945, 400, 1000, eg);
	x3 = Linear( 145, 918, 400, 880, eg);
	a = Linear( 145, 900, 400, 600, eg);
	b = Linear( 145, 950, 400, 1000, eg);

	// Subtracted EmPi0
	d1r->Sumw2();
	d1r->Scale( pa);
	d1p->Add( d1r, -1.0);
	d1p->SetLineStyle(1);
	d1p->SetMarkerStyle(20);
	d1p->SetMarkerSize(1.0);
	d1p->SetLineWidth(3);
	d1p->GetXaxis()->SetRangeUser(a,b);
	sprintf( name, "M_{miss} for E_{#gamma} = %5.1f MeV", egamma[chan]);
	d1p->SetTitle( name);
	d1p->GetXaxis()->SetTitle("M_{miss} (MeV)");
	d1p->GetXaxis()->SetTitleOffset( 1.1);
	d1p->GetXaxis()->CenterTitle();
	d1p->GetYaxis()->SetTitle("Counts");
	d1p->GetYaxis()->SetTitleOffset( 1.2);
	d1p->GetYaxis()->CenterTitle();

	ch2 = d1p->GetXaxis()->FindBin( x2);
	ch3 = d1p->GetXaxis()->FindBin( x3);
	ch4 = d1p->GetXaxis()->FindBin( x1);
	Double_t val, mval = 0;
	for ( Int_t j = ch4; j < ch2; j++) {
		val = d1p->GetBinContent(j);
		if ( val > mval) {
			mval = val;
			ch1 = j;
		}
	}
	ch1--;

	cdat_h = d1p->Integral( ch1, ch2);
	cdat_c = d1p->Integral( ch3, ch4);

	cout << " cdat_c( ";
	cout << d1p->GetBinCenter( ch3);
	cout << ", ";
	cout << d1p->GetBinCenter( ch4);
	cout << ") = ";
	cout << cdat_c << endl;

	cout << " cdat_h( ";
	cout << d1p->GetBinCenter( ch1);
	cout << ", ";
	cout << d1p->GetBinCenter( ch2);
	cout << ") = ";
	cout << cdat_h << endl;

	cout << "ratio = " << (double) cdat_h/cdat_c << endl;

	max = d1p->GetMaximum();
	max *= 1.4;
	d1p->SetMaximum(max);
	min = d1p->GetMinimum();
	min *= 3;
	d1p->SetMinimum(min);
	gStyle->SetOptStat(0);

	l3 = new TLine( a, 0, b, 0);
//	l3->SetLineColor(1);
//	l3->SetLineStyle(1);
	l3->SetLineWidth(2);

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineColor(2);
	l1->SetLineStyle(2);
	l1->SetLineWidth(4);

	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineColor(2);
	l2->SetLineStyle(2);
	l2->SetLineWidth(4);

	l4 = new TLine( x3, 0, x3, max);
	l4->SetLineColor(2);
	l4->SetLineStyle(2);
	l4->SetLineWidth(4);

// Simulation
	// Hydrogen
	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, eg);
	sim_p = new TFile( name);
//	s2h = (TH2D*) sim_p->Get( "THits_EmPi0_P");
//	TH1D *s1h = s2h->ProjectionX( "Simulation");
	s3h = (TH3D*) sim_p->Get( "EmPi0_TGG_ThetaCM");
	TH1D *s1h = s3h->Project3D("z");

	// Carbon
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	sim_c = new TFile( name);
//	s2c = (TH2D*) sim_c->Get( "THits_EmPi0_P");
//	TH1D *s1c = s2c->ProjectionX( "carb");
	s3c = (TH3D*) sim_c->Get( "EmPi0_TGG_ThetaCM");
	TH1D *s1c = s3c->Project3D("z2");

	csim_h = s1h->Integral( ch1, ch2);
	csim_c = s1c->Integral( ch3, ch4);

	s1h->SetFillStyle(3004);
	s1h->SetFillColor(4);
	s1h->SetLineColor(4);
	s1h->SetLineWidth(4);
	fact = (double) cdat_h/csim_h;
	s1h->Scale( fact);

	Double_t tot_h = s1h->Integral();

	s1c->SetFillStyle(3005);
	s1c->SetFillColor(2);
	s1c->SetLineColor(2);
	s1c->SetLineWidth(4);
	fact = (double) cdat_c/csim_c;
	s1c->Scale( fact);

	Double_t tot_c = s1c->Integral();

	cout << tot_h << " " << tot_c << " " << tot_h/tot_c << endl;

	TH1D *s1 = (TH1D*)s1h->Clone( "totsim");
	s1->Add(s1c, 1.0);

	y1 = s1c->Integral( ch4, ch2);
	y2 = s1->Integral( ch4, ch2);
	cout << " Contamination " << endl;
	cout << " y1 = " << y1 << " y2 = " << y2 << endl;
	cout << " a_p = " << (y2-y1)/y2 << endl;
	cout << " a_c = " << y1/y2 << endl;

	s1->SetFillStyle(0);
	s1->SetLineColor(3);
	s1->SetLineWidth(4);

	d1p->Draw("E1P");
	s1->Draw("same");
	s1c->Draw("same");
	s1h->Draw("same");
	l1->Draw();
	l2->Draw();
	l3->Draw();
	l4->Draw();

	Char_t contstr[32];
	sprintf( contstr, "a_{c} = %4.1f %%\n", y1/y2*100);
	pl = new TPaveLabel(0.2,0.4,0.45,0.5, contstr, "NDC");
	pl->SetBorderSize(0);
	pl->SetFillColor(0);
	pl->SetTextSize(0.45);
	pl->Draw();

	pt = new TLegend(0.2,0.6,0.4,0.85);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->AddEntry( d1p,"Data","p");
	pt->AddEntry( s1c,"Carbon","l");
	pt->AddEntry( s1h,"Hydrogen","l");
	pt->AddEntry( s1,"Sum","l");
	pt->Draw();

	sprintf( name, "plots/eps/sc_em_%d.eps", eg);
	c1->Print( name);
}

void SimCompTgg( Int_t chan)
{
	Int_t i, bin, egamma;
	Int_t ch1, ch2, ch3, ch4;
	Int_t cdat_c, cdat_h, csim_c, csim_h;
	Double_t Tgg, Tgg_c12;
	Double_t fact, junk, energy[284], eg;
	Double_t x1, x2, x3, a, b, max, min;
	Double_t pa, pa_ratio[284];
	Double_t y1, y2;
	Double_t c, d, e, f, g, h;
	Double_t t1, t2, t3, t4, t5;
	Double_t xx[1], yy[1], dxx[1], dyy[1];
	Char_t name[64], cond[] = "masscorr";

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetFillColor( 0);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	strcpy( name, "xs/chan_subt/chan_subt.out");
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> pa ;
		pa_ratio[i] = pa;
	}
	inFile.close();

	strcpy( name, "xs/tageng450.out");
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> junk >> junk >> junk >> eg ;
		energy[i] = eg;
	}
	inFile.close();

	eg = energy[chan];
	egamma = (int)(eg + 0.5);
	pa = pa_ratio[chan];

	if ( eg >= E_thr( kMP_MEV, kMPI0_MEV))
		Tgg = Tgg_Min( qT_max(eg, kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg = 180.;

	if ( eg >= E_thr( kM_C12_MEV, kMPI0_MEV))
		Tgg_c12 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg_c12 = 180.;

	// Prompt
	h3p->GetXaxis()->SetRangeUser(0,180);		// The entire ThetaCM range
	h3p->GetZaxis()->SetRangeUser(chan,chan);	// The appropriate tagger channel
	TH1D *h1p = h3p->Project3D("y");

	// Random
	h3r->GetXaxis()->SetRangeUser(0,180);
	h3r->GetZaxis()->SetRangeUser(chan,chan);
	TH1D *h1r = h3r->Project3D("y");

	a = Tgg_c12 - 10;
	b = 180;
	// Subtracted Tgg
	h1r->Sumw2();
	h1r->Scale( pa);
	h1p->Add( h1r, -1.0);
	h1p->SetLineStyle(1);
	h1p->SetLineWidth(3);
	h1p->GetXaxis()->SetRangeUser(a,b);
	Double_t max = h1p->GetMaximum();
	max *= 1.2;
	h1p->SetMaximum( max);
	sprintf( name, "#theta_{#gamma#gamma} for E_{#gamma} = %5.1f MeV", eg);
	h1p->SetTitle( name);
	h1p->Draw();
	h1p->GetXaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h1p->GetXaxis()->CenterTitle();
	l3 = new TLine( a, 0, b, 0);
	l3->SetLineWidth(1);
	l3->Draw();

	l6 = new TLine(Tgg,0,Tgg,max);
	l6->SetLineStyle(2);
	l6->SetLineWidth(4);
	l6->SetLineColor(1);
	l6->Draw();

	l7 = new TLine(Tgg_c12,0,Tgg_c12,max);
	l7->SetLineStyle(2);
	l7->SetLineWidth(4);
	l7->SetLineColor(1);
	l7->Draw();

	cdat_h = h1p->Integral();

// Simulation
	// Hydrogen
	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, egamma);
	sim_p = new TFile( name);
//	s3h = (TH3D*) sim_p->Get( "THits_TGG_ThetaCM_P");
//	s3h->GetXaxis()->SetRangeUser(0,180);		// The entire ThetaCM range
//	TH1D *s1h = s3h->Project3D("y2");
	s3h = (TH3D*) sim_p->Get( "EmPi0_TGG_ThetaCM");
	TH1D *s1h = s3h->Project3D("y2");

	// Carbon
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, egamma);
	sim_c = new TFile( name);
//	s3c = (TH3D*) sim_c->Get( "THits_TGG_ThetaCM_P");
//	s3c->GetXaxis()->SetRangeUser(0,180);		// The entire ThetaCM range
//	TH1D *s1c = s3c->Project3D("y3");
	s3c = (TH3D*) sim_c->Get( "EmPi0_TGG_ThetaCM");
	TH1D *s1c = s3c->Project3D("y3");

	s1h->SetFillStyle(3004);
	s1h->SetFillColor(4);
	s1h->SetLineColor(4);
	s1h->SetLineWidth(4);

	s1c->SetFillStyle(3005);
	s1c->SetFillColor(2);
	s1c->SetLineColor(2);
	s1c->SetLineWidth(4);

	TH1D *s1 = (TH1D*)s1h->Clone( "totsim");
	s1->Add(s1c, 1.0);

	y1 = s1->Integral();
//	y1 = s1h->Integral();

	fact = cdat_h/y1;

	s1->Scale( fact);
	s1c->Scale( fact);
	s1h->Scale( fact);

	s1h->SetFillStyle(3004);
	s1h->SetLineColor(4);
	s1h->SetLineWidth(4);

	s1h->GetXaxis()->SetRangeUser(a,b);
	sprintf( name, "Simulated #theta_{#gamma#gamma} for E_{#gamma} = %5.1f MeV",
			eg);
	s1h->SetTitle( name);
//	s1->Draw("same");
	s1c->Draw("same");
	s1h->Draw("same");

	if ( egamma < 165) t1 = 0.05;
	else t1 = 0.70;

	t2 = t1 + 0.05;
	t3 = t1 + 0.07;
	t4 = t1 + 0.14;
	t5 = t1 + 0.22;

	c = a+(b-a)*t1;
	d = a+(b-a)*t2;
	e = max*0.9;
	f = max*0.98;
	g = a+(b-a)*t3;
	h = a+(b-a)*t4;
	l4 = new TLine( c, e, d, e);
	l4->SetLineWidth(10);
	l4->SetLineColor(2);
	l4->Draw();
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	pt1->AddText("^{12}C");
	pt1->SetBorderSize(0);
	pt1->Draw();

	e = max*0.8;
	f = max*0.88;
	l4 = new TLine( c, e, d, e);
	l4->SetLineWidth(10);
	l4->SetLineColor(4);
	l4->Draw();
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	pt1->AddText("^{1}H");
	pt1->SetBorderSize(0);
	pt1->SetFillStyle(4000);
	pt1->Draw();

//	e = max*0.7;
//	f = max*0.78;
//	h = a+(b-a)*t5;
//	l4 = new TLine( c, e, d, e);
//	l4->SetLineWidth(10);
//	l4->SetLineColor(2);
//	l4->Draw();
//	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
//	pt1->AddText("^{12}C + ^{1}H");
//	pt1->SetBorderSize(0);
//	pt1->Draw();

	gStyle->SetOptStat(0);
	xx[0] = (c+d)/2;
	yy[0] = max*0.60;
	dyy[0] = max*0.04;
	dxx[0] = (b-a)*0.01;
	pm = new TGraphErrors( 1, xx, yy, dxx, dyy);
//	pm->SetMarkerSize(1);
	pm->SetLineWidth(4);
	pm->Draw();
	e = max*0.6;
	f = max*0.68;
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	pt1->AddText("Data");
	pt1->SetBorderSize(0);
	pt1->Draw();

	sprintf( name, "plots/eps/tgg_simcomp_%d.eps", egamma);
//	c1->Print( name);
}

void TggThetaCM( Int_t chan)
{
	Int_t i, egamma;
	Double_t a, b;
	Double_t pa, pa_ratio[284];
	Double_t eg, junk, energy[284], Tgg, Tgg_c12;
	Char_t name[64], cond[] = "masscorr";

	strcpy( name, "xs/chan_subt/chan_subt.out");
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> pa ;
		pa_ratio[i] = pa;
	}
	inFile.close();

	strcpy( name, "xs/tageng450.out");
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> junk >> junk >> junk >> eg ;
		energy[i] = eg;
	}
	inFile.close();

	eg = energy[chan];
	egamma = (int)(eg + 0.5);
	pa = pa_ratio[chan];

	if ( eg >= E_thr( kMP_MEV, kMPI0_MEV))
		Tgg = Tgg_Min( qT_max(eg, kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg = 180.;

	if ( eg >= E_thr( kM_C12_MEV, kMPI0_MEV))
		Tgg_c12 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg_c12 = 180.;

	a = Tgg_c12 - 10;
	b = 180;

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 700);
   c1->SetFillColor( 38);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(1,3);

	// Prompt and Random
	h3p->GetZaxis()->SetRangeUser(chan,chan);	// The appropriate tagger channel
	TH2D *h2p = h3p->Project3D("xy");
	h3r->GetZaxis()->SetRangeUser(chan,chan);	// The appropriate tagger channel
	TH2D *h2r = h3r->Project3D("xy");

	h2r->Sumw2();
	h2r->Scale( pa);
	h2p->Add( h2r, -1.0);

	c1->cd(1);
	h2p->GetXaxis()->SetRangeUser(a,b);
	h2p->SetOption("col");
	h2p->Draw();

// Simulation
	// Hydrogen
	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, egamma);
	sim_p = new TFile( name);
	s3h = (TH3D*) sim_p->Get( "THits_TGG_ThetaCM_P");
	TH2D *s2h = s3h->Project3D("xy2");
	s2h->SetTitle( "Hydrogen");
//	TH1D *s1h = s3h->Project3D("x");
//	s1h->SetTitle( "Hydrogen");

	// Carbon
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, egamma);
	sim_c = new TFile( name);
	s3c = (TH3D*) sim_c->Get( "THits_TGG_ThetaCM_P");
	TH2D *s2c = s3c->Project3D("xy3");
	s2c->SetTitle( "Carbon");
//	TH1D *s1c = s3c->Project3D("x2");
//	s1c->SetTitle( "Carbon");

	c1->cd(2);
	s2h->GetXaxis()->SetRangeUser(a,b);
	s2h->SetOption("col");
	s2h->Draw();
//	s1h->Draw();

	c1->cd(3);
	s2c->GetXaxis()->SetRangeUser(a,b);
	s2c->SetOption("col");
	s2c->Draw();
//	s1c->Draw();
}
