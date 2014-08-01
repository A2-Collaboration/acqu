gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

// Empty Data
TFile data( "histograms/ARH_empty.root");

// EmPi0
TH2D *d2p = (TH2D*)data.Get( "THits_EmPi0_P");	// Prompt
TH2D *d2r = (TH2D*)data.Get( "THits_EmPi0_R");	// Random

// Tgg
//TH3D *h3p = (TH3D*)data.Get( "THits_TGG_ThetaCM_P");	// Prompt
//TH3D *h3r = (TH3D*)data.Get( "THits_TGG_ThetaCM_R");	// Random
TH3D *h3p = (TH3D*)data.Get( "THits_TGG_ThetaCM_C_P");	// Prompt
TH3D *h3r = (TH3D*)data.Get( "THits_TGG_ThetaCM_C_R");	// Random

void EmPi0( Int_t chan)
{
	Int_t i, bin, egamma[284], eg;
	Int_t ch1, ch2, ch3, ch4;
	Int_t cdat_c, cdat_h, csim_c, csim_h;
	Int_t bmin, bmax, count;
	Double_t fact, eng, junk;
	Double_t x1, x2, x3, a, b, max, min;
	Double_t pa, pa_ratio[284];
	Double_t y1, y2, lo, hi;
	Double_t c, d, e, f, g, h;
	Double_t chi_sq, b1, db1, b2, db2;
	Double_t xx[1], yy[1], dxx[1], dyy[1];
	Char_t name[64], cond[] = "masscorr";

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetFillColor( 38);
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
		inFile >> i >> junk >> junk >> junk >> eng ;
		egamma[i] = (int)(eng + 0.5);
	}
	inFile.close();

	eg = egamma[chan];
	pa = pa_ratio[chan];

// Data
	// Prompt EmPi0
	bin = d2p->GetYaxis()->FindBin( chan);
	TH1D *d1p = d2p->ProjectionX( "Data", bin, bin);

	// Random EmPi0
	bin = d2r->GetYaxis()->FindBin( chan);
	TH1D *d1r = d2r->ProjectionX( "Random", bin, bin);

	x1 = Linear( 145, 933, 400, 880, eg);
	x2 = Linear( 145, 945, 400, 1000, eg);
	x3 = Linear( 145, 918, 400, 880, eg);

	a = Linear( 145, 900, 400, 600, eg);
	b = Linear( 145, 950, 400, 1100, eg);

	// Subtracted EmPi0
	d1r->Sumw2();
	d1r->Scale( pa);
	d1p->Add( d1r, -1.0);
	d1p->SetLineStyle(1);
	d1p->SetLineWidth(3);
	d1p->GetXaxis()->SetRangeUser(a,b);
	sprintf( name, "M_{miss} for E_{#gamma} = %3d MeV", eg);
	d1p->SetTitle( name);

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

	cdat_c = d1p->Integral( ch3, ch4);

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

	l1 = new TLine( x3, 0, x3, max);
	l1->SetLineColor(2);
	l1->SetLineStyle(2);
	l1->SetLineWidth(3);

	l2 = new TLine( x1, 0, x1, max);
	l2->SetLineColor(2);
	l2->SetLineStyle(2);
	l2->SetLineWidth(3);

	l5 = new TLine( x2, 0, x2, max);
	l5->SetLineColor(2);
	l5->SetLineStyle(2);
	l5->SetLineWidth(3);

// Simulation
	// Carbon
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	sim_c = new TFile( name);
//	s2c = (TH2D*) sim_c->Get( "THits_EmPi0_P");
//	TH1D *s1c = s2c->ProjectionX( "carb");
	s3c = (TH3D*) sim_c->Get( "EmPi0_TGG_ThetaCM");
	TH1D *s1c = s3c->Project3D("z");

	csim_c = s1c->Integral( ch3, ch4);

	s1c->SetFillStyle(3004);
	s1c->SetFillColor(4);
	s1c->SetLineColor(4);
	s1c->SetLineWidth(4);
	fact = (double) cdat_c/csim_c;
	s1c->Scale( fact);

	d1p->Draw("E1P");
	s1c->Draw("same");
	l1->Draw();
	l2->Draw();
	l3->Draw();
	l5->Draw();

	c = a+(b-a)*0.08;
	d = a+(b-a)*0.13;
	e = max*0.9;
	f = max*0.98;
	g = a+(b-a)*0.13;
	h = a+(b-a)*0.21;
	l4 = new TLine( c, e, d, e);
	l4->SetLineWidth(10);
	l4->SetLineColor(3);
	l4->Draw();
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	pt1->AddText("^{12}C");
	pt1->SetBorderSize(0);
	pt1->Draw();

	xx[0] = (c+d)/2;
	yy[0] = max*0.80;
	dyy[0] = max*0.04;
	dxx[0] = (b-a)*0.01;
	pm = new TGraphErrors( 1, xx, yy, dxx, dyy);
//	pm->SetMarkerSize(1);
	pm->SetLineWidth(4);
	pm->Draw();
	e = max*0.8;
	f = max*0.88;
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	pt1->AddText("Data");
	pt1->SetBorderSize(0);
	pt1->Draw();

//	sprintf( name, "plots/eps/Empi0_empty_%d.eps", egamma);
//	c1->Print( name);

	// Carbon windows for EmPi0 Cut-off
	lo = Linear( 145, 918, 400, 880, eg);
	hi = Linear( 145, 933, 400, 880, eg);
//	lo = Linear( 137, 921, 145, 911, eg);
//	hi = Linear( 137, 931, 145, 932, eg);

	bmin = d1p->GetXaxis()->FindBin( lo);
	bmax = d1p->GetXaxis()->FindBin( hi);
	count = 0;
	chi_sq = 0;
	for ( i = bmin; i < = bmax; i++)
	{
		b1 = d1p->GetBinContent(i);
		db1 = d1p->GetBinError(i);

		b2 = s1c->GetBinContent(i);
		db2 = s1c->GetBinError(i);
      
		//Calculates the chi-sq for this bin and crystal
		if ( ( db1 != 0) || ( db2 != 0))
		{
			chi_sq += Sqr( ( b1 - b2)/sqrt( Sqr( db1) + Sqr( db2)));
			count++;
		}
	}
	chi_sq /= count - 1;

	e = max*0.35;
	f = max*0.58;
	g = a+(b-a)*0.1;
	h = a+(b-a)*0.3;
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	sprintf( name, "#frac{#chi^{2}}{n_{dof}} = %5.2f\n", chi_sq);
	pt1->AddText( name);
	pt1->SetBorderSize(0);
	pt1->Draw();
}

void Tgg( Int_t chan)
{
	Int_t i, bin, egamma;
	Int_t d_counts, s_counts;
	Int_t min, max, count;
	Double_t Tgg_c12;
	Double_t fact, junk, energy[284], eg;
	Double_t a, b, ym, ym_d, ym_s;
	Double_t pa, pa_ratio[284];
	Double_t b1, db1, b2, db2;
	Double_t chi_sq;
	Double_t e, f, g, h;
	Char_t name[64], cond[] = "masscorr";

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700,
			500);
   c1->SetFillColor( 38);
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

	if ( eg >= E_thr( kM_C12_MEV, kMPI0_MEV))
		Tgg_c12 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg_c12 = 180.;


// Empty Data

	// Prompt
	h3p->GetXaxis()->SetRangeUser(0,180);		// The entire ThetaCM range
	h3p->GetZaxis()->SetRangeUser(chan,chan);	// The appropriate tagger channel
	TH1D *h1p = h3p->Project3D("y");

	// Random
	h3r->GetXaxis()->SetRangeUser(0,180);
	h3r->GetZaxis()->SetRangeUser(chan,chan);
	TH1D *h1r = h3r->Project3D("y");

//	a = Tgg_c12 - 40;
	a = 90;
	b = 180;

	// Subtracted
	h1r->Sumw2();
	h1r->Scale( pa);
	h1p->Add( h1r, -1.0);
	h1p->SetLineStyle(1);
	h1p->SetLineWidth(3);
	h1p->GetXaxis()->SetRangeUser(a,b);
	sprintf( name, "#theta_{#gamma#gamma} for E_{#gamma} = %5.1f MeV", eg);
	h1p->SetTitle( name);
	c1->cd(1);
	h1p->Draw();
	l3 = new TLine( a, 0, b, 0);
	l3->SetLineWidth(1);
	l3->Draw();

	d_counts = h1p->Integral();

// Simulation with carbon

	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, egamma);
	sim_c = new TFile( name);
//	s3c = (TH3D*) sim_c->Get( "THits_TGG_ThetaCM_P");
//	s3c->GetXaxis()->SetRangeUser(0,180);		// The entire ThetaCM range
//	TH1D *s1c = s3c->Project3D("y3");
	s3c = (TH3D*) sim_c->Get( "EmPi0_TGG_ThetaCM");
	TH1D *s1c = s3c->Project3D("y");

	s_counts = s1c->Integral();

	s1c->SetTitle();
	s1c->SetFillStyle(3005);
	s1c->SetFillColor(3);
	s1c->SetLineColor(3);
	s1c->SetLineWidth(4);
	fact = (double) d_counts/s_counts;
	s1c->Scale( fact);
	s1c->Draw( "same");

	Double_t ym_d = h1p->GetMaximum();
	Double_t ym_s = s1c->GetMaximum();
	if ( ym_d > ym_s ) ym = ym_d;
	else ym = ym_s;

	ym *= 1.2;
	h1p->SetMaximum( ym);

	l1 = new TLine(Tgg_c12,0,Tgg_c12,ym);
	l1->SetLineStyle(2);
	l1->SetLineWidth(4);
	l1->SetLineColor(6);
	l1->Draw();

	gStyle->SetOptStat(0);

//	sprintf( name, "plots/eps/Tgg_empty_%d.eps", egamma);
//	c1->Print( name);

	min = h1p->GetXaxis()->FindBin( a);
	max = h1p->GetXaxis()->FindBin( b);
	count = 0;
	chi_sq = 0;
	for ( i = min; i < =max; i++)
	{
		b1 = h1p->GetBinContent(i);
		db1 = h1p->GetBinError(i);

		b2 = s1c->GetBinContent(i);
		db2 = s1c->GetBinError(i);
      
		//Calculates the chi-sq for this bin and crystal
		if ( ( db1 != 0) || ( db2 != 0))
		{
			chi_sq += Sqr( ( b1 - b2)/sqrt( Sqr( db1) + Sqr( db2)));
			count++;
		}
	}
	chi_sq /= count - 1;

	e = ym*0.85;
	f = ym*0.98;
	g = a+(b-a)*0.6;
	h = a+(b-a)*0.9;
	pt1 = new TPaveText( g, e*0.95, h, f*0.95);
	sprintf( name, "#frac{#chi^{2}}{n_{dof}} = %5.2f\n", chi_sq);
	pt1->AddText( name);
	pt1->SetBorderSize(0);
	pt1->Draw();

}
