// Macro DiffXS.C
//
// This macro calculates the differential cross section from
// the proton missing mass as a function of tagger channel
// converted to incident photon energy.
//
// 30.11.2006	DLH	Added Empty Subtraction
// 08.07.2007	DLH	Improved Empty Subtraction
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

// Target Full
//TFile full( "histograms/Full/SomeFull.root");
TFile full( "histograms/Full/Full.root");
TH3D *hf3dp = (TH3D*)full.Get( "ThetaCMP_v_MMissPi0P_v_TChanP");
TH3D *hf3dr = (TH3D*)full.Get( "ThetaCMR_v_MMissPi0R_v_TChanR");
TH1D *hfsc = (TH1D*)full.Get( "SumScalers152to503");

TH1D *hdead = (TH1D*)full.Get( "SumScalers504to535");

// Target Empty
//TFile empty( "histograms/Empty/SomeEmpty.root");
TFile empty( "histograms/Empty/Empty.root");
TH3D *he3dp = (TH3D*)empty.Get( "ThetaCMP_v_MMissPi0P_v_TChanP");
TH3D *he3dr = (TH3D*)empty.Get( "ThetaCMR_v_MMissPi0R_v_TChanR");
TH1D *hesc = (TH1D*)empty.Get( "SumScalers152to503");

Int_t egamma[352];
Double_t yield, dyield, pa_rat, r_fe;
Double_t pa_ratio[352], etag[352], d_etag[352];
Double_t energy[352], denergy[352];
Double_t f_dead;

Int_t broken_ladd_ch[] = { -1};

Double_t theta[] = { 10, 30, 50, 70, 90, 110, 130, 150, 170};
Double_t dtheta[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t dom[] = { 0.379, 1.091, 1.672, 2.051, 2.182, 2.051, 1.672, 1.091,
		0.379};

void InitXS()
{
	Int_t i, eg;
	Double_t eff, deff;
	Double_t f_tagg, f_tot;
	TString name;

	name = "xs/tageng855.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff >> deff;
		energy[i] = eff;
		denergy[i] = deff/2;
		egamma[i] = (int)(eff + 0.5);
	}
	inFile.close();

	// Random Subtraction Ratios
	// NOTE: This needs to be done for both empty and full,
	//			but I am going to leave it for now.
	name = "xs/chan_subt/chan_subt_full.out";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff;
		pa_ratio[i] = eff;
	}
	inFile.close();

	// Tagging efficiencies
	name = "xs/eff/tag_eff.out";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff >> deff;
		if ( eff != 0) {
			etag[i] = eff;
			d_etag[i] = deff;
		}
		else {
			etag[i] = -1;
			d_etag[i] = -1;
		}
	}
	inFile.close();

	f_tagg = hdead->GetBinContent( 32)/hdead->GetBinContent( 31);
	f_tot = hdead->GetBinContent( 25)/hdead->GetBinContent( 26);
	f_dead = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "tagger scalers deadtime correction = " << f_dead << endl;


	Int_t* lch = broken_ladd_ch;
	while ( *lch != -1){ 
		hfsc->SetBinContent( *lch, 0);
		hesc->SetBinContent( *lch, 0);
		lch++;
	}
}

void DiffXS( TString subt_str, Int_t chan)
{
	Int_t i, chan_lo, chan_hi;
	Double_t junk;
	Double_t sum, fscalers, dfscalers, escalers, descalers, yield_tot;
	Double_t n_yld, dn_yld, xs[9], dxs[9], theta_lo, theta_hi;
	Double_t scalers, dscalers;
	Double_t fact, t;
	Double_t emin, emax, edet[11], d_edet[11], eff, deff;
	Double_t Tgg, Tgg_c12;

	TString name;

	cout << "Energy = " << energy[chan] << endl;

	emin = Linear( 145, 933, 400, 880, energy[chan]);
	emax = Linear( 145, 945, 400, 1000, energy[chan]);

	t = 4.24e23;	// This is in nuclei/cm^2 for a 10-cm cell.

	// The 1e30 converts cm^2 to microbarn.
	fact = 1e30/t;

	name = Form( "xs/eff/det_eff_new/Eff_p_chan_%d.out", chan);
	ifstream inFile( name);
	i = 0;
	if ( inFile) {
		while( !inFile.eof()) {
			inFile >> junk >> eff >> deff ;
			edet[i] = eff;
			d_edet[i++] = deff;
		}
		inFile.close();
	}
	else {
		cerr << "Error opening file " << name << endl;
		cout << "Using constant efficiency" << endl;
		for ( i = 0; i <= 9; i++) {
			edet[i] = 0.8;
			d_edet[i] = 0.05;
		}
	}

	// Scalers
	fscalers = hfsc->GetBinContent( chan+1)*f_dead;
	dfscalers = sqrt( fscalers);
	escalers = hesc->GetBinContent( chan+1)*f_dead;
	descalers = sqrt( escalers);
	r_fe = fscalers/escalers;

	chan_lo = chan;
	chan_hi = chan;

	yield_tot = 0;
	for (  i = 0; i <= 8; i++) {

		Yield( subt_str, i+1, emin, emax, chan_lo, chan_hi);

		if ( i == 3) cout << " yield = " << yield << endl;

		if ( subt_str == "empty") {
			scalers = escalers;
			dscalers = descalers;
		}
		else {
			scalers = fscalers;
			dscalers = dfscalers;
		}

		yield_tot += yield;

		n_yld = yield/scalers;
		dn_yld = 1/scalers*sqrt( pow( dyield, 2)
				+ pow( yield*dscalers/scalers, 2));
		
		xs[i] = n_yld*fact/dom[i]/edet[i+1]/etag[chan];
		dxs[i] = dn_yld*fact/dom[i]/edet[i+1]/etag[chan];
	}
	
	// Output results to a file
	name = Form( "xs/diff_xs/diff_xs_%d.out", egamma[chan]);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( i = 0; i <= 8; i++)
		outFile << theta[i] << "  " << xs[i] << "  " << dxs[i] << endl;
	outFile.close();

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
//	c1->SetFillColor( 38);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Plot Results
	gr = new TGraphErrors( 9, theta, xs, dtheta, dxs);
	name = Form( "#gammap#rightarrowp#pi^{0} Differential Cross Section"
			" E_{#gamma} = %5.1f MeV", energy[chan]);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("#theta^{*}_{#pi^{0}} (MeV)");
	gr->GetYaxis()->SetTitle("d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
//	gr->GetXaxis()->SetRangeUser(100,402);
//	gr->SetMaximum(350);
//	gr->SetMinimum(-20);
	gr->Draw( "AP");

//	c1->cd();
//	TPad *npad = new TPad( "npad", "Transparent Pad", 0, 0, 1, 1);
//	npad->SetFillStyle( 4000);
//	npad->Draw();
//	npad->cd();
//
//	// "Preliminary"
////	TPaveLabel *pl = new TPaveLabel( 100, 50, 380, 380, "PRELIMINARY");
//	TPaveLabel *pl = new TPaveLabel( 0, 0, 1, 1, "PRELIMINARY");
//	pl->SetTextAngle(30);
//	pl->SetTextColor(14);
////	pl->SetTextFont(82);
//	pl->SetBorderSize(0);
//	pl->SetFillStyle(4000);
//	pl->Draw();

	name = Form( "plots/diff_xs_%s_%d", subt_str.Data(), egamma[chan]);
//	name.Append( ".eps");
	name.Append( ".pdf");
//	c1->Print( name);

}

// This routine calculates the yield and error for a certain bin of phase
// space.
//
void Yield( TString subt_str, Int_t thetabin, Double_t emin, Double_t emax,
		Int_t chan_lo, Int_t chan_hi)
{
	Int_t xmin, xmax, ymin, ymax;
	Double_t yield_p, yield_r; Double_t fyield, dfyield, eyield, deyield;

	pa_rat = pa_ratio[chan_lo];

	xmin = hf3dp->GetXaxis()->FindBin( chan_lo);
	xmax = hf3dp->GetXaxis()->FindBin( chan_hi);
	ymin = hf3dp->GetYaxis()->FindBin( emin);
	ymax = hf3dp->GetYaxis()->FindBin( emax);

	// Full
	yield_p = hf3dp->Integral( xmin, xmax, ymin, ymax, thetabin, thetabin);
	yield_r = hf3dr->Integral( xmin, xmax, ymin, ymax, thetabin, thetabin);
	fyield = yield_p - pa_rat*yield_r;
	dfyield = sqrt( yield_p + Sqr( pa_rat)*yield_r);

/*
	cout << "Full";
	cout << "  " << xmin;
	cout << "  " << xmax;
	cout << "  " << ymin;
	cout << "  " << ymax;
	cout << "  " << thetabin;
	cout << "  " << yield_p;
	cout << "  " << yield_r;
	cout << "  " << pa_rat;
	cout << "  " << fyield;
	cout << "  " << dfyield;
	cout << endl;
*/

	// Empty
	yield_p = he3dp->Integral( xmin, xmax, ymin, ymax, thetabin, thetabin);
	yield_r = he3dr->Integral( xmin, xmax, ymin, ymax, thetabin, thetabin);
	eyield = yield_p - pa_rat*yield_r;
	deyield = sqrt( yield_p + Sqr( pa_rat)*yield_r);

/*
	cout << "Empty";
	cout << "  " << yield_p;
	cout << "  " << yield_r;
	cout << "  " << pa_rat;
	cout << "  " << eyield;
	cout << "  " << deyield;
	cout << endl;
*/

	// Subtracted
	if ( subt_str == "full") {
		yield = fyield;
		dyield = dfyield;
	}
	else if ( subt_str == "empty") {
		yield = eyield;
		dyield = deyield;
	}
	else if ( subt_str == "subt") {
		yield = fyield - r_fe*eyield;
		dyield = sqrt( Sqr( dfyield) + Sqr( r_fe*deyield));
	}
	else {
		cout << "ERROR: must use known subtraction string";
		cout << endl;
		break;
	}

}

void DiffComp( Int_t eg)
{
	Int_t i, ct;
	Double_t th, x, dx, theta[30], dtheta[30], xs[30], dxs[30];
	Double_t xx[1], dxx[1], yy[1], dyy[1], max;
	Double_t x1, x2, y1, y2;
	TString name;

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
//	c1->SetFillColor( 38);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// My results
	i = 0;
	name = Form( "xs/diff_xs/diff_xs_%d.out", eg);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	ifstream inFile( name);
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
	name = Form( "Preliminary #gammap#rightarrowp#pi^{0} "
			"Differential Cross Section for E_{#gamma} = %d MeV", eg);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("#theta^{cm}_{#pi^{0}} (deg)");
	gr->GetYaxis()->SetTitle("d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	max *= 1.4;
	gr->GetYaxis()->SetRangeUser(0,max);
	gr->Draw( "AP");

	// Previous results
	if ( eg == 156) name = Form( "xs/previous/dxs_%d.dat", eg-1);
	else if ( eg == 166) name = Form( "xs/previous/dxs_%d.dat", eg+1);
	else name = Form( "xs/previous/dxs_%d.dat", eg);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
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

	// Plot Results
	gr1 = new TGraphErrors( ct, theta, xs, dtheta, dxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Psame");

	pt = new TLegend(0.6,0.15,0.8,0.30);
	pt->SetTextSize(0.04);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr1, "Schmidt - TAPS", "p");
	pt->Draw();

	name = Form( "plots/dxs_comp_%d", eg);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}
