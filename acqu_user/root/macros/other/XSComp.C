// Macro XSComp.C
//
// This macro calculates the differential and total cross sections from the
// proton missing mass as a function of tagger channel converted to incident
// photon energy.
//
// Adapted from my old version from the June/July 2004 data with CB-TAPS round
// #1.
//
// NOTE:  The tagger channels run from 0 to 351.
//
// DLH		2009.12.10		First "New" Version for December 2008 data.
// DLH		2010.03.01		Changed to include differential cross sections.
// DLH		2010.03.18		Compton version.
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

#include "xs/msep_cuts.h"

typedef struct {

	Int_t egamma;
	Double_t energy;
	Double_t denergy;
	Double_t etag;
	Double_t d_etag;
	Double_t edet[10];
	Double_t d_edet[10];
	Double_t pa_ratio;
	Double_t xs[10];
	Double_t dxs[10];

} TData;
TData tcd[352];
TData tcdbin[11];

typedef struct {
	Int_t th;
	Double_t theta;
	Double_t dtheta;
	Double_t dom;
} TBins;
TBins tbin[10];

Double_t f_dead;
Int_t broken_ladd_ch[] = { -1};

TString prompt, random, scalers, dtscalers;
prompt = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP";
random = "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR";
scalers = "SumScalers152to503";
dtscalers = "SumScalers504to535";

// Target Full Compton Histograms
TFile full( "histograms/Compton/Full/Full.root");
TH3D *h3dp = (TH3D*)full.Get( prompt);
TH3D *h3dr = (TH3D*)full.Get( random);
TH1D *hsc = (TH1D*)full.Get( scalers);
TH1D *hdead = (TH1D*)full.Get( dtscalers);

Int_t tchan[] = { 286, 282, 277, 273, 269, 265, 261, 257, 253, 249, 245};

void InitXS()
{
	UInt_t i, eg;
	Double_t eff, deff;
	Double_t f_tagg, f_tot;
	TString name;

	gROOT->ProcessLine( ".L ReadParams.C");

	ReadTagEng( "xs/tageng855.dat");

//	ReadSubt( "xs/chan_subt/chan_subt_full.out");
	for ( i = 0; i <= 351; tcd[i++].pa_ratio = 0.0833);

	ReadDetEff( "compton");
	ReadTagEff( "xs/tag_eff/TaggEff.out");

	// Full target scaler deadtime correction.
	f_tagg = hdead->GetBinContent( 32)/hdead->GetBinContent( 31);
	f_tot = hdead->GetBinContent( 25)/hdead->GetBinContent( 26);
	f_dead = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "full tagger scalers deadtime correction = " << f_dead << endl;

	// Zero scalers for broken ladder channels.
	Int_t* lch = broken_ladd_ch;
	while ( *lch != -1) {
		hsc->SetBinContent( *lch, 0);
		lch++;
	}

	tbin[0].th = 0;
	tbin[0].theta = 0;
	tbin[0].dtheta = 0;
	tbin[0].dom = 1;
	for ( i = 1; i < 10; i++) {
		Double_t th, dth;
		th = 10*(2*i-1);
		dth = 20;
		tbin[i].th = th;
		tbin[i].theta = th;
		tbin[i].dtheta = dth;
		tbin[i].dom = (2*kPI)*(cos((th-dth/2)*kD2R)-cos((th+dth/2)*kD2R));
	}
}

void TotXS( UInt_t chan_lo, UInt_t chan_hi)
{
	Int_t i, count, chan;
	Double_t enn[352], denn[352], xsn[352], dxsn[352];
	Double_t eg, xs, junk, egg[352], xxs[352];
	TString name;

	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Calculate the channel cross sections
	// For total cross sections, the theta bin is 0 for 0-180 deg.
	for ( i = chan_lo; i <= chan_hi; ChanXS( i++, 0));

	name = "xs/compton/xs/tot_xs.out";
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	count = 0;
	for ( chan = chan_lo; chan <= chan_hi; chan++) {
		enn[count] = tcd[chan].energy;
		denn[count] = 0;
		xsn[count] = tcd[chan].xs[0];
		dxsn[count] = tcd[chan].dxs[0];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count++]
			<< endl ;
	}
	outFile.close();

	// Plot Results
	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
	gr->SetMarkerStyle( 20);
	gr->SetMarkerSize( 1.0);
	gr->SetLineWidth( 3);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle( "#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( tcd[chan_hi].energy - 10,
			tcd[chan_lo].energy + 10);
	gr->SetMaximum( 5);
	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	pt = new TLegend( 0.3, 0.7, 0.6, 0.85);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->SetTextSize( 0.04);

	pt->AddEntry( gr, "This Work", "p");

	// Theory
	name = "xs/compton/theory/totxs_params.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> eg >> xs >> junk >> junk >> junk;
		egg[i] = eg;
		xxs[i++] = xs;
	}
	inFile.close();
	count = i-1;

	// Plot Results
	gr1 = new TGraph( count, egg, xxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Lsame");

	pt->AddEntry( gr1, "Pasquini", "l");

	pt->Draw();

	name = "plots/Compton/xstot_compton";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void TotXSBin()
{
	Int_t ebin, i, count, chan;
	Double_t eg[20], deg[20], xs[20], dxs[20];
	Double_t egam, xxs, junk;
	TString name;

	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Calculate the bin cross sections
	// For total cross sections, the theta bin is 0 for 0-180 deg.
	for ( ebin = 0; ebin <= 10; BinXS( ebin++, 0));

	name = "xs/compton/xs/tot_xs_ebin.out";
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	count = 0;
	for ( ebin = 0; ebin <= 10; ebin++) {
		eg[count] = tcdbin[ebin].energy;
		deg[count] = 0;
		xs[count] = tcdbin[ebin].xs[0];
		dxs[count] = tcdbin[ebin].dxs[0];
		outFile << eg[count] << "  " << xs[count] << "  " << dxs[count++]
			<< endl ;
	}
	outFile.close();

	// Plot Results
	gr = new TGraphErrors( count-1, eg, xs, deg, dxs);
	gr->SetTitle( "Preliminary #gamma p #rightarrow #gamma p Total Cross Section");
	gr->SetMarkerStyle( 20);
	gr->SetMarkerSize( 1.0);
	gr->SetLineWidth( 3);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle( "#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( 220, 310);
	gr->SetMaximum( 3);
	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	pt = new TLegend( 0.3, 0.7, 0.6, 0.85);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->SetTextSize( 0.04);

	pt->AddEntry( gr, "This Work", "p");

	// Theory
	name = "xs/compton/theory/totxs_params.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> egam >> xxs >> junk >> junk >> junk;
		eg[i] = egam;
		xs[i++] = xxs;
	}
	inFile.close();
	count = i-1;

	// Plot Results
	gr1 = new TGraph( count, eg, xs);
	gr1->SetLineWidth( 2);
//	gr1->SetLineStyle( 2);
//	gr1->SetLineColor( 2);
	gr1->Draw( "Lsame");

	pt->AddEntry( gr1, "Pasquini", "l");

	pt->Draw();

	name = "plots/Compton/xstot_ebin";
	name.Append( ".pdf");
	c1->Print( name);
}

void ChanXS( UInt_t chan, UInt_t theta_bin) 
{
	UInt_t i, j;
	UInt_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t t, fact, emin, emax, sum, n_yld, dn_yld;
	Double_t theta_lo, theta_hi, mm_lo, mm_hi;
	Double_t theta_lo, theta_hi, chan_lo, chan_hi;
	Double_t yield_p, yield_r;
	Double_t yield, dyield, scal, dscal;
	Double_t r_fe, pa;
	TString name;

	// For brevity...
	i = chan;
	j = theta_bin;

	// Target thickness in nuclei/cm^2 for the 10-cm cell.
	t = 4.242e23;

	// Converts cm^2 to microbarn.
	fact = 1e30/t;

	// Prompt-Accidental subtraction ratio
	pa = tcd[i].pa_ratio;

	// Limits over which to integrate 3D histogram.
	//		- One tagger channel.
	// 	- CM theta.
	// 	- MMiss region.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = i+1;
	xmax = i+1;

	// Theta CM limits
	if ( j == 0)
	{
		theta_lo = 0;
		theta_hi = 180;
	}
	else
	{
		theta_lo = (j-1)*20;
		theta_hi = j*20;
	}
	ymin = h3dp->GetYaxis()->FindBin( theta_lo);
	ymax = h3dp->GetYaxis()->FindBin( theta_hi);
	ymax -= 1;	// To avoid double counting...
	
	// Missing mass integration limits
	mm_lo = 900;
	mm_hi = 940;
	zmin = h3dp->GetZaxis()->FindBin( mm_lo);
	zmax = h3dp->GetZaxis()->FindBin( mm_hi);

	// Full target
	yield_p = h3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = h3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield = yield_p - pa*yield_r;
	dyield = sqrt( yield_p - Sqr( pa)*yield_r);
	scal = hsc->GetBinContent( xmin)/f_dead;
	dscal = sqrt( scal);

	if ( ( scal != 0) && ( tcd[i].edet[j] > 0)) {

		n_yld = yield/scal;
		dn_yld = 1/scal*sqrt( Sqr( dyield) + Sqr( yield*dscal/scal));

		tcd[i].xs[j] = n_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;
		tcd[i].dxs[j] = dn_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;

		name = Form( "Ch = %d  Eg = %5.1f  edet = %4.2f  etag = %4.2f"
			"  dom = %5.3f  Y = %6.1f +/- %4.1f  S = %6.2e"
			"  xs = %4.2f +/- %4.2f", i, tcd[i].energy, tcd[i].edet[j],
			tcd[i].etag, tbin[j].dom, yield, dyield, scal,
			tcd[i].xs[j], tcd[i].dxs[j]);
		cout << name << endl;
	}
	else
	{
		tcd[i].xs[j] = 0;
		tcd[i].dxs[j] = 0;
	}
}

Bool_t BinXS( UInt_t ebin, UInt_t theta_bin) 
{
	UInt_t i, j, chan_lo, chan_hi;
	Int_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t t, fact, emin, emax, sum, n_yld, dn_yld;
	Double_t theta_lo, theta_hi, mm_lo, mm_hi;
	Double_t yield_p, yield_r;
	Double_t scal, dscal;
	Double_t r_fe, yield, dyield, pa;
	Double_t temin, temax;
	TString name;
	TFile *file;

	// Middle channel bin and bin limits
	i = tchan[ebin];
	chan_lo = i-2;
	chan_hi = i+2;

	temin = tcd[chan_hi].energy - tcd[chan_hi].denergy;
	temax = tcd[chan_lo].energy + tcd[chan_lo].denergy;

	// This is just for brevity...
	j = theta_bin;

	// Prompt-Accidental subtraction ratio
	pa = tcd[i].pa_ratio;

	file = full;
	tcdbin[ebin].energy = WeightedPar( file, "energy", chan_lo, chan_hi);
	tcdbin[ebin].etag = WeightedPar( file, "tageff", chan_lo, chan_hi);
//	tcdbin[ebin].edet[j] = WeightedPar( file, "deteff", chan_lo, chan_hi, j);
	tcdbin[ebin].edet[j] = tcd[i].edet[j];

	t = 4.242e23;		// This is in nuclei/cm^2 for the 10-cm cell.

	// The 1e30 converts cm^2 to microbarn.
	fact = 1e30/t;

	// Limits over which to integrate 3D histogram.
	//		- Multiple tagger channels.
	// 	- CM theta bin.
	// 	- MMiss region.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = chan_lo;
	xmax = chan_hi;

	// Theta CM limits
	if ( j == 0)
	{
		theta_lo = 0;
		theta_hi = 180;
	}
	else
	{
		theta_lo = (j-1)*20;
		theta_hi = j*20;
	}
	ymin = h3dp->GetYaxis()->FindBin( theta_lo);
	ymax = h3dp->GetYaxis()->FindBin( theta_hi);
	// Adjust upper limit so we don't double count...
	ymax -= 1;
	
	// Missing mass integration limits
	mm_lo = 900;
	mm_hi = msep_cut[ebin][j];
	zmin = h3dp->GetZaxis()->FindBin( mm_lo);
	zmax = h3dp->GetZaxis()->FindBin( mm_hi);

	// Full target
	yield_p = h3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = h3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield = yield_p - pa*yield_r;
	dyield = sqrt( yield_p - Sqr( pa)*yield_r);

	scal = hsc->Integral( xmin, xmax)/f_dead;
	dscal = sqrt( scal);

	if ( ( scal != 0) && ( tcdbin[ebin].edet > 0)) {

		n_yld = yield/scal;
		dn_yld = 1/scal*sqrt( Sqr( dyield) + Sqr( yield*dscal/scal));

		tcdbin[ebin].xs[j]
			= n_yld*fact/tcdbin[ebin].edet[j]/tcdbin[ebin].etag/tbin[j].dom;
		tcdbin[ebin].dxs[j]
			= dn_yld*fact/tcdbin[ebin].edet[j]/tcdbin[ebin].etag/tbin[j].dom;

		if ( j == 0) {
			name = Form( "Channel = %d (%d-%d)  Egamma = %5.1f (%5.1f-%5.1f)"
					"  edet = %4.2f  etag = %4.2f", i, chan_lo, chan_hi,
					tcdbin[ebin].energy, temin, temax, tcdbin[ebin].edet[j],
					tcdbin[ebin].etag);
			cout << name << endl;
			name = Form( "  dom = %5.3f  Y = %6.1f +/- %4.1f  S = %6.2e"
					"  xs = %4.2f +/- %4.2f", tbin[j].dom, yield, dyield, scal,
					tcdbin[ebin].xs[j], tcdbin[ebin].dxs[j]);
			cout << name << endl;
		}
		else {
			name = Form( "Channel = %d (%d-%d)  Egamma = %5.1f (%5.1f-%5.1f)", i,
					chan_lo, chan_hi, tcdbin[ebin].energy, temin, temax);
			cout << name << endl;
			name = Form( "  CM Theta = %d (%d-%d)  edet = %4.2f  etag = %4.2f",
					tbin[j].th, (int) theta_lo, (int) theta_hi, tcdbin[ebin].edet[j],
					tcdbin[ebin].etag);
			cout << name << endl;
			name = Form( "  dom = %5.3f  Y = %6.1f +/- %4.1f  S = %6.2e"
					"  xs = %5.3f +/- %5.3f", tbin[j].dom, yield, dyield, scal,
					tcdbin[ebin].xs[j], tcdbin[ebin].dxs[j]);
			cout << name << endl;
		}
		if ( tcdbin[ebin].xs[j] >= 0.001) return( kTRUE);
		else return( kFALSE);
	}
	else return( kFALSE);
}

void DiffXSBin( UInt_t ebin)
{
	UInt_t start, stop;
	UInt_t j, k, chan_lo, chan_hi;
	Double_t temin, temax;
	Double_t tageff, deteff;
	Double_t th[9], dth[9], xs[9], dxs[9];
	TString name;
	TFile* file;

	chan_lo = tchan[ebin] - 2;
	chan_hi = tchan[ebin] + 2;

	temin = tcd[chan_hi].energy - tcd[chan_hi].denergy;
	temax = tcd[chan_lo].energy + tcd[chan_lo].denergy;

	if ( ebin <= 6) start = 6;
	else start = 5;
	stop = 8;

	name = Form( "xs/compton/xs/diff_xs_ebin%d.out", ebin);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	k = 0;
	for (  j = start; j <= stop; j++) {
		if ( BinXS( ebin, j) == kTRUE)
		{
			name = Form( "%3d  %6.4f  %5.4f", tbin[j].th, tcdbin[ebin].xs[j],
				tcdbin[ebin].dxs[j]);
			outFile << name << endl;
			th[k] = tbin[j].theta;
			dth[k] = 0;
			xs[k] = tcdbin[ebin].xs[j];
			dxs[k++] = tcdbin[ebin].dxs[j];
		}
	}
	outFile.close();

	cout << "Weighted Energy = " << tcdbin[ebin].energy << endl;

	// Output results to a file
	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Plot Results
	gr = new TGraphErrors( k, th, xs, dth, dxs);
	name = Form( "#gammap#rightarrow#gammap Differential Cross Section"
			" E_{#gamma} = %5.1f (%5.1f - %5.1f) MeV", tcdbin[ebin].energy,
			temin, temax);
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
	gr->SetMinimum( 0);
	gr->Draw( "AP");

/*
	c1->cd();
	TPad *npad = new TPad( "npad", "Transparent Pad", 0, 0, 1, 1);
	npad->SetFillStyle( 4000);
	npad->Draw();
	npad->cd();

	// "Preliminary"
//	TPaveLabel *pl = new TPaveLabel( 100, 50, 380, 380, "PRELIMINARY");
	TPaveLabel *pl = new TPaveLabel( 0, 0, 1, 1, "PRELIMINARY");
	pl->SetTextAngle(30);
	pl->SetTextColor(14);
//	pl->SetTextFont(82);
	pl->SetBorderSize(0);
	pl->SetFillStyle(4000);
	pl->Draw();
*/

	name = Form( "plots/Compton/diff_xs_ebin%d", ebin);
//	name.Append( ".eps");
	name.Append( ".pdf");
//	c1->Print( name);
}

void DiffComp( UInt_t ebin)
{
	UInt_t i, ct, eg;
	Double_t th, cth, x, dx, cmth[100], dth[100], xs[100], dxs[100];
	Double_t max;
	TString name;

	eg = 200 + ebin*10;

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Theory results
	name = Form( "xs/compton/theory/compton_%d.dat", eg);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> th >> cth >> x;
		cmth[i] = cth;
		xs[i++] = x;
		cout << i-1;
		cout << " " << th;
		cout << " " << cth;
		cout << " " << x;
		cout << endl;
	}
	inFile.close();
	ct = i-1;

	if ( ebin <= 4) max = 0.18;
	else if ( ebin <= 6) max = 0.25;
	else if ( ebin <= 7) max = 0.3;
	else if ( ebin <= 8) max = 0.4;
	else if ( ebin <= 10) max = 0.5;

	// Plot Theory Results
	gr = new TGraph( ct, cmth, xs);
	name = Form( "Preliminary #gammap#rightarrow#gammap "
			"Differential Cross Section for E_{#gamma} = %d MeV", eg);
	gr->SetTitle( name);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta^{cm}_{#gamma} (deg)");
	gr->GetYaxis()->SetTitle( "d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->SetLineWidth( 2);
//	gr->SetLineColor( 2);
	gr->SetMinimum( 0.0);
	gr->SetMaximum( max);
	gr->Draw( "AL");

	// My results
	i = 0;
	name = Form( "xs/compton/xs/diff_xs_ebin%d.out", ebin);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	ifstream inFile( name);
	while( !inFile.eof()) {
		inFile >> th >> x >> dx;
		cout << i;
		cout << " " << th;
		cout << " " << x;
		cout << " " << dx;
		cout << endl;
		if ( dx != 0) {
			cmth[i] = th;
			dth[i] = 0;
			xs[i] = x;
			dxs[i++] = dx;
		}
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr1 = new TGraphErrors( ct, cmth, xs, dth, dxs);
	gr1->SetLineWidth( 2);
//	gr1->SetLineColor( 4);
//	gr1->SetMarkerColor( 4);
	gr1->SetMarkerSize( 1.0);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Psame");

	pt = new TLegend( 0.4, 0.65, 0.6, 0.80);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "Pasquini", "l");
	pt->AddEntry( gr1, "This Work", "p");
	pt->Draw();

	name = Form( "plots/Compton/dxs_comp_%d", eg);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void DiffXSChan( UInt_t i)
{
	UInt_t i, j, k, eg;
	Double_t tageff, deteff;
	Double_t tth[9], dtth[9], xxs[9], dxxs[9];
	TString name;

	eg = tcd[i].egamma;

	cout << "Energy = " << eg << endl;

	// Output results to a file
	name = Form( "xs/compton/xs/diff_xs_%d.out", eg);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	k = 0;
	for (  j = 1; j <= 9; j++) {
		ChanXS( i, j);
		name = Form( "%3d  %6.4f  %5.4f", tbin[j].th, tcd[i].xs[j],
				tcd[i].dxs[j]);
		outFile << name << endl;
		tth[k] = tbin[j].theta;
		dtth[k] = 0;
		xxs[k] = tcd[i].xs[j];
		dxxs[k++] = tcd[i].dxs[j];
	}
	outFile.close();

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Plot Results
	gr = new TGraphErrors( 9, tth, xxs, dtth, dxxs);
	name = Form( "#gammap#rightarrow#gammap Differential Cross Section"
			" E_{#gamma} = %5.1f MeV", tcd[i].energy);
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
	gr->SetMinimum( 0);
	gr->Draw( "AP");

/*
	c1->cd();
	TPad *npad = new TPad( "npad", "Transparent Pad", 0, 0, 1, 1);
	npad->SetFillStyle( 4000);
	npad->Draw();
	npad->cd();

	// "Preliminary"
//	TPaveLabel *pl = new TPaveLabel( 100, 50, 380, 380, "PRELIMINARY");
	TPaveLabel *pl = new TPaveLabel( 0, 0, 1, 1, "PRELIMINARY");
	pl->SetTextAngle(30);
	pl->SetTextColor(14);
//	pl->SetTextFont(82);
	pl->SetBorderSize(0);
	pl->SetFillStyle(4000);
	pl->Draw();
*/

	name = Form( "plots/Compton/diff_xs_chan%d", i);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}
