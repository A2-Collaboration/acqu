// Macro XSPi0.C
//
// This macro calculates the total cross section from the proton missing mass
// as a function of tagger channel converted to incident photon energy.
//
// Adapted from my old version from the June/July 2004 data with CB-TAPS round
// #1.
//
// NOTE:  The tagger channels run from 0 to 351.
//
// DLH		2009.12.10		First "New" Version for December 2008 data.
// DLH		2010.06.11		A few minor modifications.
//

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

typedef struct {

	Int_t egamma;
	Double_t energy;
	Double_t denergy;
	Double_t etag;
	Double_t d_etag;
	Double_t edet[21];
	Double_t d_edet[21];
	Double_t pa_ratio;
	Double_t xs[21];
	Double_t dxs[21];
	Double_t ed[21];
	Double_t d_ed[21];

} TData;
TData tcd[352];

typedef struct {
	Int_t th;
	Double_t theta;
	Double_t dtheta;
	Double_t dom;
} TBins;
TBins tbin[10];

typedef struct {
	Double_t cth;
	Double_t dcth;
	Double_t dom;
} CTHBins;
CTHBins cthbin[21];

typedef struct {
	Double_t y[10];
	Double_t dy[10];
} YData;
YData yld[352];

Double_t f_dead_f, f_dead_e;
Int_t broken_ladd_ch[] = { -1};

TString prompt, random, scalers, dtscalers;

// Same histograms as for asymmetries
// Mmiss and Tgg cuts can't be changed
//
// Mmiss and Tgg cuts
//prompt = "PhiCMCut2P_v_ThetaCMCut2P_v_TChanCut2P";
//random = "PhiCMCut2R_v_ThetaCMCut2R_v_TChanCut2R";
//
// Mmiss only
//prompt = "PhiCMCut1P_v_ThetaCMCut1P_v_TChanCut1P";
//random = "PhiCMCut1R_v_ThetaCMCut1R_v_TChanCut1R";

// No Phi dependence (aren't used for asymmetries)
// Dynamic cuts on Mmiss or Emiss can be done in this macro!
//
// Missing mass for yield
//prompt = "ThetaCMP_v_MMissP_v_TChanP";
//random = "ThetaCMR_v_MMissR_v_TChanR";
//
// Missing energy for yield
//prompt = "ThetaCMP_v_EMissP_v_TChanP";
//random = "ThetaCMR_v_EMissR_v_TChanR";
prompt = "CosThetaCMP_v_EMissP_v_TChanP";
random = "CosThetaCMR_v_EMissR_v_TChanR";

// Scalers
scalers = "SumScalers152to503";
dtscalers = "SumScalers504to535";

// Target Full
TFile full( "histograms/Pi0/Full/All_20MeV.root");
TH3D *hf3dp = (TH3D*)full.Get( prompt);
TH3D *hf3dr = (TH3D*)full.Get( random);
TH1D *hfsc = (TH1D*)full.Get( scalers);
TH1D *hdeadf = (TH1D*)full.Get( dtscalers);

// Target Empty
TFile empty( "histograms/Pi0/Empty/All_20MeV.root");
TH3D *he3dp = (TH3D*)empty.Get( prompt);
TH3D *he3dr = (TH3D*)empty.Get( random);
TH1D *hesc = (TH1D*)empty.Get( scalers);
TH1D *hdeade = (TH1D*)empty.Get( dtscalers);

void InitXS()
{
	UInt_t i, eg, n;
	Double_t eff, deff;
	Double_t f_tagg, f_tot;
	Double_t lo, hi, inc, th;
	TString name;

	gROOT->ProcessLine( ".L ReadParams.C");

	ReadTagEng( "xs/tageng855.dat");

//	ReadSubt( "xs/chan_subt/chan_subt_full.out");
	for ( i = 0; i <= 351; tcd[i++].pa_ratio = 0.0833);

//	ReadDetEff( "pi0", kTRUE);
//	ReadTagEff( "xs/tag_eff/TaggEff.out", kTRUE);
	ReadDetEff( "pi0", kFALSE);
	ReadTagEff( "xs/tag_eff/tageff_sergey.dat", kFALSE);

	// Full target scaler deadtime correction.
	f_tagg = hdeadf->GetBinContent( 32)/hdeadf->GetBinContent( 31);
	f_tot = hdeadf->GetBinContent( 25)/hdeadf->GetBinContent( 26);
	f_dead_f = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "full tagger scalers deadtime correction = " << f_dead_f << endl;

	// Empty target scaler deadtime correction.
	f_tagg = hdeade->GetBinContent( 32)/hdeade->GetBinContent( 31);
	f_tot = hdeade->GetBinContent( 25)/hdeade->GetBinContent( 26);
	f_dead_e = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "empty tagger scalers deadtime correction = " << f_dead_e << endl;

	// Zero scalers for broken ladder channels.
	Int_t* lch = broken_ladd_ch;
	while ( *lch != -1) {
		hfsc->SetBinContent( *lch, 0);
		hesc->SetBinContent( *lch, 0);
		lch++;
	}

	// For dxs as a function of cm theta
	// Bin '0' is for the total xs
	tbin[0].th = 0;
	tbin[0].theta = 0;
	tbin[0].dtheta = 0;
	tbin[0].dom = 1;
	lo = 0;
	hi = 180;
	n = 9;
	inc = (hi - lo)/n;
	for ( i = 1; i < n+1; i++) {
		th = lo - inc/2 + inc*i;
		tbin[i].th = th;
		tbin[i].theta = th;
		tbin[i].dtheta = inc;
		tbin[i].dom = (2*kPI)*(cos((th-inc/2)*kD2R)-cos((th+inc/2)*kD2R));
	}

	// For dxs as a function of cm cos(theta)
	// Bin '0' is for the total xs
	cthbin[0].cth = 0;
	cthbin[0].dcth = 0;
	cthbin[0].dom = 1;
	lo = -1;
	hi = 1;
	n = 20;
	inc = (hi - lo)/n;
	for ( i = 1; i < n+1; i++) {
		cthbin[i].cth = lo - inc/2 + inc*i;
		cthbin[i].dcth = inc;
		cthbin[i].dom = 4*kPI/n;
	}
}

void TotXS( TString subt_str)
{
	Int_t i, count, chan;
	Double_t eng, x, dx;
	Double_t enn[352], denn[352], xsn[352], dxsn[352];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Calculate the channel cross sections
	// For total cross sections, the theta bin is 0 for 0-180 deg.
	for ( i = 0; i <= 351; ChanXS( subt_str, i++, 0));

	name = Form( "xs/pi0/xs/tot_xs_%s.out", subt_str.Data());
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	count = 0;
	for ( chan = 0; chan <= 351; chan++) {
		enn[count] = tcd[chan].energy;
//		denn[count] = tcd[chan].denergy;
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
//	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 20);
	gr->SetMarkerSize( 1.0);
	gr->SetLineWidth( 3);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle("#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( 100, 855);
	gr->SetMaximum( 400);
	gr->SetMinimum( -5);
	gr->Draw( "AP");

	pt = new TLegend( 0.6, 0.2, 0.85, 0.4);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->SetTextSize( 0.04);

	pt->AddEntry( gr, "This Work", "p");

	/*
	l2 = new TLine(130,-0.5,130,15);
	l2->SetLineStyle(1);
	l2->SetLineWidth(2);
	l2->Draw();
	l3 = new TLine(130,15,180,15);
	l3->SetLineStyle(1);
	l3->SetLineWidth(2);
	l3->Draw();
	l4 = new TLine(180,-0.5,180,15);
	l4->SetLineStyle(1);
	l4->SetLineWidth(2);
	l4->Draw();
	l5 = new TLine(130,-0.5,180,-0.5);
	l5->SetLineStyle(1);
	l5->SetLineWidth(2);
	l5->Draw();

	l6 = new TLine(130,0,117,185);
	l6->SetLineStyle(2);
	l6->SetLineWidth(2);
	l6->Draw();
	l7 = new TLine(180,0,250,185);
	l7->SetLineStyle(2);
	l7->SetLineWidth(2);
	l7->Draw();
*/

	// Roman's Data
	i = 0;
	name = "xs/pi0/previous/tot_xs_rl.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();

	rl = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	rl->SetMarkerColor( 2);
	rl->SetMarkerStyle( 21);
	rl->SetMarkerSize( 1.0);
	rl->SetLineWidth( 3);
	rl->SetLineColor( 2);
	rl->Draw("P");

	pt->AddEntry( rl, "Leukel - TAPS", "p");
	pt->Draw();

//	sp1 = new TPad( "sp1", "Threshold ", 0.13, 0.475, 0.51, 0.875);
	sp1 = new TPad( "sp1", "Threshold ", 0.48, 0.475, 0.86, 0.875);
	sp1->Draw();
//	sp1->SetFillColor( 40);
	sp1->SetFillColor( 29);
	sp1->SetGrid();
	sp1->GetFrame()->SetFillColor( 21);
	sp1->GetFrame()->SetBorderSize( 12);
	sp1->SetRightMargin(0.01);
	sp1->SetLeftMargin(0.04);
	sp1->SetBottomMargin(0.05);
	sp1->cd();
	gr2 = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
//	gr2->SetTitle( "Threshold Region");
	gr2->SetTitle();
//	gr2->SetMarkerColor( 4);
	gr2->SetMarkerStyle( 20);
	gr2->SetMarkerSize( 1.0);
	gr2->SetLineWidth( 3);
	gr2->GetXaxis()->SetRangeUser(130,170);
	gr2->SetMaximum(5.5);
	gr2->SetMinimum(-0.5);
	gr2->Draw( "AP");
	l1 = new TLine(145,-0.5,145,5.5);
	l1->SetLineStyle(2);
	l1->SetLineColor(2);
	l1->SetLineWidth(3);
	l1->Draw();

	// Axel's Data
	i = 0;
	name = "xs/pi0/previous/tot_xs_as.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();
	as = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	as->SetMarkerColor( 4);
	as->SetMarkerStyle( 22);
	as->SetMarkerSize( 1.0);
	as->SetLineWidth(3);
	as->SetLineColor(2);
	as->Draw("P");

	pt->AddEntry( as, "Schmidt - TAPS", "p");

	TPaveText *pt1 = new TPaveText(145.5,4.2,162.5,4.9);
	pt1->AddText(" p#pi^{0} Threshold");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(29);
	pt1->Draw();

	name = "plots/Pi0/xstot";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

//
// ChanXS
//
// Calculates channel xs or dxs from 3D histogram with the following variables:
//
// 	X: Tagger Channel
//		Y: Pi0 CM Theta
//		Z: Pi0 CM Phi
//
//	NOTE: Theta bin '0' is all theta, i.e. used for total cross section.  Target
//	string can be either "full", "empty", or "subt".
//
void ChanXS( TString subt_str, UInt_t i, UInt_t j) 
{
	Int_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t t, fact, n_yld, dn_yld;
	Double_t yield_p, yield_r;
	Double_t fyield, dfyield, fscalers, dfscalers;
	Double_t eyield, deyield, escalers, descalers;
	Double_t r_fe, yield, dyield;
	Double_t x1, x2, eg;

	t = 4.242e23;		// This is in nuclei/cm^2 for the 10-cm cell.

	// The 1e30 converts cm^2 to microbarn.
	fact = 1e30/t;

	// Limits over which to integrate 3D histogram:
	//		- One tagger channel.
	// 	- CM theta bin.
	// 	- All phi.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = i+1;
	xmax = i+1;

//	// CM Theta (this assumes 9 theta bins)
//	if ( j == 0) {
//		ymin = 1;
//		ymax = 9;
//	}
//	else {
//		ymin = j;
//		ymax = j;
//	}
//
//	// CM Phi (this assumes 18 phi bins)
//	zmin = 1;
//	zmax = 18;

	//		- Mmiss
	// 	- CM theta bin.

	// Mmiss
	eg = tcd[i].energy;
	x1 = Linear( 145, 933, 400, 880, eg);
	x2 = Linear( 145, 945, 400, 1000, eg);
	ymin = hf3dp->GetYaxis()->FindBin( x1);
	ymax = hf3dp->GetYaxis()->FindBin( x2);

	// CM theta (this assumes 9 theta bins)
	if ( j == 0) {
		zmin = 1;
		zmax = 9;
	}
	else {
		zmin = j;
		zmax = j;
	}

	// Full target
	yield_p = hf3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = hf3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	fyield = yield_p - tcd[i].pa_ratio*yield_r;
	dfyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	fscalers = hfsc->GetBinContent( xmin)/f_dead_f;
	dfscalers = sqrt( fscalers);

	// Empty
	yield_p = he3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = he3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	eyield = yield_p - tcd[i].pa_ratio*yield_r;
	deyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	escalers = hesc->GetBinContent( xmin)/f_dead_e;
	descalers = sqrt( escalers);

	if ( ( fscalers != 0) && ( escalers != 0)) {

		if ( subt_str == "full") {
			n_yld = fyield/fscalers;
			dn_yld = 1/fscalers*sqrt( Sqr( dfyield)
					+ Sqr( fyield*dfscalers/fscalers));
			yield = fyield;
			dyield = dfyield;
		}
		else if ( subt_str == "empty") {
			n_yld = eyield/escalers;
			dn_yld = 1/escalers*sqrt( Sqr( deyield)
					+ Sqr( eyield*descalers/escalers));
			yield = eyield;
			dyield = deyield;
		}
		else if ( subt_str == "subt") {
			n_yld = eyield/escalers;
			r_fe = fscalers/escalers;
			yield = fyield - r_fe*eyield;
			dyield = sqrt( Sqr( dfyield) + Sqr( r_fe*deyield));

			n_yld = yield/fscalers;
			dn_yld = 1/fscalers*sqrt( Sqr( dyield)
					+ Sqr( yield*dfscalers/fscalers));
		}
		else {
			cout << "ERROR: must use known subtraction string";
			cout << endl;
			break;
		}

		tcd[i].xs[j] = n_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;
		tcd[i].dxs[j] = dn_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;

		cout << "Chan = " << i;
		cout << "  Eg = " << tcd[i].energy;
		cout << "  edet = " << tcd[i].edet[j];
		cout << "  etag = " << tcd[i].etag;
		cout << "  dom = " << tbin[j].dom;
		cout << "  Y = " << yield;
		cout << " +/- " << dyield;
		cout << "  S = " << fscalers;
		cout << " +/- " << dfscalers;
		cout << "  xs = " << tcd[i].xs[j];
		cout << " +/- " << tcd[i].dxs[j];
		cout << endl;
	}
}

//
// ChanXS_CTH
//
// Calculates channel xs or dxs from 3D histogram with the following variables:
//
// 	X: One Tagger Channel
//		Y: EMiss
//		Z: Pi0 CM Cos Theta
//
//	NOTE: Theta bin '0' is all theta, i.e. used for total cross section.  Target
//	string can be either "full", "empty", or "subt".
//
void ChanXS_CTH( TString subt_str, UInt_t i, UInt_t j) 
{
	Int_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t t, fact, n_yld, dn_yld;
	Double_t yield_p, yield_r;
	Double_t fyield, dfyield, fscalers, dfscalers;
	Double_t eyield, deyield, escalers, descalers;
	Double_t r_fe, yield, dyield;
	Double_t eg, x1, x2;

	t = 4.242e23;		// This is in nuclei/cm^2 for the 10-cm cell.

	// The 1e30 converts cm^2 to microbarn.
	fact = 1e30/t;

	// Limits over which to integrate 3D histogram:
	//		- One tagger channel.
	//		- Emiss
	// 	- CM cos(theta) bin.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = i+1;
	xmax = i+1;

	// EMiss
	eg = tcd[i].energy;
	x1 = Linear( 144, -1, 190, -15, eg);
	x2 = Linear( 144, 5, 190, 20, eg);
	ymin = hf3dp->GetYaxis()->FindBin( x1);
	ymax = hf3dp->GetYaxis()->FindBin( x2);

	// CM cos(theta) (this assumes 20 theta bins)
	if ( j == 0) {
		zmin = 1;
		zmax = 20;
	}
	else {
		zmin = j;
		zmax = j;
	}

	// Full target
	yield_p = hf3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = hf3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	fyield = yield_p - tcd[i].pa_ratio*yield_r;
	dfyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	fscalers = hfsc->GetBinContent( xmin)/f_dead_f;
	dfscalers = sqrt( fscalers);

	// Empty
	yield_p = he3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = he3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	eyield = yield_p - tcd[i].pa_ratio*yield_r;
	deyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	escalers = hesc->GetBinContent( xmin)/f_dead_e;
	descalers = sqrt( escalers);

	if ( ( fscalers != 0) && ( escalers != 0)) {

		if ( subt_str == "full") {
			n_yld = fyield/fscalers;
			dn_yld = 1/fscalers*sqrt( Sqr( dfyield)
					+ Sqr( fyield*dfscalers/fscalers));
			yield = fyield;
			dyield = dfyield;
		}
		else if ( subt_str == "empty") {
			n_yld = eyield/escalers;
			dn_yld = 1/escalers*sqrt( Sqr( deyield)
					+ Sqr( eyield*descalers/escalers));
			yield = eyield;
			dyield = deyield;
		}
		else if ( subt_str == "subt") {
			n_yld = eyield/escalers;
			r_fe = fscalers/escalers;
			yield = fyield - r_fe*eyield;
			dyield = sqrt( Sqr( dfyield) + Sqr( r_fe*deyield));

			n_yld = yield/fscalers;
			dn_yld = 1/fscalers*sqrt( Sqr( dyield)
					+ Sqr( yield*dfscalers/fscalers));
		}
		else {
			cout << "ERROR: must use known subtraction string";
			cout << endl;
			break;
		}

		tcd[i].xs[j] = n_yld*fact/tcd[i].edet[j]/tcd[i].etag/cthbin[j].dom;
		tcd[i].dxs[j] = dn_yld*fact/tcd[i].edet[j]/tcd[i].etag/cthbin[j].dom;

		cout << "Chan = " << i;
		cout << "  Eg = " << tcd[i].energy;
		cout << "  edet = " << tcd[i].edet[j];
		cout << "  etag = " << tcd[i].etag;
		cout << "  dom = " << cthbin[j].dom;
		cout << "  Y = " << yield;
		cout << " +/- " << dyield;
		cout << "  S = " << fscalers;
		cout << " +/- " << dfscalers;
		cout << "  xs = " << tcd[i].xs[j];
		cout << " +/- " << tcd[i].dxs[j];
		cout << endl;
	}
}

void ThreshXS( TString subt_str)
{
	Int_t i, count, min, max, chan;
	Double_t eng, x, dx;
	Double_t enn[352], denn[352], xsn[352], dxsn[352];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetFillColor( 42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	min = 295;
	max = 320;
	left = 140;
	right = 210;
//	min = 305;
//	max = 313;
//	left = 135;
//	right = 155;

	count = 0;
	name = Form( "xs/pi0/xs/thresh_xs_%s.out", subt_str.Data());
	ofstream outFile( name);
	for ( chan = min; chan <= max; chan++) {
		ChanXS( subt_str, chan, 0);
		enn[count] = tcd[chan].energy;
		denn[count] = 0;
		xsn[count] = tcd[chan].xs[0];
		dxsn[count] = tcd[chan].dxs[0];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count]
			<< endl ;
		count++;
	}
	outFile.close();

	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 0.5);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle("#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( left, right);
//	gr->SetMaximum( 6);
//	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	// Label for my data
//	xx[0] = 300;
//	yy[0] = 100;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 4);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,90,345,110);
//	pt2->AddText("This Work");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	// Label for Axel's data
//	xx[0] = 300;
//	yy[0] = 60;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 6);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,50,370,70);
//	pt2->AddText("Schmidt - TAPS");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	// Axel's Data
	i = 0;
	ifstream inFile( "xs/pi0/previous/tot_xs_as.dat");
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();
	as = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	as->SetMarkerColor( 6);
	as->SetMarkerStyle( 21);
	as->SetMarkerSize( 0.7);
	as->Draw("P");

//	TPaveText *pt1 = new TPaveText(147,11.5,165,13.5);
//	pt1->AddText(" p#pi^{0} Threshold");
//	pt1->SetBorderSize(0);
//	pt1->SetFillColor(40);
//	pt1->SetFillStyle(4000);
//	pt1->Draw();

	name = "plots/Pi0/xstot_thresh";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void ThreshXS_CTH( TString subt_str)
{
	Int_t i, count, min, max, chan;
	Double_t eng, x, dx;
	Double_t left, right;
	Double_t enn[352], denn[352], xsn[352], dxsn[352];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetFillColor( 42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	min = 295;
	max = 320;
	left = 140;
	right = 210;
//	min = 305;
//	max = 313;
//	left = 135;
//	right = 155;

	count = 0;
	name = Form( "xs/pi0/xs/thresh_xs_%s.out", subt_str.Data());
	ofstream outFile( name);
	for ( chan = min; chan <= max; chan++) {
		ChanXS_CTH( subt_str, chan, 0);
		enn[count] = tcd[chan].energy;
		denn[count] = 0;
		xsn[count] = tcd[chan].xs[0];
		dxsn[count] = tcd[chan].dxs[0];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count]
			<< endl ;
		count++;
	}
	outFile.close();

	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 0.5);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle( "#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( left, right);
//	gr->SetMaximum( 6);
//	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	// Label for my data
//	xx[0] = 300;
//	yy[0] = 100;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 4);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,90,345,110);
//	pt2->AddText("This Work");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	// Label for Axel's data
//	xx[0] = 300;
//	yy[0] = 60;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 6);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,50,370,70);
//	pt2->AddText("Schmidt - TAPS");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	// Axel's Data
	i = 0;
	ifstream inFile( "xs/pi0/previous/tot_xs_as.dat");
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();
	as = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	as->SetMarkerColor( 6);
	as->SetMarkerStyle( 21);
	as->SetMarkerSize( 0.7);
	as->Draw("P");

//	TPaveText *pt1 = new TPaveText(147,11.5,165,13.5);
//	pt1->AddText(" p#pi^{0} Threshold");
//	pt1->SetBorderSize(0);
//	pt1->SetFillColor(40);
//	pt1->SetFillStyle(4000);
//	pt1->Draw();

	name = "plots/Pi0/xstot_thresh";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void DeltaXS( TString subt_str)
{
	Int_t i, count, min, max, chan;
	Double_t eng, x, dx;
	Double_t enn[200], denn[200], xsn[200], dxsn[200];
	Double_t eng1[50], deng1[50], xs1[50], dxs1[50];
	Double_t xx[10], yy[10];
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetFillColor( 42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	min = 160;
	max = 286;

	count = 0;
	name = Form( "xs/pi0/xs/delta_xs_%s.out", subt_str.Data());
	ofstream outFile( name);
	for ( chan = min; chan <= max; chan++) {
		ChanXS( subt_str, chan, 0);
		enn[count] = tcd[chan].energy;
//		denn[count] = tcd[chan].denergy;
		denn[count] = 0;
		xsn[count] = tcd[chan].xs[0];
		dxsn[count] = tcd[chan].dxs[0];
		outFile << enn[count] << "  " << xsn[count] << "  " << dxsn[count]
			<< endl ;
		count++;
	}
	outFile.close();

	gr = new TGraphErrors( count-1, enn, xsn, denn, dxsn);
	gr->SetTitle( "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Section");
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 0.5);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle("#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser(200,500);
	gr->SetMaximum( 350);
	gr->SetMinimum( 50);
	gr->Draw( "AP");

	// Roman's Data
	i = 0;
	ifstream inFile( "xs/pi0/previous/tot_xs_rl.dat");
	while( !inFile.eof()) {
		inFile >> eng >> x >> dx;
		eng1[i] = eng;
		deng1[i] = 0;
		xs1[i] = x;
		dxs1[i++] = dx;
	}
	inFile.close();

	rl = new TGraphErrors( i-1, eng1, xs1, deng1, dxs1);
	rl->SetMarkerColor( 2);
	rl->SetMarkerStyle( 21);
	rl->SetMarkerSize( 0.5);
	rl->Draw("P");

	// Label for my data
//	xx[0] = 300;
//	yy[0] = 100;
//	TPolyMarker *pm = new TPolyMarker( 1, xx, yy);
//	pm->SetMarkerColor( 4);
//	pm->SetMarkerStyle( 21);
//	pm->SetMarkerSize( 1.0);
//	pm->Draw();
//	TPaveText *pt2 = new TPaveText(305,90,345,110);
//	pt2->AddText("This Work");
//	pt2->SetBorderSize(0);
//	pt2->SetFillColor(42);
//	pt2->Draw();

	name = "plots/Pi0/xstot_delta";
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void DiffXS( TString subt_str, UInt_t i)
{
	UInt_t j, k;
	Double_t tth[9], dtth[9], xxs[9], dxxs[9];
	TString name;

	cout << "Energy = " << tcd[i].energy << endl;

	// Output results to a file
	name = Form( "xs/pi0/xs/diff_xs_%s_%d.out", subt_str.Data(), tcd[i].egamma);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	k = 0;
	for (  j = 1; j <= 9; j++) {
		ChanXS( subt_str, i, j);
		name = Form( "%3d  %6.4f  %5.4f", tbin[j].th, tcd[i].xs[j],
				tcd[i].dxs[j]);
		outFile << name << endl;
		tth[k] = tbin[j].theta;
		dtth[k] = 0;
		xxs[k] = tcd[i].xs[j];
		dxxs[k++] = tcd[i].dxs[j];
	}
	outFile.close();

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Plot Results
	gr = new TGraphErrors( 9, tth, xxs, dtth, dxxs);
	name = Form( "#gammap#rightarrowp#pi^{0} Differential Cross Section"
			" E_{#gamma} = %5.1f MeV", tcd[i].energy);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta^{*}_{#pi^{0}} (deg)");
	gr->GetYaxis()->SetTitle( "d#sigma/d#Omega (#mub/sr)");
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

	name = Form( "plots/Pi0/diff_xs_%s_%d", subt_str.Data(), tcd[i].egamma);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);

}

void DiffXS_CTH( UInt_t i, Bool_t plot = kFALSE)
{
	UInt_t j, k;
	Double_t cth[21], dcth[21], xs[21], dxs[21];
	Double_t eg, deg;
	TString name;

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	cout << "Energy = " << eg << endl;

	// Output results to a file
	name = Form( "xs/pi0/xs/dxs_cth_chan%d.out", i);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	k = 0;
	for (  j = 1; j <= 20; j++) {
		ChanXS_CTH( "subt", i, j);
		name = Form( "%6.3f  %6.4f  %5.4f", cthbin[j].cth, tcd[i].xs[j],
				tcd[i].dxs[j]);
		outFile << name << endl;
		cth[k] = cthbin[j].cth;
		dcth[k] = 0;
		xs[k] = tcd[i].xs[j];
		dxs[k++] = tcd[i].dxs[j];
	}
	outFile.close();

	if ( plot == kTRUE) {

		if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
		c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
		c1->SetGrid();
		c1->GetFrame()->SetFillColor( 21);
		c1->GetFrame()->SetBorderSize( 12);

		// Plot Results
		gr = new TGraphErrors( 20, cth, xs, dcth, dxs);
		name = Form( "#gammap#rightarrowp#pi^{0} "
				"Differential Cross Section for E_{#gamma} = %5.1f #pm %3.1f MeV",
				eg, deg);
		gr->SetTitle( name);
		gr->SetMarkerColor( 4);
		gr->SetMarkerStyle( 21);
		gr->SetLineWidth( 2);
		gr->SetLineColor( 4);
		gr->GetXaxis()->SetTitleOffset( 1.0);
		gr->GetYaxis()->SetTitleOffset( 1.0);
		gr->GetXaxis()->SetTitleSize( 0.045);
		gr->GetYaxis()->SetTitleSize( 0.045);
		gr->GetXaxis()->SetTitle( "cos #theta_{#pi}");
		gr->GetYaxis()->SetTitle( "d#sigma/d#Omega_{#pi} (#mub/sr)");
		gr->GetXaxis()->SetLabelSize( 0.03);
		gr->GetYaxis()->SetLabelSize( 0.03);
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->CenterTitle();
		gr->GetXaxis()->SetRangeUser( -1, 1);
		gr->SetMinimum( 0);
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

		name = Form( "plots/Pi0/dxs_cth_chan%d", i);
		name.Append( ".pdf");
		c1->Print( name);
	}

}

void DiffComp( Int_t eg)
{
	Int_t i, ct;
	Double_t th, x, dx, tth[30], dtth[30], xxs[30], dxxs[30];
	Double_t xx[1], dxx[1], yy[1], dyy[1], max;
	Double_t x1, x2, y1, y2;
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
//	c1->SetFillColor( 38);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// My results
	i = 0;
	name = Form( "xs/pi0/xs/diff_xs_%d.out", eg);
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
		tth[i] = th;
		dtth[i] = 0;
		xxs[i] = x;
		if ( max < xxs[i]) max = xxs[i];
		dxxs[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr = new TGraphErrors( ct, tth, xxs, dtth, dxxs);
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
	gr->GetXaxis()->SetTitle( "#theta^{cm}_{#pi^{0}} (deg)");
	gr->GetYaxis()->SetTitle( "d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	max *= 1.4;
	gr->GetYaxis()->SetRangeUser( 0, max);
	gr->Draw( "AP");

	// Previous results
//	if ( eg == 156) name = Form( "xs/pi0/previous/dxs_%d.dat", eg-1);
//	else if ( eg == 166) name = Form( "xs/pi0/previous/dxs_%d.dat", eg-1);
//	else if ( eg == 168) name = Form( "xs/pi0/previous/dxs_%d.dat", eg-1);
	if ( eg == 156) name = "xs/pi0/previous/dxs_155.3MeV.dat";
	else if ( eg == 166) name = "xs/pi0/previous/dxs_165.6MeV.dat";
	else if ( eg == 168) name = "xs/pi0/previous/dxs_167.8MeV.dat";
	else name = Form( "xs/pi0/previous/dxs_%d.dat", eg);
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
		tth[i] = th;
		dtth[i] = 0;
		xxs[i] = x;
		dxxs[i++] = dx;
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr1 = new TGraphErrors( ct, tth, xxs, dtth, dxxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Psame");

	pt = new TLegend( 0.6, 0.15, 0.8, 0.30);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr1, "Schmidt - TAPS", "p");
	pt->Draw();

	name = Form( "plots/Pi0/dxs_comp_%d", eg);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

UInt_t chan_rl[] = { 257, 261, 265, 269, 273, 277, 282, 286};

void DiffXS_CTH_all()
{
	UInt_t i;

	for ( i = 290; i <= 308; i++) DiffXS_CTH( i, kFALSE);

	for ( i = 0; i < 8; i++) DiffXS_CTH( chan_rl[i], kFALSE);

}

void DiffCompSP_all()
{
	UInt_t i;

	for ( i = 0; i < 8; i++) DiffCompSP( chan_rl[i]);

	for ( i = 290; i <= 308; i++) DiffCompSP( i);

}

void DiffCompSP( UInt_t chan)
{
	UInt_t i;
	Double_t costh, cs, dcs;
	Double_t cth[30], dcth[30], xs[30], dxs[30];
	Double_t eg, deg;
	TString name;

	eg = tcd[chan].energy;
	deg = tcd[chan].denergy;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	// Legend
	if ( chan > 306) pt = new TLegend( 0.3, 0.15, 0.5, 0.28);
	else pt = new TLegend( 0.6, 0.15, 0.8, 0.28);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);

	// My results
	name = Form( "xs/pi0/xs/dxs_cth_chan%d.out", chan);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> costh >> cs >> dcs;
		cth[i] = costh;
		dcth[i] = 0;
		xs[i] = cs;
		dxs[i++] = dcs;
	}
	inFile.close();
	gr = new TGraphErrors( i-1, cth, xs, dcth, dxs);
	name = Form( "#gammap#rightarrowp#pi^{0} "
			"Differential Cross Section for E_{#gamma} = %5.1f #pm %3.1f MeV",
			eg, deg);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetXaxis()->SetTitleSize( 0.04);
	gr->GetYaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleSize( 0.04);
	gr->GetXaxis()->SetTitle( "cos #theta");
	gr->GetYaxis()->SetTitle( "d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRangeUser( -1, 1);
	gr->SetMinimum( 0);
	gr->Draw( "AP");
	c1->Update();
	Double_t max = gPad->GetUymax();
	max *= 1.1;
	gr->SetMaximum( max);
	pt->AddEntry( gr, "Hornidge", "p");

	// Sergey's Results
	name = Form( "xs/pi0/xs_sergey/dxs_cth_chan%d.dat", chan);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( i = 0; i < 8; i++) name.ReadLine( inFile);
	i = 0;
	while( !inFile.eof()) {
		name.ReadLine( inFile);
		sscanf( name.Data(), "%lf%lf%lf", &costh, &cs, &dcs);
		cth[i] = costh;
		dcth[i] = 0;
		xs[i] = cs;
		dxs[i++] = dcs;
	}
	inFile.close();
	gr1 = new TGraphErrors( i-1, cth, xs, dcth, dxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Psame");
	pt->AddEntry( gr1, "Prakhov", "p");

	// TAPS results
	Double_t eg_taps;
	if ( ( chan >= 299) && ( chan <= 308)) {
		if ( chan == 308) eg_taps = 146.2;
		else if ( chan == 307) eg_taps = 148.5;
		else if ( chan == 306) eg_taps = 150.8;
		else if ( chan == 305) eg_taps = 154.2;
		else if ( chan == 304) eg_taps = 156.5;
		else if ( chan == 303) eg_taps = 158.8;
		else if ( chan == 302) eg_taps = 161.0;
		else if ( chan == 301) eg_taps = 163.3;
		else if ( chan == 300) eg_taps = 165.6;
		else if ( chan == 299) eg_taps = 167.8;
		name = Form( "xs/pi0/previous/dxs_%5.1fMeV.dat", eg_taps);
		ifstream inFile( name);
		if ( !inFile.is_open()) {
			cout << "Error opening file ";
			cout << name;
			cout << endl;
			break;
		}
		i = 0;
		while( !inFile.eof()) {
			inFile >> costh >> cs >> dcs;
			cth[i] = cos( costh*kD2R);
			dcth[i] = 0;
			xs[i] = cs;
			dxs[i++] = dcs;
		}
		inFile.close();
		gr2 = new TGraphErrors( i-1, cth, xs, dcth, dxs);
		gr2->SetMarkerSize( 1.2);
		gr2->SetLineWidth( 2);
		gr2->SetMarkerStyle( 22);
		gr2->Draw( "Psame");
		pt->AddEntry( gr2, "Schmidt - TAPS", "p");
	}
	else if ( ( chan == 257) || ( chan == 261) || ( chan == 265)
			|| ( chan == 269) || ( chan == 273) || ( chan == 277)
			|| ( chan == 282) || ( chan == 285)) {
		if ( chan == 257) eg_taps = 270;
		else if ( chan == 261) eg_taps = 260;
		else if ( chan == 265) eg_taps = 250;
		else if ( chan == 269) eg_taps = 240;
		else if ( chan == 273) eg_taps = 230;
		else if ( chan == 277) eg_taps = 220;
		else if ( chan == 282) eg_taps = 210;
		else if ( chan == 285) eg_taps = 202;
		name = Form( "xs/pi0/previous/dxs_%3.0f.dat", eg_taps);
		ifstream inFile( name);
		if ( !inFile.is_open()) {
			cout << "Error opening file ";
			cout << name;
			cout << endl;
			break;
		}
		i = 0;
		while( !inFile.eof()) {
			inFile >> costh >> cs >> dcs;
			cth[i] = cos( costh*kD2R);
			dcth[i] = 0;
			xs[i] = cs;
			dxs[i++] = dcs;
		}
		inFile.close();
		gr2 = new TGraphErrors( i-1, cth, xs, dcth, dxs);
		gr2->SetMarkerSize( 1.2);
		gr2->SetLineWidth( 2);
		gr2->SetMarkerStyle( 22);
		gr2->Draw( "Psame");
		pt->AddEntry( gr2, "Leukel - TAPS", "p");
	}

	pt->Draw();

	// PDF Output
	name = Form( "plots/Pi0/dxs_comp_chan%d", chan);
	name.Append( ".pdf");
	c1->Print( name);
}

void ThreshYield( UInt_t chan_lo = 300, UInt_t chan_hi = 314,
		Bool_t thr_flag = kFALSE)
{
	UInt_t count, chan;
	Double_t left, right;
	Double_t eg[352], deg[352];
	Double_t fy[352], dfy[352];
	Double_t ey[352], dey[352];
	Double_t sy[352], dsy[352];
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Threshold Yield", 200, 10, 700, 500);
//	c1->SetFillColor( 42);
//	c1->SetGrid();
//	c1->GetFrame()->SetFillColor( 21);
//	c1->GetFrame()->SetBorderSize( 12);

	// Full
	count = 0;
	for ( chan = chan_lo; chan <= chan_hi; chan++) {
		ChanYield( "full", chan, 0);
		eg[count] = tcd[chan].energy;
		deg[count] = tcd[chan].denergy;
		fy[count] = yld[chan].y[0];
		dfy[count] = yld[chan].dy[0];
		count++;
	}

	// Scaled Empty
	count = 0;
	for ( chan = chan_lo; chan <= chan_hi; chan++) {
		ChanYield( "empty", chan, 0);
		eg[count] = tcd[chan].energy;
		deg[count] = tcd[chan].denergy;
		ey[count] = yld[chan].y[0];
		dey[count] = yld[chan].dy[0];
		count++;
	}

	// Subtracted
	count = 0;
	for ( chan = chan_lo; chan <= chan_hi; chan++) {
		ChanYield( "subt", chan, 0);
		eg[count] = tcd[chan].energy;
		deg[count] = tcd[chan].denergy;
		sy[count] = yld[chan].y[0];
		dsy[count] = yld[chan].dy[0];
		count++;
	}

	// Full
	gr = new TGraphErrors( count-1, eg, fy, deg, dfy);
	gr->SetTitle( "Target Window Contributions");
	gr->SetMarkerColor( 4);
	gr->SetLineColor( 4);
	gr->SetLineWidth( 2);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle( "#pi^{0} Yield (counts)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetXaxis()->SetRange( chan_lo-1, chan_hi+1);
//	gr->SetMaximum( 6);
//	gr->SetMinimum( -2000);
	gr->Draw( "AP");
	c1->Update();
	Double_t low = gPad->GetUymin();
	Double_t left = gPad->GetUxmin();
	Double_t right = gPad->GetUxmax();

	l2 = new TLine( left, 0, right, 0);
//	l2->SetLineWidth( 2);
	l2->Draw();

	// Scaled Empty
	gr1 = new TGraphErrors( count-1, eg, ey, deg, dey);
	gr1->SetMarkerColor( 2);
	gr1->SetLineColor( 2);
	gr1->SetLineWidth( 2);
	gr1->SetMarkerStyle( 21);
	gr1->SetMarkerSize( 1);
	gr1->Draw( "Psame");

	// Subtracted
	gr2 = new TGraphErrors( count-1, eg, sy, deg, dsy);
	gr2->SetMarkerStyle( 21);
	gr2->SetMarkerSize( 1);
	gr2->SetLineWidth( 2);
	gr2->Draw( "Psame");

	pt = new TLegend( 0.2, 0.6, 0.4, 0.75);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->SetTextSize( 0.04);

	pt->AddEntry( gr, "Full", "p");
	pt->AddEntry( gr1, "Scaled Empty", "p");
	pt->AddEntry( gr2, "Subtracted", "p");

	pt->Draw();

	if ( thr_flag == kTRUE) {
		// p pi0 Threshold
		l = new TArrow( 144.7, -1000, 144.7, low, 0.02);
		l->SetLineWidth( 2);
		l->SetFillColor( 1);
		l->Draw();
		pl = new TPaveLabel( 145, -700, 148, -1200, "p Threshold");
		pl->SetBorderSize( 0);
		pl->SetFillColor( 0);
		pl->SetTextAlign( 12);
		pl->Draw();

		// 12C pi0 Threshold
		l1 = new TArrow( 135.8, -1000, 135.8, low, 0.02);
		l1->SetLineWidth( 2);
		l1->SetFillColor( 1);
		l1->Draw();
		pl1 = new TPaveLabel( 132, -700, 135.6, -1200, "^{12}C Threshold");
		pl1->SetBorderSize( 0);
		pl1->SetFillColor( 0);
		pl1->SetTextAlign( 32);
		pl1->Draw();
	}

//	name = Form( "plots/Pi0/thresh_yield_%d-%d", chan_lo, chan_hi);
	name = Form( "plots/Pi0/thresh_yield_cut_%d-%d", chan_lo, chan_hi);
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);

	for ( UInt_t i = 0; i < count; i++) {
		if ( fy[i] != 0) {
			name = Form( "%5.1f  %5.3f", eg[i], ey[i]/fy[i]);
			cout << name << endl;
		}
	}
}

//
// ChanYield
//
// Calculates channel yield (only!) from 3D histogram with the following
// variables:
//
// 	X: Tagger Channel
//		Y: Proton Missing Energy
//		Z: Pi0 CM Theta
//
//	NOTE: Theta bin '0' is all theta, i.e. used for total cross section.  Target
//	string can be either "full" or "empty".
//
//	Note also that this 3D histogram has no cuts on it other than that for the
//	pi0 invariant mass.
//
void ChanYield( TString subt_str, UInt_t i, UInt_t j) 
{
	Int_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t yield_p, yield_r;
	Double_t fyield, dfyield, fscalers, dfscalers;
	Double_t eyield, deyield, escalers, descalers;
	Double_t r_fe, yield, dyield;
	Double_t x1, x2, eg;

	// Limits over which to integrate 3D histogram:
	//		- One tagger channel.
	// 	- EMiss window
	// 	- CM theta bin.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = i+1;
	xmax = i+1;

	// EMiss
	eg = tcd[i].energy;
	x1 = Linear( 144, -1, 190, -15, eg);
	x2 = Linear( 144, 5, 190, 20, eg);
//	x1 = -5;
//	x2 = 25;
	ymin = hf3dp->GetYaxis()->FindBin( x1);
	ymax = hf3dp->GetYaxis()->FindBin( x2);

	// CM Theta (this assumes 9 theta bins)
	if ( j == 0) {
		zmin = 1;
		zmax = 9;
	}
	else {
		zmin = j;
		zmax = j;
	}

	// Full target
	yield_p = hf3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = hf3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	fyield = yield_p - tcd[i].pa_ratio*yield_r;
	dfyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	fscalers = hfsc->GetBinContent( xmin)/f_dead_f;
	dfscalers = sqrt( fscalers);

	// Empty
	yield_p = he3dp->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	yield_r = he3dr->Integral( xmin, xmax, ymin, ymax, zmin, zmax);
	eyield = yield_p - tcd[i].pa_ratio*yield_r;
	deyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	escalers = hesc->GetBinContent( xmin)/f_dead_e;
	descalers = sqrt( escalers);

	if ( ( fscalers != 0) && ( escalers != 0)) {

		r_fe = fscalers/escalers;

		if ( subt_str == "full") {
			yield = fyield;
			dyield = dfyield;
		}
		else if ( subt_str == "empty") {
			yield = r_fe*eyield;
			dyield = r_fe*deyield;
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

		yld[i].y[j] = yield;
		yld[i].dy[j] = dyield;

		cout << "Chan = " << i;
		cout << "  " << yield_p;
		cout << "  " << yield_r;
		cout << "  " << tcd[i].pa_ratio;
		cout << "  Eg = " << tcd[i].energy;
		cout << "  Y = " << yld[i].y[j];
		cout << " +/- " << yld[i].dy[j];
		cout << endl;
	}
}
