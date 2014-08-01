/*
 *		AsymCompton.C
 *
 *		New asymmetry macro using more intelligent organization.
 *
 *		2010.03.04		DLH		First Version
 *		2010.03.26		DLH		Modified for Compton
 *
 *		Functions:
 *			1) InitAsym
 *			2) GetHistFile
 *			3) Get1DHist
 *			4) AsymBin
 *			6) Asymmetry
 *			7) AsymComp
 *			8) AsymCompBin
 *			10) AsymComp2
 *
 */

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

typedef struct {
	Int_t egamma;
	Double_t energy;
	Double_t denergy;
	Double_t etag;
	Double_t edet[10];
} TData;
TData tcd[352];

Int_t tchan[] = { 286, 282, 277, 273, 269, 265, 261, 257, 253, 249, 245};

// Full
TFile fperpFile( "histograms/Compton/Full/Perp.root");
TFile fparaFile( "histograms/Compton/Full/Para.root");

TH3D* hS;
TH1D* hsc;
TH1D* hA;

Double_t c_pp, c_pa;

// InitAsym
//
// This just reads in the tagger energies.
//
void InitAsym()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");
}

// 
// GetHistFile
//
// Sets pointers to the scaler and random-subtracted asymmetry 3D histograms
// for a certain input file (full/empty and perp/para).  Results are NOT
// dead-time corrected.
//
void GetHistFile( TFile* file)
{
	Double_t pa;
	TString prompt, random, scalers;

	TH3D* hP;
	TH3D* hR;

	pa = 0.0833;

	prompt =
		"PhotonPhiCMProtOAMMP_v_PhotonThetaCMProtOAMMP_v_TChanPhotProtOAMMP";
	random =
		"PhotonPhiCMProtOAMMR_v_PhotonThetaCMProtOAMMR_v_TChanPhotProtOAMMR";
	scalers = "SumScalers152to503";

	hP = (TH3D*) file->Get( prompt);
	hR = (TH3D*) file->Get( random);
	hsc = (TH1D*) file->Get( scalers);

	hS = (TH3D*) hP->Clone( "sub");
	hS->Sumw2();
	hS->Add( hR, -pa);
}

//
// Get1DHist
//
// Sets the pointer to a single 1D histogram of phiCM.  Everything is weighted
// to Full Para scalers.  Histograms are corrected for dead time.
//
void Get1DHist( UInt_t bin_lo, UInt_t bin_hi, UInt_t tb_lo, UInt_t tb_hi)
{
	Double_t f_dead_f, factor;

	TFile* file;

	TH1D* h1_perp_f;
	TH1D* hsc_perp_f;
	TH1D* h1_para_f;
	TH1D* hsc_para_f;

	TH1D* diff1d;
	TH1D* sum1d;

	// Full Perp
	file = fperpFile;
	f_dead_f = DeadTimeSF( file);
	GetHistFile( file);
	hS->GetXaxis()->SetRange( bin_lo, bin_hi);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_perp_f = (TH1D*) hS->Project3D( "z");
	hsc_perp_f = (TH1D*) hsc->Clone( "sc_perp_f");

	c_pp = h1_perp_f->Integral();

	// Full Para
	file = fparaFile;
	f_dead_f = DeadTimeSF( file);
	GetHistFile( file);
	hS->GetXaxis()->SetRange( bin_lo, bin_hi);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_para_f = (TH1D*) hS->Project3D( "z3");
	hsc_para_f = (TH1D*) hsc->Clone( "sc_para_f");

	c_pa = h1_para_f->Integral();

	factor = HistSF( hsc_para_f, hsc_perp_f, bin_lo, bin_hi);
	diff1d = (TH1D*) h1_para_f->Clone( "diff1d");
	diff1d->Add( h1_perp_f, -factor);
	sum1d = (TH1D*) h1_para_f->Clone( "sum1d");
	sum1d->Add( h1_perp_f, factor);

	hA = (TH1D*) diff1d->Clone( "asym1d");
	hA->Divide( sum1d);

	h1_perp_f->Reset();
	hsc_perp_f->Reset();
	h1_para_f->Reset();
	hsc_para_f->Reset();
	diff1d->Reset();
	sum1d->Reset();

	cout << " perp = " << c_pp;
	cout << "  para = " << c_pa;
	cout << "  sum = " << c_pp + c_pa;
	cout << endl;
}

//
// AsymBin( Th_cm, Eg, Tgt)
//
// This plots p_gamma*Sigma for a specific CM theta bin and tagger-channel
// range.
//
void AsymBin( UInt_t tbin = 0, UInt_t chan_lo = 300, UInt_t chan_hi = 300)
{
	UInt_t bin_lo, bin_hi;
	UInt_t tb_lo, tb_hi;
	UInt_t theta;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	TString name;
	TFile* file;

	TH1D* asym;

	theta = 20*tbin - 10;

	// Bins start from 1, but tagger channels from 0!
	bin_lo = chan_lo + 1;
	bin_hi = chan_hi + 1;

	if ( tbin == 0)
	{
		tb_lo = 1;
		tb_hi = 9;
	}
	else
	{
		tb_lo = tbin;
		tb_hi = tbin;
	}

	Get1DHist( bin_lo, bin_hi, tb_lo, tb_hi);
	name = Form( "asym_tbin%d", tbin);
	asym = (TH1D*) hA->Clone( name);

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 700, 500);
	asym->Draw();

	file = fparaFile;
	eg = WeightedPar( file, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;
	if ( tbin == 0)
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = 0-180 deg",
				eg, deg);
	else
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = %d #pm 10 deg",
				eg, deg, theta);
	asym->SetTitle( name);
	asym->GetXaxis()->SetTitle( "#phi (deg)");
	asym->GetXaxis()->CenterTitle();
	asym->GetYaxis()->SetTitle( "p_{#gamma}#Sigma cos(2#phi + #phi_{0})");
	asym->GetYaxis()->CenterTitle();

	TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);
	asym->Fit( "f1", "REMQ");
	par[0] = f1->GetParameter( 0);
	par[1] = f1->GetParError( 0);
	par[2] = f1->GetParameter( 1);
	par[3] = f1->GetParError( 1);
	red_chisq = f1->GetChisquare()/17;
	if ( par[2] < 0) {
		par[0] *= -1;
		par[2] *= -1;
	}
	cout << par[0];
	cout << " +/- " << par[1];
	cout << "   " << par[2];
	cout << " +/- " << par[3];
	cout << "   " << red_chisq;
	cout << endl;

	pt = new TPaveText( 0.5, 0.60, 0.8, 0.85, "NDC");
	pt->SetTextAlign( 12);
	pt->SetTextSize( 0.05);
	pt->SetTextFont( 132);
	pt->SetBorderSize( 0);
	pt->SetFillColor( 0);
	name = Form( "p_{#gamma}#Sigma = %6.4f #pm %6.4f", par[0], par[1]);
	pt->AddText( name);
	name = Form( "#phi_{0} = %5.1f #pm %4.1f deg", par[2]/kD2R, par[3]/kD2R);
	pt->AddText( name);
	name = Form( "#chi^{2}/n_{d.o.f.} = %4.2f", red_chisq);
	pt->AddText( name);
	pt->Draw();

	name = Form( "plots/Compton/AsymBin_%d-%d_t%d", chan_lo, chan_hi, tbin);
	name.Append( ".pdf");
//	name.Append( ".eps");
	c1->Print( name);
}

//
// AsymEbin( Th_cm, Eg, Tgt)
//
// This plots p_gamma*Sigma for a specific CM theta bin and tagger-channel
// range.
//
void AsymEbin( UInt_t ebin = 8, UInt_t tbin = 0)
{
	UInt_t i, chan_lo, chan_hi, bin_lo, bin_hi;
	UInt_t tb_lo, tb_hi;
	UInt_t theta;
	Double_t eg, deg;
	Double_t par[4], red_chisq;
	TString name;
	TFile* file;

	TH1D* asym;

	// Middle tagger channel and bin limits
	i = tchan[ebin];
	chan_lo = i-2;
	chan_hi = i+2;
	// Bins start from 1, but tagger channels from 0!
	bin_lo = chan_lo + 1;
	bin_hi = chan_hi + 1;

	// CM Theta
	theta = 20*tbin - 10;
	if ( tbin == 0)
	{
		tb_lo = 1;
		tb_hi = 9;
	}
	else
	{
		tb_lo = tbin;
		tb_hi = tbin;
	}

	Get1DHist( bin_lo, bin_hi, tb_lo, tb_hi);
	name = Form( "asym_tbin%d", tbin);
	asym = (TH1D*) hA->Clone( name);

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 700, 500);
	asym->Draw();

	file = fparaFile;
	eg = WeightedPar( file, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;
	if ( tbin == 0)
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = 0-180 deg",
				eg, deg);
	else
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = %d #pm 10 deg",
				eg, deg, theta);
	asym->SetTitle( name);
	asym->GetXaxis()->SetTitle( "#phi (deg)");
	asym->GetXaxis()->CenterTitle();
	asym->GetYaxis()->SetTitle( "p_{#gamma}#Sigma cos(2#phi + #phi_{0})");
	asym->GetYaxis()->CenterTitle();

	TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);
	asym->Fit( "f1", "REMQ");
	par[0] = f1->GetParameter( 0);
	par[1] = f1->GetParError( 0);
	par[2] = f1->GetParameter( 1);
	par[3] = f1->GetParError( 1);
	red_chisq = f1->GetChisquare()/17;
	if ( par[2] < 0) {
		par[0] *= -1;
		par[2] *= -1;
	}
	cout << par[0];
	cout << " +/- " << par[1];
	cout << "   " << par[2];
	cout << " +/- " << par[3];
	cout << "   " << red_chisq;
	cout << endl;

	pt = new TPaveText( 0.5, 0.60, 0.8, 0.85, "NDC");
	pt->SetTextAlign( 12);
	pt->SetTextSize( 0.05);
	pt->SetTextFont( 132);
	pt->SetBorderSize( 0);
	pt->SetFillColor( 0);
	name = Form( "p_{#gamma}#Sigma = %6.4f #pm %6.4f", par[0], par[1]);
	pt->AddText( name);
	name = Form( "#phi_{0} = %5.1f #pm %4.1f deg", par[2]/kD2R, par[3]/kD2R);
	pt->AddText( name);
	name = Form( "#chi^{2}/n_{d.o.f.} = %4.2f", red_chisq);
	pt->AddText( name);
	pt->Draw();

	name = Form( "plots/Compton/AsymBin_%d-%d_t%d", chan_lo, chan_hi, tbin);
	name.Append( ".pdf");
//	name.Append( ".eps");
	c1->Print( name);
}

//
// Asymmetry( Eg, Tgt)
//
// Plots the subtracted asymmetry as a function of CM theta for a specific
// tagger-channel range.  Note that the degree of polarization is only very
// roughly estimated.
//
void Asymmetry( UInt_t chan_lo = 300, UInt_t chan_hi = 300)
{
	UInt_t i, j;
	UInt_t col;
	UInt_t bin_lo, bin_hi;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	Double_t theta[9], dtheta[9], as[9], das[9];
	TString name;
	TFile* file;

	TH1D* asym;

	Double_t polgrad;

	if ( tcd[chan_lo].energy <= 175)
		polgrad = (Linear( 140, 0.4, 175, 0.7, tcd[chan_lo].energy)
			+ Linear( 140, 0.4, 170, 0.65, tcd[chan_lo].energy))/2;
	else if ( tcd[chan_lo].energy <= 200)
		polgrad = (Linear( 175, 0.7, 200, 0.2, tcd[chan_lo].energy)
			+ Linear( 140, 0.4, 170, 0.65, tcd[chan_lo].energy))/2;
	else if ( ( tcd[chan_lo].energy <= 320) &&
			( tcd[chan_hi].energy >= 220)) polgrad = 0.2;
	else polgrad = 1;

	cout << "polgrad = " << polgrad << endl;

	col = 4;

	bin_lo = chan_lo + 1;
	bin_hi = chan_hi + 1;
	file = fparaFile;
	eg = WeightedPar( file, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	name = Form( "xs/compton/asym/asymmetry_%d-%d.out", chan_lo, chan_hi);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}

	j = 0;
	for ( i = 3; i <= 7; i++)
	{
		theta[j] = 10 + 20*i;
		dtheta[j] = 0;

		Get1DHist( bin_lo, bin_hi, i+1, i+1);
		name = Form( "asym_%d_%d", chan_lo, i+1);
		asym = (TH1D*) hA->Clone( name);
		hA->Reset();

		TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);
		asym->Fit( "f1", "REMQ0");
		par[0] = f1->GetParameter(0);
		par[1] = f1->GetParError(0);
		par[2] = f1->GetParameter(1);
		par[3] = f1->GetParError(1);
		red_chisq = f1->GetChisquare()/17;
		if ( par[2] < 0) {
			par[0] *= -1;
			par[2] *= -1;
		}

		as[j] = par[0]/polgrad;
		das[j] = par[1]/polgrad;

		name = Form( "%3d  %6.3f %4.3f  %4.2f %3.1f  %5.3f", (int) theta[j],
				as[j], das[j], par[2], par[3], red_chisq);
		cout << name << endl ;
		outFile << name << endl ;
		j++;
	}
	outFile.close();

	c1 = new TCanvas( "c1", "Asymmetry", 100, 10, 700, 500);

	gr = new TGraphErrors( j, theta, as, dtheta, das);
	name = Form( "Asymmetry for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
	gr->SetLineWidth(2);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
//	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->GetXaxis()->SetLimits( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetLineColor( col);
	gr->SetMarkerColor( col);
	gr->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

/*
	pl = new TPaveLabel( 0.70, 0.75, 0.9, 0.85, subt, "NDC");
	pl->SetTextAlign( 12);
	pl->SetTextSize( 0.8);
	pl->SetTextFont( 132);
	pl->SetTextColor( col);
	pl->SetBorderSize( 0);
	pl->SetFillColor( 0);
	pl->Draw();
*/

	name = Form( "plots/Compton/Asymmetry_%d-%d", chan_lo, chan_hi);
	name.Append( ".pdf");
//	name.Append( ".eps");
	c1->Print( name);
}

//
// AsymComp
//
// Compares asymmetry at approx. 160 MeV to ChPT, DMT2001, and the Schmidt
// measurement.
//
void AsymComp( UInt_t include = 0, Bool_t save = kFALSE)
{
	UInt_t chan_lo, chan_hi;
	Double_t eg, deg;
	Double_t theta[19], dtheta[19], as[19], das[19];
	Double_t junk;
	TString name;
	TFile* file;

	chan_lo = 299;
	chan_hi = 308;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymComp", 400, 10, 700, 500);

	// Weighting to para file.
	file = fparaFile;
	eg = WeightedPar( file, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	name = "xs/asym/asymmetry_subt_299-308.out";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 9; i++)
	{
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
	}
	inFile.close();

	gr = new TGraphErrors( 9, theta, as, dtheta, das);
	name = Form( "Asymmetry for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.2);
	gr->SetLineWidth(2);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetMaximum( 0.3);
	gr->SetMinimum( -0.2);
	gr->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	name = "xs/previous/asym_160.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 9; i++)
	{
		inFile >> theta[i] >> as[i] >> das[i];
		dtheta[i] = 0;
	}
	inFile.close();

	gr2 = new TGraphErrors( 9, theta, as, dtheta, das);
	gr2->SetMarkerStyle( 20);
	gr2->SetMarkerSize( 1.2);
	gr2->SetLineWidth( 2);
	gr2->SetMarkerColor( 4);
	gr2->SetLineColor( 4);
	if ( include >= 1) gr2->Draw( "Psame");

	name = "xs/asym/asymmetry_full_299-308.out";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 9; i++)
	{
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
	}
	inFile.close();
	gr3 = new TGraphErrors( 9, theta, as, dtheta, das);
	gr3->SetMarkerStyle( 22);
	gr3->SetMarkerSize( 1.2);
	gr3->SetLineWidth( 2);
	gr3->SetMarkerColor( 2);
	gr3->SetLineColor( 2);
	if ( include >= 2) gr3->Draw( "Psame");

	// Theory DMT2001
	name = "xs/previous/asym_dmt_160.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 19; i++)
	{
		inFile >> theta[i] >> as[i];
		dtheta[i] = 0;
	}
	inFile.close();

	gr4 = new TGraph( 19, theta, as);
	gr4->SetLineWidth( 4);
	gr4->SetLineColor( 6);
	gr4->Draw( "Lsame");

	// Theory ChPT
	name = "xs/previous/asym_chpt_160.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 19; i++)
	{
		inFile >> theta[i] >> as[i];
		dtheta[i] = 0;
	}
	inFile.close();

	gr5 = new TGraph( 19, theta, as);
	gr5->SetLineWidth( 4);
	gr5->SetLineColor( 5);
	gr5->Draw( "Lsame");

	pt = new TLegend( 0.4, 0.15, 0.6, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr4, "DMT 2001", "l");
	pt->AddEntry( gr5, "ChPT", "l");
	if ( include >= 1) pt->AddEntry( gr2, "Schmidt - TAPS", "p");
	if ( include >= 2) pt->AddEntry( gr3, "This Work - Unsubtracted", "p");
	pt->Draw();

	name = Form( "plots/Compton/AsymComp%d", include);
//	name.Append( ".pdf");
	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

//
// AsymCompBin
//
// Compares the asymmetries to Pasquini for different energy bins
// near threshold.
//
void AsymCompBin( UInt_t bin = 0, Bool_t save = kFALSE)
{
	UInt_t egamma, chan_lo, chan_hi;
	Double_t eg, deg;
	Double_t theta[19], dtheta[19], as[19], das[19];
	Double_t junk;
	TString name;
	TFile* file;

	if ( bin == 0) {
		chan_lo = 265;
		chan_hi = 273;
		egamma = 240;
	}
	else if ( bin == 1) {
		chan_lo = 257;
		chan_hi = 264;
		egamma = 260;
	}
	else if ( bin == 2) {
		chan_lo = 249;
		chan_hi = 256;
		egamma = 280;
	}

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymComp", 400, 10, 700, 500);

	file = fparaFile;
	eg = WeightedPar( file, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	// Theory - Pasquini
	name = Form( "xs/compton/theory/asym_%d.dat", egamma);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 19; i++)
	{
		inFile >> theta[i] >> as[i];
		dtheta[i] = 0;
	}
	inFile.close();

	gr = new TGraph( 19, theta, as);
	name = Form( "Asymmetry for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	gr->SetTitle( name);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
//	gr->SetMaximum( 0.25);
//	gr->SetMinimum( -0.1);
	gr->SetMaximum( 1);
	gr->SetMinimum( -1);
	gr->SetLineWidth( 4);
	gr->SetLineColor( 6);
	gr->Draw( "AL");

	// Our Data
	name = Form( "xs/compton/asym/asymmetry_%d-%d.out", chan_lo, chan_hi);
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( UInt_t i = 0; i < 4; i++)
	{
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
	}
	inFile.close();

	gr2 = new TGraphErrors( 4, theta, as, dtheta, das);
	gr2->SetMarkerStyle( 21);
	gr2->SetMarkerSize( 1.2);
	gr2->SetLineWidth( 2);
	gr2->Draw( "Psame");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	pt = new TLegend( 0.7, 0.70, 0.9, 0.85);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr2, "This Work", "p");
	pt->AddEntry( gr, "Pasquini", "l");
	pt->Draw();

	name = Form( "plots/Compton/AsymCompTheory_%d", egamma);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}
