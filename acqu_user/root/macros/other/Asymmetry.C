/*
 *		Asymmetry.C
 *
 *		New asymmetry macro using more intelligent organization.
 *
 *		2010.03.04		DLH		First Version
 *		2010.06.14		DLH		Added channel polarization and cleaned up
 *										tagger channel vs. tagger bin, i.e. removed
 *										all instances of tagger bin (bin_lo and bin_hi).
 *
 *		Functions:
 *
 *			1) InitAsym:			Reads in relevant tagger energies and weighted
 *										channel polarization.
 *
 *			2) GetHistFile:		Gets 3D random-subtracted full/empty/subt
 *										histogram and tagger scalers.
 *
 *			3) Get1DHist:			Gets 1D full/empty/subt phi histogram for
 *										specific tagger-channel range and theta bin.
 *
 *			4) AsymBin:				Plots p_g*Sigma for full, empty, OR subt, eg
 *										range, and theta bin.
 *
 *			5) AsymBinSubt:		Plots full, empty, AND subt p_g*Sigma for eg
 *										range and theta bin.
 *
 *			6) Asymmetry:			Plots the full, empty OR subt Sigma for an eg
 *										range as a function of theta.
 *
 *			7) AsymComp:			Compares Sigma at roughly 160 MeV to ChPT,
 *										DMT2001, and the Schmidt measurement.
 *
 *			8) AsymCompBin:		Compares the asymmetry to DMT2001 and ChPT for
 *										different energy bins near threshold.
 *
 *			9) AsymCompDelta:		Compares the asymmetry to DMT2001 one energy
 *										bin in the Delta region.
 *
 *			10) AsymCompTAPS:		Compares the full AND subt asymmetries to the
 *										TAPS asymmetries for one bin in the threshold
 *										OR Delta region.
 *
 *			11) EmptyAsymCheck:	Check of p_gamma*Sigma for empty target data.
 *										This was necessary because the 2D histogram
 *										Divide() method did not seem to work properly.
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
	Double_t d_etag;
//	Double_t etag[4];
//	Double_t d_etag[4];
	Double_t edet[10];
	Double_t pol[4];
} TData;
TData tcd[352];

typedef struct {
	UInt_t ch_lo;
	UInt_t ch_hi;
	Double_t eg;
	Double_t deg;
} SPBins;
SPBins spbin[18];

typedef struct {
	Double_t th;
	Double_t dth;
	Double_t cth;
	Double_t dcth;
	Double_t as;
	Double_t das;
} ASYM;
ASYM Asym[9];

// Full
TFile fperpFile( "histograms/Pi0/Full/Perp_20MeV.root");
TFile fparaFile( "histograms/Pi0/Full/Para_20MeV.root");

// Empty
TFile eperpFile( "histograms/Pi0/Empty/Perp_20MeV.root");
TFile eparaFile( "histograms/Pi0/Empty/Para_20MeV.root");

TH3D* hS;
TH1D* hsc;
TH1D* hA;

// InitAsym
//
// This just reads in the tagger energies.
//
void InitAsym()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");
	ReadPolData( "xs/pol_data.dat");
	ReadSergeyBins( "xs/sergey_bins.dat");
//	ReadTaggEffW( "xs/eff/tageffw.dat");
}

// 
// GetHistFile
//
// Sets pointers to the scaler and random-subtracted asymmetry 3D histograms
// for a certain input file (full/empty and perp/para).  Results are NOT
// dead-time corrected.
//
void GetHistFile( TFile* file, Bool_t cthflag = kFALSE)
{
	Double_t pa;
	TString prompt, random, scalers;

	TH3D* hP;
	TH3D* hR;

	pa = 0.0833;

	if ( cthflag == kFALSE ) {
//		prompt = "PhiCMCut2P_v_ThetaCMCut2P_v_TChanCut2P";
//		random = "PhiCMCut2R_v_ThetaCMCut2R_v_TChanCut2R";
		prompt = "PhiCMCut1P_v_ThetaCMCut1P_v_TChanCut1P";
		random = "PhiCMCut1R_v_ThetaCMCut1R_v_TChanCut1R";
	}
	else {
		prompt = "PhiCMCut1P_v_CosThetaCMCut1P_v_TChanCut1P";
		random = "PhiCMCut1R_v_CosThetaCMCut1R_v_TChanCut1R";
	}
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
// Sets the pointer to a single 1D histogram of phiCM.  for either "full",
// "empty", or "subt", and everything is weighted to Full Para scalers.
// Histograms are corrected for dead time and weighted polarization.
//
void Get1DHist( UInt_t chan_lo, UInt_t chan_hi, TString subt, UInt_t tb_lo,
		UInt_t tb_hi, Bool_t cthflag = kFALSE)
{
	UInt_t fpe, fpa, epe, epa;
	Double_t f_dead[4], polgrad[4], te_corr[4], factor;
	Double_t fact[3], factor;

	// Tags for full perp, full para, empty perp, and empty para
	fpe = 0;
	fpa = 1;
	epe = 2;
	epa = 3;

	TFile* file;

	TH1D* h1_f_perp;
	TH1D* hsc_f_perp;
	TH1D* h1_f_para;
	TH1D* hsc_f_para;
	TH1D* h1_e_perp;
	TH1D* hsc_e_perp;
	TH1D* h1_e_para;
	TH1D* hsc_e_para;

	TH1D* diff1d;
	TH1D* sum1d;

	// Full Perp
	file = fperpFile;
	GetHistFile( file, cthflag);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_f_perp = (TH1D*) hS->Project3D( "z");
	hsc_f_perp = (TH1D*) hsc->Clone( "sc_f_perp");
	f_dead[fpe] = DeadTimeSF( hsc_f_perp);
	polgrad[fpe] = WeightedPol( hsc_f_perp, chan_lo, chan_hi, fpe);
//	te_corr[fpe] = TaggEffCorr( hsc_f_perp, chan_lo, chan_hi, fpe);

	// Empty Perp
	file = eperpFile;
	GetHistFile( file, cthflag);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_e_perp = (TH1D*) hS->Project3D( "z2");
	hsc_e_perp = (TH1D*) hsc->Clone( "sc_e_perp");
	f_dead[epe] = DeadTimeSF( hsc_e_perp);
	polgrad[epe] = WeightedPol( hsc_e_perp, chan_lo, chan_hi, epe);
//	te_corr[epe] = TaggEffCorr( hsc_e_perp, chan_lo, chan_hi, epe);

	// Corrected Perp
	//
	// NOTE: we correct Empty to Full using:
	// 			1) sum scalers
	// 			2) deadtime scalers
	// 			3) weighted polarization
	if ( subt == "subt") {
		factor = HistSF( hsc_f_perp, hsc_e_perp, chan_lo, chan_hi);
		fact[0] = f_dead[fpe]/f_dead[epe];
		fact[1] = polgrad[fpe]/polgrad[epe];
//		fact[2] = te_corr[fpe]/te_corr[epe];
		factor *= fact[0]; 
		factor /= fact[1]; 

		h1_f_perp->Add( h1_e_perp, -factor);

//		cout << "perp";
//		cout << " " << fact[0];
//		cout << " " << fact[1];
//		cout << " " << fact[2];
//		cout << " " << factor;
//		cout << endl;
	}

	// Full Para
	file = fparaFile;
	GetHistFile( file, cthflag);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_f_para = (TH1D*) hS->Project3D( "z3");
	hsc_f_para = (TH1D*) hsc->Clone( "sc_f_para");
	f_dead[fpa] = DeadTimeSF( hsc_f_para);
	polgrad[fpa] = WeightedPol( hsc_f_para, chan_lo, chan_hi, fpa);
//	te_corr[fpa] = TaggEffCorr( hsc_f_para, chan_lo, chan_hi, fpa);

	// Empty Para
	file = eparaFile;
	GetHistFile( file, cthflag);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_e_para = (TH1D*) hS->Project3D( "z4");
	hsc_e_para = (TH1D*) hsc->Clone( "sc_e_para");
	f_dead[epa] = DeadTimeSF( hsc_e_para);
	polgrad[epa] = WeightedPol( hsc_e_para, chan_lo, chan_hi, epa);
//	te_corr[epa] = TaggEffCorr( hsc_e_para, chan_lo, chan_hi, epa);

	// Corrected Para
	//
	// NOTE: we correct Empty to Full using:
	// 			1) sum scalers
	// 			2) deadtime scalers
	// 			3) weighted polarization
	if ( subt == "subt") {
		factor = HistSF( hsc_f_para, hsc_e_para, chan_lo, chan_hi);
		fact[0] = f_dead[fpa]/f_dead[epa];
		fact[1] = polgrad[fpa]/polgrad[epa];
//		fact[2] = te_corr[fpa]/te_corr[epa];
		factor *= fact[0]; 
		factor /= fact[1]; 

		h1_f_para->Add( h1_e_para, -factor);

//		cout << "para";
//		cout << " " << fact[0];
//		cout << " " << fact[1];
//		cout << " " << fact[2];
//		cout << " " << factor;
//		cout << endl;
	}

	// Generate Difference, Sum, and Difference/Sum
	//
	// NOTE: we correct Perp to Para using:
	// 			1) sum scalers
	// 			2) deadtime scalers
	// 			3) weighted polarization
	if ( subt != "empty") {

		factor = HistSF( hsc_f_para, hsc_f_perp, chan_lo, chan_hi);
		fact[0] = f_dead[fpa]/f_dead[fpe];
		fact[1] = polgrad[fpa]/polgrad[fpe];
//		fact[2] = te_corr[fpa]/te_corr[fpe];
		factor *= fact[0]; 
		factor /= fact[1]; 

		diff1d = (TH1D*) h1_f_para->Clone( "diff1d");
		diff1d->Add( h1_f_perp, -factor);
		sum1d = (TH1D*) h1_f_para->Clone( "sum1d");
		sum1d->Add( h1_f_perp, factor);

//		cout << "subt/full";
//		cout << " " << fact[0];
//		cout << " " << fact[1];
//		cout << " " << fact[2];
//		cout << " " << factor;
//		cout << endl;
	}
	else {

		factor = HistSF( hsc_e_para, hsc_e_perp, chan_lo, chan_hi);
		fact[0] = f_dead[epa]/f_dead[epe];
		fact[1] = polgrad[epa]/polgrad[epe];
//		fact[2] = te_corr[epa]/te_corr[epe];
		factor *= fact[0]; 
		factor /= fact[1]; 

		diff1d = (TH1D*) h1_e_para->Clone( "diff1d");
		diff1d->Add( h1_e_perp, -factor);
		sum1d = (TH1D*) h1_e_para->Clone( "sum1d");
		sum1d->Add( h1_e_perp, factor);

//		cout << "empty";
//		cout << " " << fact[0];
//		cout << " " << fact[1];
//		cout << " " << fact[2];
//		cout << " " << factor;
//		cout << endl;
	}
	hA = (TH1D*) diff1d->Clone( "asym1d");
	hA->Divide( sum1d);

	h1_f_perp->Reset();
	hsc_f_perp->Reset();
	h1_e_perp->Reset();
	hsc_e_perp->Reset();
	h1_f_para->Reset();
	hsc_f_para->Reset();
	h1_e_para->Reset();
	hsc_e_para->Reset();
	diff1d->Reset();
	sum1d->Reset();
}

//
// AsymBin( Th_cm, Eg, Tgt)
//
// This plots the individual full, empty, or subtracted p_gamma*Sigma for a
// specific CM theta bin and tagger-channel range.
//
void AsymBin( UInt_t tbin = 0, UInt_t chan_lo = 300, UInt_t chan_hi = 300,
		TString subt = "subt")
{
	UInt_t tb_lo, tb_hi;
	UInt_t theta;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	TString name;

	TH1D* asym;
	TH1D* hscal;

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

	Double_t polgrad;

	Get1DHist( chan_lo, chan_hi, subt, tb_lo, tb_hi);
	name = Form( "asym_tbin%d", tbin);
	asym = (TH1D*) hA->Clone( name);

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 700, 500);
	asym->Draw();

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
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

	name = Form( "plots/Pi0/AsymBin_%s_%d-%d_t%d", subt.Data(), chan_lo, chan_hi,
			tbin);
	name.Append( ".pdf");
//	name.Append( ".eps");
	c1->Print( name);
}

//
// AsymBinSubt( Th_cm, Eg)
//
// Plots the full, empty, AND subtracted p_gamma*Sigma all at once for a theta
// CM bin and tagger-channel range.
//
void AsymBinSubt( UInt_t tbin = 0, UInt_t chan_lo = 300, UInt_t chan_hi = 300)
{
	UInt_t col;
	UInt_t tb_lo, tb_hi;
	UInt_t theta;
	Double_t eg, deg;
	Double_t par[4];
	TString subt, name;

	TH1D* hscal;

	TFile* file;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 900, 400);
	c1->Divide(3,1);

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

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);

	for ( UInt_t i = 0; i < 3; i++)
	{
		c1->cd( i+1);

		if ( i == 0) {
			subt = "full";
			col = 4;
		}
		else if ( i == 1) {
			subt = "empty";
			col = 2;
		}
		else if ( i == 2) {
			subt = "subt";
			col = 1;
		}

		Get1DHist( chan_lo, chan_hi, subt, tb_lo, tb_hi);
		name = Form( "asym_%s_%d_%d", subt.Data(), chan_lo, tbin);
		if ( asym = (TH1D*)(gROOT->FindObject( name)) ) delete asym;
		TH1D *asym = (TH1D*) hA->Clone( name);

		if ( i == 0) {
			if ( tbin == 0)
				name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = 0-180 deg",
						eg, deg);
			else
				name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV  #theta = %d"
						"#pm 10 deg", eg, deg, theta);
		}
		else name = "";

		asym->SetTitle( name);
		asym->GetXaxis()->SetTitle( "#phi (deg)");
		asym->GetXaxis()->CenterTitle();
		asym->GetYaxis()->SetTitle( "p_{#gamma}#Sigma cos(2#phi + #phi_{0})");
		asym->GetYaxis()->CenterTitle();
		asym->SetLineColor( col);
		asym->SetMarkerColor( col);
		asym->Draw();
		f1->SetLineColor( col);
		asym->Fit( "f1", "REMQ");

		par[0] = f1->GetParameter( 0);
		par[1] = f1->GetParError( 0);
		par[2] = f1->GetParameter( 1);
		par[3] = f1->GetParError( 1);
		par[4] = f1->GetChisquare()/17;
		if ( par[2] < 0) {
			par[0] *= -1;
			par[2] *= -1;
		}
		cout << par[0];
		cout << " +/- " << par[1];
		cout << "   " << par[2];
		cout << " +/- " << par[3];
		cout << "   " << par[4];
		cout << endl;

		pt = new TPaveText( 0.45, 0.74, 0.74, 0.89, "NDC");
		pt->SetTextAlign( 12);
		pt->SetTextSize( 0.05);
		pt->SetTextFont( 132);
		pt->SetBorderSize( 0);
		pt->SetFillColor( 0);
		name = Form( "p_{#gamma}#Sigma = %6.4f #pm %6.4f", par[0], par[1]);
		pt->AddText( name);
		name = Form( "#phi_{0} = %5.1f #pm %4.1f deg", par[2]/kD2R,
			par[3]/kD2R);
		pt->AddText( name);
		name = Form( "#chi^{2}/n_{d.o.f.} = %4.2f", par[4]);
		pt->AddText( name);
		pt->Draw();

		pl = new TPaveLabel( 0.31, 0.15, 0.5, 0.22, subt, "NDC");
		pl->SetTextAlign( 12);
		pl->SetTextSize( 0.8);
		pl->SetTextFont( 132);
		pl->SetTextColor( col);
		pl->SetBorderSize( 0);
		pl->SetFillColor( 0);
		pl->Draw();
	}

	name = Form( "plots/Pi0/AsymBinSubt_%d-%d_t%d", chan_lo, chan_hi, tbin);
	name.Append( ".pdf");
//	name.Append( ".eps");
	c1->Print( name);
}

//
// Asymmetry( Eg, Tgt)
//
// Plots the full, empty, or subtracted asymmetry as a function of CM theta for
// a specific tagger-channel range.
//
void Asymmetry( UInt_t chan_lo = 300, UInt_t chan_hi = 300,
		TString subt = "subt")
{
	UInt_t i;
	UInt_t col;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	Double_t theta[9], dtheta[9], as[9], das[9];
	Double_t polgrad;
	TString name;

	TH1D* asym;
	TH1D* hscal;

	if ( subt == "full") col = 4;
	else if ( subt == "empty") col = 2;
	else if ( subt == "subt") col = 1;

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	// Everything is weighted to para...
	if ( subt != "empty") polgrad = WeightedPol( hscal, chan_lo, chan_hi, 1);
	else polgrad = WeightedPol( hscal, chan_hi, chan_hi, 3);

	cout << "polgrad = " << polgrad << endl;

	name = Form( "xs/pi0/asym/asymmetry_%s_%d-%d.out", subt.Data(), chan_lo,
			chan_hi);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}

	for ( i = 0; i <= 8; i++)
		name = Form( "asym_%s_%d_%d", subt.Data(), chan_lo, i+1);

	for ( i = 0; i <= 8; i++)
	{
		theta[i] = 10 + 20*i;
		dtheta[i] = 0;

		Get1DHist( chan_lo, chan_hi, subt, i+1, i+1);
		name = Form( "asym_%s_%d_%d", subt.Data(), chan_lo, i+1);
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

		as[i] = par[0]/polgrad;
		das[i] = par[1]/polgrad;

		name = Form( "%3d  %6.3f %4.3f  %4.2f %3.1f  %5.3f", (int) theta[i],
				as[i], das[i], par[2], par[3], red_chisq);
		cout << name << endl ;
		outFile << name << endl ;
	}
	outFile.close();

	c1 = new TCanvas( "c1", "Asymmetry", 100, 10, 700, 500);

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
	gr->SetLineColor( col);
	gr->SetMarkerColor( col);
	gr->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	pl = new TPaveLabel( 0.70, 0.75, 0.9, 0.85, subt, "NDC");
	pl->SetTextAlign( 12);
	pl->SetTextSize( 0.8);
	pl->SetTextFont( 132);
	pl->SetTextColor( col);
	pl->SetBorderSize( 0);
	pl->SetFillColor( 0);
	pl->Draw();

	name = Form( "plots/Pi0/Asymmetry_%s_%d-%d", subt.Data(), chan_lo,
			chan_hi);
	name.Append( ".pdf");
//	name.Append( ".eps");
	c1->Print( name);
}

//
// AsymmetryCTH( bin, draw)
//
// Plots the full, empty, or subtracted asymmetry as a function of
// cos(theta) for a specific tagger-channel range.
//
// Bins run from 0-17.
//
void AsymmetryCTH( UInt_t bin = 0, Bool_t draw = kFALSE)
{
	UInt_t i, chan_lo, chan_hi, cthbins;
	UInt_t col;
	Double_t eg, deg;
	Double_t par[6], red_chisq;
	Double_t cth[9], dcth[9], as[9], das[9];
	Double_t polgrad;
	TString subt, name;

	if ( bin < 3) {
		chan_lo = 304 - 2*bin;
		chan_hi = chan_lo + 1;
	}
	else {
		chan_lo = 302 - bin;
		chan_hi = chan_lo;
	}

	if ( bin < 4) cthbins = 7;
	else cthbins = 9;

	TH1D* asym;
	TH1D* hscal;

	subt = "subt";
	if ( subt == "full") col = 4;
	else if ( subt == "empty") col = 2;
	else if ( subt == "subt") col = 1;

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	// Everything is weighted to para...
	if ( subt != "empty") polgrad = WeightedPol( hscal, chan_lo, chan_hi, 1);
	else polgrad = WeightedPol( hscal, chan_hi, chan_hi, 3);

	cout << "polgrad = " << polgrad << endl;

	name = Form( "xs/pi0/asym/asymmetry_cth_%d.out", bin);
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}

	for ( i = 0; i <= cthbins; i++)
		name = Form( "asym_%s_%d_%d", subt.Data(), chan_lo, i+1);

	Double_t interval = 2/Double_t( cthbins);
	UInt_t cthb_lo, cthb_hi;
	for ( i = 0; i < cthbins; i++)
	{
		cth[i] = -1 + interval*(i + 0.5);
		dcth[i] = 0;

		cthb_lo = 1 + 63/cthbins*i;
		cthb_hi = 63/cthbins*(i + 1);

//		name = Form( "%1d %6.3f %2d %2d\n", i, cth[i], cthb_lo[i], cthb_hi[i]);
//		cout << name;

		Get1DHist( chan_lo, chan_hi, subt, cthb_lo, cthb_hi, kTRUE);
		name = Form( "asym_%s_%d_%d", subt.Data(), chan_lo, i+1);
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

		as[i] = par[0]/polgrad;
		das[i] = par[1]/polgrad;

		Asym[i].cth = cth[i];
		Asym[i].dcth = dcth[i];
		Asym[i].as = as[i];
		Asym[i].das = das[i];

		name = Form( "%6.3f  %6.3f %4.3f  %4.2f %3.1f  %5.3f", cth[i],
				as[i], das[i], par[2], par[3], red_chisq);
		cout << name << endl ;
		outFile << name << endl ;
	}
	outFile.close();

	if ( draw == kTRUE) {

		c1 = new TCanvas( "c1", "Asymmetry", 100, 10, 700, 500);

		gr = new TGraphErrors( cthbins, cth, as, dcth, das);
		name = Form( "Asymmetry for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
		gr->SetTitle( name);
		gr->SetMarkerStyle( 21);
		gr->SetMarkerSize( 1.2);
		gr->SetLineWidth(2);
		gr->GetXaxis()->SetTitleOffset( 1.1);
		gr->GetYaxis()->SetTitleOffset( 0.8);
		gr->GetYaxis()->SetTitleSize( 0.05);
		gr->GetXaxis()->SetTitle( "cos #theta");
		gr->GetYaxis()->SetTitle( "#Sigma");
		gr->GetXaxis()->SetLabelSize( 0.03);
		gr->GetYaxis()->SetLabelSize( 0.03);
		gr->GetXaxis()->SetRangeUser( -1, 1);
		if ( bin > 1) gr->SetMinimum( 0);
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->CenterTitle();
		gr->SetLineColor( col);
		gr->SetMarkerColor( col);
		gr->Draw( "AP");

		l1 = new TLine( -1, 0, 1, 0);
		l1->Draw();

//		pl = new TPaveLabel( 0.70, 0.75, 0.9, 0.85, subt, "NDC");
//		pl->SetTextAlign( 12);
//		pl->SetTextSize( 0.8);
//		pl->SetTextFont( 132);
//		pl->SetTextColor( col);
//		pl->SetBorderSize( 0);
//		pl->SetFillColor( 0);
//		pl->Draw();

		name = Form( "plots/Pi0/Asymmetry_cth_%d", bin);
		name.Append( ".pdf");
//		name.Append( ".eps");
		c1->Print( name);
	}
}

//
// AsymComp
//
// Compares asymmetry at approx. 160 MeV to ChPT, DMT2001, and the Schmidt
// measurement.
//
void AsymComp( UInt_t include = 0, Bool_t save = kFALSE)
{
	UInt_t i, chan_lo, chan_hi;
	Double_t eg, deg;
	Double_t theta[100], dtheta[100], as[100], das[100];
	Double_t junk;
	TString name;

	TH1D* hscal;

//	chan_lo = 299;
	chan_lo = 300;
	chan_hi = 308;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymComp", 400, 10, 700, 500);

	// Weighting to para file.
	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	name = "xs/pi0/asym/asymmetry_subt_299-308.out";
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
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
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

	// Previous Data - Axel TAPS
	name = "xs/pi0/previous/asym_160.dat";
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
	gr2->SetMarkerSize( 1.2);
	gr2->SetLineWidth( 2);
	gr2->SetMarkerColor( 4);
	gr2->SetLineColor( 4);
	if ( include >= 1) gr2->Draw( "Psame");

	name = "xs/pi0/asym/asymmetry_full_299-308.out";
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
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr3 = new TGraphErrors( i-1, theta, as, dtheta, das);
	gr3->SetMarkerStyle( 22);
	gr3->SetMarkerSize( 1.2);
	gr3->SetLineWidth( 2);
	gr3->SetMarkerColor( 2);
	gr3->SetLineColor( 2);
	if ( include >= 2) gr3->Draw( "Psame");

	// Theory DMT2001
	name = "xs/pi0/theory/asym_dmt_160MeV.dat";
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
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr4 = new TGraph( i-1, theta, as);
	gr4->SetLineWidth( 4);
	gr4->SetLineColor( 6);
	gr4->Draw( "Lsame");

	// Theory ChPT
	name = "xs/pi0/theory/asym_chpt_160MeV.dat";
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
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr5 = new TGraph( i-1, theta, as);
	gr5->SetLineWidth( 4);
	gr5->SetLineColor( 5);
	gr5->Draw( "Lsame");

	pt = new TLegend( 0.4, 0.15, 0.6, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "CB-TAPS Preliminary", "p");
	pt->AddEntry( gr4, "DMT 2001", "l");
	pt->AddEntry( gr5, "ChPT", "l");
	if ( include >= 1) pt->AddEntry( gr2, "Schmidt - TAPS", "p");
	if ( include >= 2) pt->AddEntry( gr3, "CB-TAPS Preliminary - Unsubtracted",
			"p");
	pt->Draw();

	name = Form( "plots/Pi0/AsymComp%d", include);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

//
// AsymCompBin
//
// Compares the asymmetries to DMT2001 and ChPT for different energy bins
// near threshold.
//
void AsymCompBin( UInt_t bin = 0, Bool_t save = kFALSE)
{
	UInt_t i, egamma, chan_lo, chan_hi;
	Double_t eg, deg;
	Double_t theta[100], dtheta[100], as[100], das[100];
	Double_t junk;
	Double_t x[2], y[2];
	TString name;

	TH1D* hscal;

	if ( bin == 0) {
		chan_lo = 305;
		chan_hi = 308;
		egamma = 150;
	}
	else if ( bin == 1) {
		chan_lo = 301;
		chan_hi = 304;
		egamma = 160;
	}
	else if ( bin == 2) {
		chan_lo = 297;
		chan_hi = 300;
		egamma = 170;
	}
	else if ( bin == 3) {
		chan_lo = 293;
		chan_hi = 296;
		egamma = 180;
	}
	else if ( bin == 4) {
		chan_lo = 251;
		chan_hi = 255;
		egamma = 280;
	}

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymComp", 400, 10, 700, 500);

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	name = Form( "xs/pi0/asym/asymmetry_subt_%d-%d.out", chan_lo, chan_hi);
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
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
	name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	gr->SetTitle( name);
	gr->SetMarkerStyle( 21);
	gr->SetMarkerSize( 1.5);
	gr->SetLineWidth( 3);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 0.80);
	gr->GetYaxis()->SetTitleSize( 0.06);
	gr->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr->GetYaxis()->SetTitle( "#Sigma");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	if ( bin != 4) gr->SetMaximum( 0.4);
	else gr->SetMaximum( 0.60);
	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	if ( bin != 4) tex = new TLatex( 15, 0.075, "PRELIMINARY");              
	else tex = new TLatex( 15, 0.35, "PRELIMINARY");              
	tex->SetTextColor( 17);
	tex->SetTextSize( 0.15);
//	tex->SetTextAngle( 15);
//	tex->SetLineWidth( 2);
	tex->Draw();

	gr->Draw( "P");
	gPad->RedrawAxis();

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	if ( bin == 1) {
		// Previous Data - Axel TAPS
		name = "xs/pi0/previous/asym_160.dat";
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

		gr5 = new TGraphErrors( i-1, theta, as, dtheta, das);
		gr5->SetMarkerStyle( 20);
		gr5->SetMarkerSize( 1.5);
		gr5->SetLineWidth( 3);
		gr5->SetMarkerColor( 4);
		gr5->SetLineColor( 4);
		gr5->Draw( "Psame");
	}

	// Theory DMT2001
	name = Form( "xs/pi0/theory/asym_dmt_%dMeV.dat", egamma);
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
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr2 = new TGraph( i-1, theta, as);
	gr2->SetLineWidth( 4);
	gr2->SetLineColor( 2);
	gr2->Draw( "Lsame");

	// Theory
	if ( bin != 4) {
		name = Form( "xs/pi0/theory/asym_chpt_%dMeV.dat", egamma);
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
			dtheta[i] = 0;
			i++;
		}
		inFile.close();

		gr3 = new TGraph( i-1, theta, as);
		gr3->SetLineWidth( 4);
		gr3->SetLineColor( 6);
		gr3->SetLineStyle( 2);
		gr3->Draw( "Lsame");

		x[0] = 0.54;
		x[1] = 0.88;
		y[0] = 0.67;
		y[1] = 0.88;
	}
	else {
		// Previous Data -- Leukel TAPS
		name = "xs/pi0/previous/asym_280.dat";
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

		gr4 = new TGraphErrors( i-1, theta, as, dtheta, das);
		gr4->SetMarkerStyle( 20);
		gr4->SetMarkerSize( 1.5);
		gr4->SetLineWidth( 3);
		gr4->SetMarkerColor( 4);
		gr4->SetLineColor( 4);
		gr4->Draw( "Psame");
	
		x[0] = 0.3;
		x[1] = 0.7;
		y[0] = 0.3;
		y[1] = 0.5;
	}
	
	pt = new TLegend( x[0], y[0], x[1], y[1]);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "Hornidge CB-TAPS", "p");
	if ( bin == 1) pt->AddEntry( gr5, "Schmidt TAPS", "p");
	if ( bin == 4) pt->AddEntry( gr4, "Leukel TAPS", "p");
	pt->AddEntry( gr2, "DMT 2001", "l");
	if ( bin != 4) pt->AddEntry( gr3, "ChPT", "l");
	pt->Draw();

	name = Form( "plots/Pi0/AsymCompTheory_%d", egamma);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

//
// AsymCompDelta
//
// Compares the asymmetries to DMT2001 one energy bin in the Delta region,
//	tagger channels 251-255.
//
void AsymCompDelta( Bool_t save = kFALSE)
{
	UInt_t i, chan_lo, chan_hi;
	Double_t eg, deg;
	Double_t theta[100], dtheta[100], as[100], das[100];
	Double_t junk;
	TString name;

	TH1D* hscal;

	chan_lo = 251;
	chan_hi = 255;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymComp", 400, 10, 700, 500);

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	name = Form( "xs/pi0/asym/asymmetry_subt_%d-%d.out", chan_lo, chan_hi);
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
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr = new TGraphErrors( i-1, theta, as, dtheta, das);
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
//	gr->SetMaximum( 0.3);
	gr->SetMinimum( -0.1);
	gr->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	// Previous Data -- Leukel TAPS
	name = "xs/pi0/previous/asym_280.dat";
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
	gr2->SetMarkerSize( 1.2);
	gr2->SetLineWidth( 2);
	gr2->SetMarkerColor( 4);
	gr2->SetLineColor( 4);
	gr2->Draw( "Psame");

	// Theory DMT2001 at 280 MeV
	name = "xs/pi0/theory/asym_dmt_280MeV.dat";
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
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr3 = new TGraph( i-1, theta, as);
	gr3->SetLineWidth( 4);
	gr3->SetLineColor( 6);
	gr3->Draw( "Lsame");

	pt = new TLegend( 0.4, 0.30, 0.6, 0.45);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "CB-TAPS Preliminary", "p");
	pt->AddEntry( gr2, "Leukel - TAPS", "p");
	pt->AddEntry( gr3, "DMT 2001", "l");
	pt->Draw();

	name = "plots/Pi0/AsymCompDelta";
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

//
// AsymCompTAPS
//
// Compares the full AND subt asymmetries to previous TAPS measurements in
// either the threshold region (bin 0) or the Delta region (bin 1).
//
void AsymCompTAPS( UInt_t bin, Bool_t save = kFALSE)
{
	UInt_t i, chan_lo, chan_hi, egam, imax;
	Double_t eg, deg;
	Double_t theta[100], dtheta[100], as[100], das[100];
	Double_t junk;
	TString name, tgt;

	TH1D* hscal;

	if ( bin == 0) {
		chan_lo = 299;
		chan_hi = 308;
		egam = 160;
		imax = 19;
	}
	else if ( bin == 1) {
		chan_lo = 251;
		chan_hi = 255;
		egam = 280;
		imax = 17;
	}

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "AsymComp", 400, 10, 700, 500);

	hscal = (TH1D*) fparaFile.Get( "SumScalers152to503");
	eg = WeightedPar( hscal, "energy", chan_lo, chan_hi);
	deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	tgt = "full";
	name = Form( "xs/pi0/asym/asymmetry_%s_%d-%d.out", tgt.Data(), chan_lo,
			chan_hi);
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
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
		i++;
	}
	inFile.close();

	gr1 = new TGraphErrors( i-1, theta, as, dtheta, das);
	name = Form( "Asymmetry for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	gr1->SetTitle( name);
	gr1->SetMarkerStyle( 21);
	gr1->SetMarkerSize( 1.2);
	gr1->SetMarkerColor( 4);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 4);
	gr1->GetXaxis()->SetTitleOffset( 1.1);
	gr1->GetYaxis()->SetTitleOffset( 0.8);
	gr1->GetYaxis()->SetTitleSize( 0.05);
	gr1->GetXaxis()->SetTitle( "#theta_{cm} (deg)");
	gr1->GetYaxis()->SetTitle( "#Sigma");
	gr1->GetXaxis()->SetLabelSize( 0.03);
	gr1->GetYaxis()->SetLabelSize( 0.03);
	gr1->GetXaxis()->SetRangeUser( 0, 180);
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	if ( bin == 0) {
		gr1->SetMaximum( 0.3);
		gr1->SetMinimum( -0.2);
	}
	gr1->Draw( "AP");

	l1 = new TLine( 0, 0, 180, 0);
	l1->Draw();

	tgt = "subt";
	name = Form( "xs/pi0/asym/asymmetry_%s_%d-%d.out", tgt.Data(), chan_lo,
			chan_hi);
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
		inFile >> theta[i] >> as[i] >> das[i] >> junk >> junk >> junk;
		dtheta[i] = 0;
		i++;
	}
	inFile.close();
	gr2 = new TGraphErrors( i-1, theta, as, dtheta, das);
	gr2->SetMarkerStyle( 22);
	gr2->SetMarkerSize( 1.2);
	gr2->SetLineWidth( 2);
	gr2->SetMarkerColor( 1);
	gr2->SetLineColor( 1);
	gr2->Draw( "Psame");

	name = Form( "xs/pi0/previous/asym_%d.dat", egam);
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

	gr3 = new TGraphErrors( i-1, theta, as, dtheta, das);
	gr3->SetMarkerStyle( 20);
	gr3->SetMarkerSize( 1.2);
	gr3->SetLineWidth( 2);
	gr3->SetMarkerColor( 2);
	gr3->SetLineColor( 2);
	gr3->Draw( "Psame");

	pt = new TLegend( 0.25, 0.15, 0.75, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr1, "CB-TAPS Preliminary - Full Target", "p");
	pt->AddEntry( gr2, "CB-TAPS Preliminary - Subtracted", "p");
	pt->AddEntry( gr3, "TAPS", "p");
	pt->Draw();

	name = Form( "plots/Pi0/AsymComp%d", egam);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

//
// EmptyAsymCheck( Th_cm, Eg, Tgt)
//
// Check of p_gamma*Sigma for empty target data.
//
void EmptyAsymCheck( UInt_t chan_lo, UInt_t chan_hi, UInt_t tbin)
{
	UInt_t theta;
	UInt_t tb_lo, tb_hi;
	Double_t f_dead, factor;
	TFile* file;

	TH1D* asym1;
	TH1D* asym2;
	TH1D* asym3;
	TH1D* asym4;
	TH1D* asym5;

	TH2D* h2_e_perp;
	TH1D* hsc_e_perp;
	TH2D* h2_e_para;
	TH1D* hsc_e_para;

	TH2D* diff2d;
	TH2D* sum2d;

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

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 1000, 700);
	c1->Divide( 3, 2);
	c1->cd( 4);

	theta = 20*tbin - 10;

	// Empty Perp
	c1->cd( 1);
	file = eperpFile;
	GetHistFile( file);
	f_dead = DeadTimeSF( hsc);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	h2_e_perp = (TH2D*) hS->Project3D( "yz");
	hsc_e_perp = (TH1D*) hsc->Clone( "sc_e_perp");
	hsc_e_perp->Scale( f_dead);
	asym1 = (TH1D*) h2_e_perp->ProjectionX( "asym1", tb_lo, tb_hi);
	asym1->Draw();

	// Empty Para
	c1->cd( 2);
	file = eparaFile;
	GetHistFile( file);
	f_dead = DeadTimeSF( hsc);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	h2_e_para = (TH2D*) hS->Project3D( "yz1");
	hsc_e_para = (TH1D*) hsc->Clone( "sc_e_para");
	hsc_e_para->Scale( f_dead);
	asym2 = (TH1D*) h2_e_para->ProjectionX( "asym2", tb_lo, tb_hi);
	asym2->Draw();

	factor = HistSF( hsc_e_para, hsc_e_perp, chan_lo, chan_hi);
	cout << "factor = " << factor << endl;

	c1->cd( 3);
	diff2d = (TH2D*) h2_e_para->Clone( "diff2d");
	diff2d->Add( h2_e_perp, -factor);
	asym3 = (TH1D*) diff2d->ProjectionX( "asym3", tb_lo, tb_hi);
	asym3->Draw();

	c1->cd( 4);
	sum2d = (TH2D*) h2_e_para->Clone( "sum2d");
	sum2d->Add( h2_e_perp, factor);
	asym4 = (TH1D*) sum2d->ProjectionX( "asym4", tb_lo, tb_hi);
	asym4->Draw();

	c1->cd( 5);
	TH2D* asym2d;
	asym2d = (TH2D*) diff2d->Clone( "asym2d");
	asym2d->Divide( sum2d);
	asym5 = (TH1D*) asym2d->ProjectionX( "asym5", tb_lo, tb_hi);
	asym5->Draw();

	c1->cd( 6);
	TH1D* asym1d;
	asym1d = (TH1D*) asym3->Clone( "asym1d");
	asym1d->Divide( asym4);
	asym1d->Draw();

	Double_t x[6];

	// Para
	x[0] = asym2->GetBinContent( 3);
	// Perp
	x[1] = asym1->GetBinContent( 3);
	// Diff
	x[2] = x[0] - x[1]*factor;
	// Sum
	x[3] = x[0] + x[1]*factor;
	// Asym
	x[4] = x[2]/x[3];

	cout << x[0];
	cout << " " << x[1];
	cout << " " << x[2];
	cout << " " << x[3];
	cout << " " << x[4];
	cout << endl;
}

//
// PerpPara
//
void PlotPerpPara( UInt_t chan_lo, UInt_t chan_hi)
{
	UInt_t tb_lo, tb_hi;
	UInt_t fpe, fpa;
	Double_t f_dead[4], polgrad[4], te_corr[4], factor;
	Double_t fact[3], factor;
	TString name;

	tb_lo = 1;
	tb_hi = 9;

	// Tags for full perp, full para
	fpe = 0;
	fpa = 1;

	TFile* file;

	TH1D* h1_f_perp;
	TH1D* hsc_f_perp;
	TH1D* h1_f_para;
	TH1D* hsc_f_para;

	// Full Perp
	file = fperpFile;
	GetHistFile( file);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_f_perp = (TH1D*) hS->Project3D( "z");
	hsc_f_perp = (TH1D*) hsc->Clone( "sc_f_perp");
	f_dead[fpe] = DeadTimeSF( hsc_f_perp);
	polgrad[fpe] = WeightedPol( hsc_f_perp, chan_lo, chan_hi, fpe);

	// Full Para
	file = fparaFile;
	GetHistFile( file);
	hS->GetXaxis()->SetRange( chan_lo+1, chan_hi+1);
	hS->GetYaxis()->SetRange( tb_lo, tb_hi);
	h1_f_para = (TH1D*) hS->Project3D( "z3");
	hsc_f_para = (TH1D*) hsc->Clone( "sc_f_para");
	f_dead[fpa] = DeadTimeSF( hsc_f_para);
	polgrad[fpa] = WeightedPol( hsc_f_para, chan_lo, chan_hi, fpa);

	factor = HistSF( hsc_f_para, hsc_f_perp, chan_lo, chan_hi);
	fact[0] = f_dead[fpa]/f_dead[fpe];
	fact[1] = polgrad[fpa]/polgrad[fpe];
	factor *= fact[0]; 
	factor /= fact[1]; 

	Double_t eg = WeightedPar( hsc_f_perp, "energy", chan_lo, chan_hi);
	Double_t deg = ((tcd[chan_lo].energy + tcd[chan_lo].denergy) -
		(tcd[chan_hi].energy - tcd[chan_hi].denergy))/2;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "Data", 400, 10, 700, 500);
	h1_f_perp->Scale( factor);
	name = Form( "#phi-Yield for E_{#gamma} = %5.1f #pm %4.1f MeV", eg, deg);
	h1_f_perp->SetTitle( name);
	h1_f_perp->Draw();
	h1_f_perp->GetXaxis()->SetTitle( "#phi (deg)");
	h1_f_perp->GetXaxis()->CenterTitle();
	h1_f_perp->GetYaxis()->SetTitle( "Counts");
	h1_f_perp->GetYaxis()->CenterTitle();
	h1_f_perp->SetLineWidth( 4);
	h1_f_perp->SetLineColor( 2);
	h1_f_perp->SetMarkerStyle( 20);
	h1_f_perp->SetMarkerColor( 2);
	h1_f_perp->SetMarkerSize( 2);
	h1_f_perp->SetMinimum( 0);
	h1_f_perp->SetMaximum( 6000);

	h1_f_para->Draw("same");
	h1_f_para->SetLineWidth( 4);
	h1_f_para->SetLineColor( 4);
	h1_f_para->SetMarkerStyle( 21);
	h1_f_para->SetMarkerColor( 4);
	h1_f_para->SetMarkerSize( 2);

	pt = new TLegend( 0.4, 0.15, 0.6, 0.35);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( h1_f_perp, "Perpendicular", "p");
	pt->AddEntry( h1_f_para, "Parallel", "p");
	pt->Draw();

	name = "plots/Pi0/PlotPerpPara";
	name.Append( ".pdf");
	c1->Print( name);
}

//
// PlotCompCTH
//
// Compares the asymmetries from Sergey Prakhov's analysis with mine
// for his bins.
//
void PlotCompCTH( UInt_t start = 0, Bool_t save == kFALSE)
{
	UInt_t i, j, bin, cthbins;
	Double_t eg, deg;
	Double_t cth[20], dcth[20], theta[20], as[20], das[20];
	Double_t junk;
	Double_t x[2], y[2];
	TString name;

	gStyle->SetOptStat( 0);
	c1 = new TCanvas( "c1", "PlotCompCTH", 400, 10, 1000, 500);
	c1->Divide( 3, 2);

	j = 1;
	for ( UInt_t bin = start; bin < start + 6; bin++) {

		if ( bin < 4) cthbins = 7;
		else cthbins = 9;

		c1->cd( j++);

		if ( bin < 12) {
			x[0] = 0.55;
			y[0] = 0.70;
		}
		else {
			x[0] = 0.40;
			y[0] = 0.35;
		}
		x[1] = x[0] + 0.25;
		y[1] = y[0] + 0.13;

		// Get my results
		AsymmetryCTH( bin, kFALSE);
		for ( i = 0; i < cthbins; i++) {
			cth[i] = Asym[i].cth;
			dcth[i] = Asym[i].dcth;
			as[i] = Asym[i].as;
			das[i] = Asym[i].das;
		}

		// Plot my results
		gr = new TGraphErrors( cthbins, cth, as, dcth, das);
		name = Form( "E_{#gamma} = %5.1f #pm %4.1f MeV", spbin[bin].eg,
				spbin[bin].deg);
		gr->SetTitle( name);
		gr->SetMarkerStyle( 21);
		gr->SetMarkerSize( 1.5);
		gr->SetLineWidth( 3);
		gr->GetXaxis()->SetTitleOffset( 1.1);
		gr->GetYaxis()->SetTitleOffset( 0.80);
		gr->GetYaxis()->SetTitleSize( 0.06);
		gr->GetXaxis()->SetTitle( "cos #theta");
		gr->GetYaxis()->SetTitle( "#Sigma");
		gr->GetXaxis()->SetLabelSize( 0.03);
		gr->GetYaxis()->SetLabelSize( 0.03);
		gr->GetXaxis()->SetRangeUser( -1, 1);
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->CenterTitle();
		gr->SetMaximum( 0.35);
		gr->SetMinimum( -0.15);
		gr->Draw( "AP");

		// Read in Sergey's data
		name = Form( "xs/pi0/asym_sergey/asym_%5.1fMeV.dat", spbin[bin].eg);
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
			cth[i] = cos( theta[i]*kD2R);
			dcth[i] = 0;
			i++;
		}
		inFile.close();

		// Plot Sergey's data
		gr2 = new TGraphErrors( i-1, cth, as, dcth, das);
		gr2->SetMarkerStyle( 20);
		gr2->SetMarkerColor( 2);
		gr2->SetMarkerSize( 1.5);
		gr2->SetLineColor( 2);
		gr2->SetLineWidth( 3);
		gr2->Draw( "Psame");

		l1 = new TLine( -1, 0, 1, 0);
		l1->Draw();
	
		pt = new TLegend( x[0], y[0], x[1], y[1]);
		pt->SetTextSize( 0.04);
		pt->SetFillColor( 0);
		pt->SetBorderSize( 0);
		pt->AddEntry( gr, "Hornidge", "p");
		pt->AddEntry( gr2, "Prakhov", "p");
		pt->Draw();
	}
	c1->Update();

	name = Form( "plots/Pi0/AsymCompSergey_%d-%d", start, start+5);
	name.Append( ".pdf");
	if ( save == kTRUE) c1->Print( name);
}
