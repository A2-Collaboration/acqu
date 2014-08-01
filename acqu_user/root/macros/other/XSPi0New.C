// Macro XSPi0New.C
//
// NOTE:  The tagger channels run from 0 to 351.
//
// DLH		2010.06.11	New version to check Alex's cross section numbers.
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
	Double_t edet[10];
	Double_t d_edet[10];
	Double_t pa_ratio;
	Double_t xs[10];
	Double_t dxs[10];

} TData;
TData tcd[352];

typedef struct {
	Int_t th;
	Double_t theta;
	Double_t dtheta;
	Double_t dom;
} TBins;
TBins tbin[10];

Double_t f_dead_f, f_dead_e;
Int_t broken_ladd_ch[] = { -1};

TString prompt, random, scalers, dtscalers;
prompt = "THR_MMissP_v_TChanP";
random = "THR_MMissR_v_TChanR";
scalers = "SumScalers152to503";
dtscalers = "SumScalers504to535";

// Target Full
TFile full( "histograms/Pi0/Full/ycheck.root");
TH2D *hf2dp = (TH2D*)full.Get( prompt);
TH2D *hf2dr = (TH2D*)full.Get( random);
TH1D *hfsc = (TH1D*)full.Get( scalers);
TH1D *hdeadf = (TH1D*)full.Get( dtscalers);

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

	ReadDetEff( "pi0");
	ReadTagEff( "xs/tag_eff/TaggEff.out");

	// Full target scaler deadtime correction.
	f_tagg = hdeadf->GetBinContent( 32)/hdeadf->GetBinContent( 31);
	f_tot = hdeadf->GetBinContent( 25)/hdeadf->GetBinContent( 26);
	f_dead_f = f_tagg/f_tot;
	cout << "f_tagg = " << f_tagg;
	cout << "  f_tot = " << f_tot << endl;
	cout << "full tagger scalers deadtime correction = " << f_dead_f << endl;

	// Zero scalers for broken ladder channels.
	Int_t* lch = broken_ladd_ch;
	while ( *lch != -1) {
		hfsc->SetBinContent( *lch, 0);
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

void ChanXS( UInt_t i) 
{
	UInt_t xmin, xmax, ymin, ymax, j;
	Double_t t, fact, n_yld, dn_yld;
	Double_t yield_p, yield_r;
	Double_t fyield, dfyield, fscalers, dfscalers;
	Double_t yield, dyield;

	j = 0;		// Total cross section, theta bin zero.

	t = 4.242e23;		// This is in nuclei/cm^2 for the 10-cm cell.

	// The 1e30 converts cm^2 to microbarn.
	fact = 1e30/t;

	// Limits over which to integrate 3D histogram:
	//		- One tagger channel.
	// 	- MMiss region.
	//  NOTE: Bins start at 1 and not 0!!!!

	// Tagger Channel (e.g., channel 0 is bin 1)
	xmin = i+1;
	xmax = i+1;

	// Missing Mass from 900-1000 MeV.
	ymin = hf2dp->GetYaxis()->FindBin( 900);
	ymax = hf2dp->GetYaxis()->FindBin( 1000);

	// Full target
	yield_p = hf2dp->Integral( xmin, xmax, ymin, ymax);
	yield_r = hf2dr->Integral( xmin, xmax, ymin, ymax);
	fyield = yield_p - tcd[i].pa_ratio*yield_r;
	dfyield = sqrt( yield_p - Sqr( tcd[i].pa_ratio)*yield_r);
	fscalers = hfsc->GetBinContent( xmin)/f_dead_f;
	dfscalers = sqrt( fscalers);

	if ( fscalers != 0) {

		n_yld = fyield/fscalers;
		dn_yld = 1/fscalers*sqrt( Sqr( dfyield)
				+ Sqr( fyield*dfscalers/fscalers));
		yield = fyield;
		dyield = dfyield;

		tcd[i].xs[j] = n_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;
		tcd[i].dxs[j] = dn_yld*fact/tcd[i].edet[j]/tcd[i].etag/tbin[j].dom;

		cout << "Chan = " << i;
		cout << "  " << yield_p;
		cout << "  " << yield_r;
		cout << "  " << tcd[i].pa_ratio;
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
