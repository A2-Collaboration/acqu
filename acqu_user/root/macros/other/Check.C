#include "includes/physics.h"
#include "includes/functions.h"

void CheckAll( TString subt = "full")
{
	TString name;

	DeadTimeCheck( subt);
	YieldCheck( subt);
	PolCheck( subt);

	TA2ThreshPi0* tpi0 = (TA2ThreshPi0*)gAN->GetPhysics();
	tpi0->ZeroAll();
	gAN->ZeroSumScalers();

	name = Form( "Zeroing Sum Scalers and Physics histograms...");
	cout << name << endl;
}

void YieldCheck( TString subt = "full")
{
	Int_t xb_lo, xb_hi, yb_lo, yb_hi;
	Int_t tchan_lo, tchan_hi;
	Int_t run;
	Double_t mmiss_lo, mmiss_hi;
	Double_t pa, yield, scalers;
	Double_t dyield;
	Double_t y_sc, dy_sc;
	Double_t f_tagg, f_tot, f_dead;
	Double_t norm;
	TString name;

	pa = 0.0833;

	// For Full target
	tchan_lo = 228;
	tchan_hi = 245;
	mmiss_lo = 890;
	mmiss_hi = 970;
	norm = 1E5;

	// For Empty target
	if ( subt == "empty") {
		tchan_lo = 290;
		tchan_hi = 312;
		mmiss_lo = 910;
		mmiss_hi = 950;
		norm = 1E7;
	}

	TH1D hsc;
	TH1D hdead;
	TH2D *hp2d;
	TH2D *hr2d;

	// Dead Time Scalers
	name = "SumScalers504to535";
	hdead = (TH1D)gROOT->Get( name);

	f_tagg = hdead.GetBinContent( 32)/hdead.GetBinContent( 31);
	f_tot = hdead.GetBinContent( 25)/hdead.GetBinContent( 26);
	f_dead = f_tagg/f_tot;

	// Tagger Scalers
	name = "SumScalers152to503";
	hsc = (TH1D)gROOT->Get( name);

	xb_lo = hsc.GetXaxis()->FindBin( tchan_lo);
	xb_hi = hsc.GetXaxis()->FindBin( tchan_hi);
	scalers = hsc.Integral( xb_lo, xb_hi);
	scalers /= f_dead;

	// Prompt
	name = "THR_MMissP_v_TChanP";
	hp2d = (TH2D)gROOT->Get( name);

	// Random
	name = "THR_MMissR_v_TChanR";
	hr2d = (TH2D)gROOT->Get( name);

	// Subtracted
	hp2d->Sumw2();
	hp2d->Add( hr2d, -pa);

	xb_lo = hp2d->GetXaxis()->FindBin( tchan_lo);
	xb_hi = hp2d->GetXaxis()->FindBin( tchan_hi);
	yb_lo = hp2d->GetYaxis()->FindBin( mmiss_lo);
	yb_hi = hp2d->GetYaxis()->FindBin( mmiss_hi);

	dyield = 0;
	yield = hp2d->IntegralAndError( xb_lo, xb_hi, yb_lo, yb_hi, dyield);

	y_sc = yield/scalers*norm;
	dy_sc = dyield/scalers*norm;

/*
	hp2d->GetXaxis()->SetRange( tchan_lo, tchan_hi);
	TH1D *hp = hp2d->ProjectionY( "MMiss");

	xb_lo = hp->GetXaxis()->FindBin( mmiss_lo);
	xb_hi = hp->GetXaxis()->FindBin( mmiss_hi);
	yield = hp->Integral( xb_lo, xb_hi);
	cout << yield << endl;

	hp->Draw( "E");
*/

	run = GetRunDave();

	name = "xs/yield_chk_";
	name.Append( subt);
	name.Append( ".out");

	ofstream outFile( name, ios::app);
	name = Form( "%5d %7.2e %7.2e %7.2e %5.3f %5.3f", run, yield, dyield,
			scalers, y_sc, dy_sc);
	outFile << name << endl;

	name = Form( "Run %d", run);
	cout << name << endl;

	name = Form( "Target %s", subt.Data());
	cout << name << endl;

	name = Form( "Corrected Scalers: %7.2e", scalers);
	cout << name << endl;

	name = Form( "Yield %7.2e +/- %7.2e", yield, dyield);
	cout << name << endl;

	name = Form( "Yield/Scalers*Factor is %5.3f +/- %5.3f", y_sc, dy_sc);
	cout << name << endl;

}

void DeadTimeCheck( TString subt = "full")
{
	Int_t run;
	Double_t f_tagg, f_tot, f_dead;
	TString name;

	TH1D hdead;

	// Dead Time Scalers
	name = "SumScalers504to535";
	hdead = (TH1D)gROOT->Get( name);

	f_tagg = hdead.GetBinContent( 32)/hdead.GetBinContent( 31);
	f_tot = hdead.GetBinContent( 25)/hdead.GetBinContent( 26);
	f_dead = f_tagg/f_tot;

	run = GetRunDave();

	name = "xs/deadtime_chk_";
	name.Append( subt);
	name.Append( ".out");

	ofstream outFile( name, ios::app);
	name = Form( "%5d %6.4f %6.4f %6.4f", run, f_tagg, f_tot, f_dead);
	outFile << name << endl;

	name = Form( "Run %d", run);
	cout << name;

	name = Form( " Target %s", subt.Data());
	cout << name;

	name = Form( "  Tagg / Tot / Dead = %6.4f / %6.4f / %6.4f", f_tagg, f_tot,
			f_dead);
	cout << name << endl;

//	gAN->ZeroSumScalers();
}

void PolCheck( TString subt = "full")
{
	Int_t tchan_lo, tchan_hi, run;
	Int_t xb_lo, xb_hi;
	Double_t f_tagg, f_tot, f_dead, scalers, ref_scalers;
	Double_t pa, factor;
	Double_t par[4], red_chisq;
	TString name;

//	TCanvas *c1 = new TCanvas ( "c1", "Pol Check", 200, 20, 700, 1000);
//	c1->Divide(1,3);

	pa = 0.0833;

	// For Full target
	tchan_lo = 228;
	tchan_hi = 245;

	// For Empty target
	if ( subt == "empty") {
		tchan_lo = 290;
		tchan_hi = 312;
	}

// Current Analysis Run

	// Dead Time Scalers
	name = "SumScalers504to535";
	TH1D *hdead = (TH1D*)gROOT->Get( name);

	// Tagger Scalers
	name = "SumScalers152to503";
	TH1D *htag = (TH1D*)gROOT->Get( name);

	f_tagg = hdead->GetBinContent( 32)/hdead->GetBinContent( 31);
	f_tot = hdead->GetBinContent( 25)/hdead->GetBinContent( 26);
	f_dead = f_tagg/f_tot;

	xb_lo = htag->GetXaxis()->FindBin( tchan_lo);
	xb_hi = htag->GetXaxis()->FindBin( tchan_hi);
	scalers = htag->Integral( xb_lo, xb_hi);
	scalers /= f_dead;

	// Prompt
	name = "THR_PhiCMP_v_TChanP";
	TH2D *hp2d = (TH2D*)gROOT->Get( name);

	// Random
	name = "THR_PhiCMR_v_TChanR";
	TH2D *hr2d = (TH2D*)gROOT->Get( name);

	// Subtracted
	hp2d->Sumw2();
	hp2d->Add( hr2d, -pa);

	xb_lo = hp2d->GetXaxis()->FindBin( tchan_lo);
	xb_hi = hp2d->GetXaxis()->FindBin( tchan_hi);
	hp2d->GetXaxis()->SetRange( xb_lo, xb_hi);
	TH1D *hp = hp2d->ProjectionY( "phiCM");

//	c1->cd(1);
//	hp->Draw("E");

/*
	cout << f_tagg;
	cout << " " << f_tot;
	cout << " " << f_dead;
	cout << " " << scalers;
	cout << endl;
*/

// Reference Run

	run = 20858;									// Full
	if ( subt == "empty") run = 20966;		// Empty

	name = Form( "data/refhist_%s_%d.root", subt.Data(), run);
	TFile refrun( name);

	name = Form( "Reference Run %d", run);
	cout << name << endl;

	// Dead Time Scalers
	name = "SumScalers504to535";
	TH1D *hdead_ref = (TH1D*)refrun.Get( name);

	// Tagger Scalers
	name = "SumScalers152to503";
	TH1D *htag_ref = (TH1D*)refrun.Get( name);

	f_tagg = hdead_ref->GetBinContent( 32)/hdead_ref->GetBinContent( 31);
	f_tot = hdead_ref->GetBinContent( 25)/hdead_ref->GetBinContent( 26);
	f_dead = f_tagg/f_tot;

	xb_lo = htag_ref->GetXaxis()->FindBin( tchan_lo);
	xb_hi = htag_ref->GetXaxis()->FindBin( tchan_hi);
	ref_scalers = htag_ref->Integral( xb_lo, xb_hi);
	ref_scalers /= f_dead;

/*
	cout << f_tagg;
	cout << " " << f_tot;
	cout << " " << f_dead;
	cout << " " << scalers;
	cout << endl;
*/

	factor = scalers/ref_scalers;

//	cout << factor;
//	cout << endl;

	// Prompt
	name = "THR_PhiCMP_v_TChanP";
	TH2D *hp2d_ref = (TH2D*)refrun.Get( name);

	// Random
	name = "THR_PhiCMR_v_TChanR";
	TH2D *hr2d_ref = (TH2D*)refrun.Get( name);

	// Subtracted
	hp2d_ref->Sumw2();
	hp2d_ref->Add( hr2d_ref, -pa);

	xb_lo = hp2d_ref->GetXaxis()->FindBin( tchan_lo);
	xb_hi = hp2d_ref->GetXaxis()->FindBin( tchan_hi);
	hp2d_ref->GetXaxis()->SetRange( xb_lo, xb_hi);
	TH1D *hp_ref = hp2d_ref->ProjectionY( "phiCMref");

	hp_ref->Scale( factor);

//	c1->cd(2);
//	hp_ref->Draw("E");

	TH1D *sum = (TH1D*)hp->Clone();
	TH1D *diff = (TH1D*)hp->Clone();

	sum->Add( hp_ref); 
	diff->Add( hp_ref, -1);
	diff->Divide( sum);

//	c1->cd(3);
//	diff->Draw( "E");

	TF1 *f1 = new TF1( "f1", "[0]*cos(2*x*0.01745+[1])", -180, 180);
	diff->Fit( "f1", "REMQ0");
//	diff->Fit( "f1", "REMQ");
	par[0] = f1->GetParameter(0);
	par[1] = f1->GetParError(0);
	par[2] = f1->GetParameter(1);
	par[3] = f1->GetParError(1);
	red_chisq = f1->GetChisquare()/17;
	if ( par[2] < 0) {
		par[0] *= -1;
		par[2] *= -1;
	}
/*
	cout << par[0];
	cout << " +/- " << par[1];
	cout << "   " << par[2];
	cout << " +/- " << par[3];
	cout << "   " << red_chisq;
	cout << endl;
*/

	run = GetRunDave();
	name = "xs/pol_chk_";
	name.Append( subt);
	name.Append( ".out");

	ofstream outFile( name, ios::app);
	name = Form( "%5d %6.4f %6.4f %6.4f", run, par[0], par[1], red_chisq);
	outFile << name << endl;

	name = Form( "Run %d", run);
	cout << name << endl;

	name = Form( "Target %s", subt.Data());
	cout << name << endl;

	name = Form( "p_gamma*Sigma = %6.4f +/- %6.4f", par[0], par[1]);
	cout << name << endl;

	name = Form( "phi = %6.4f +/- %6.4f", par[2], par[3]);
	cout << name << endl;

	name = Form( "Reduced Chi_sq = %6.4f", red_chisq);
	cout << name << endl;

	delete hp;
	delete hp_ref;
	delete diff;
	delete sum;
	delete f1;

//	c1->Print("plots/polcheck.pdf");
}

// Returns the file number for the current run.  Note that it gives you the
// NEXT run if you use it with EndFile!!
UInt_t GetRunDave()
{
	Int_t i1, i2, run;
	TString name;

	if ( gAR->IsOnline())
	{
		name = gAR->GetFileName();

		i1 = name.Index("_");
		name.Remove(0,i1+1);
		i2 = name.Index(".");
		name.Remove( i2);
		run = name.Atoi();
	}
	else {

		Int_t n = 0;
		while ( gAR->GetTreeFileList( n) != NULL)
			name = gAR->GetTreeFileList( n++);

		i1 = name.Index("_");
		name.Remove(0,i1+1);
		i2 = name.Index(".");
		name.Remove( i2);
		run = name.Atoi();
	}

	return( run);
}
