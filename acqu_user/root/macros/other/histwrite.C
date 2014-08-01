gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void HistWrite( TString subt = "full")
{
	Int_t run;
	TString name;

	// Tagger Scalers
	name = "SumScalers152to503";
	TH1D *htag = (TH1D*)gROOT->Get( name);

	// DeadTime Scalers
	name = "SumScalers504to535";
	TH1D *hdead = (TH1D*)gROOT->Get( name);

	// Prompt
	name = "THR_PhiCMP_v_TChanP";
	TH1D *hpro = (TH1D*)gROOT->Get( name);

	// Random
	name = "THR_PhiCMR_v_TChanR";
	TH1D *hran = (TH1D*)gROOT->Get( name);

	TObjArray hlist( 0);
	hlist.Add( hpro);
	hlist.Add( hran);
	hlist.Add( htag);
	hlist.Add( hdead);

	run = GetRunDave();
	name = Form( "data/refhist_%s_%d.root", subt.Data(), run);
	TFile f1( name, "recreate", "histcheck");
	hlist.Write();
	f1.Close();
}
