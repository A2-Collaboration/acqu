// Stuff to do at the end of an analysis run.

#include "includes/physics.h"
#include "includes/functions.h"

void FinishBasic( TString name)
{
	TString path;

	cout << "End-of-Run macro executing";
	cout << endl;

	path.Append( name);
	path.Append( ".root");
	TFile f( path, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	path.Prepend( "All histograms saved to \"");
	path.Append( "\"");
	cout << path;
	cout << endl;
}

void FinishData( TString exp, TString tgt, TString name)
{
	TString path;

	cout << "End-of-Run macro executing";
	cout << endl;

	path = "histograms/" + exp + "/" + tgt + "/";

	path.Append( name);
	path.Append( ".root");
	TFile f( path, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	path.Prepend( "All histograms saved to \"");
	path.Append( "\"");
	cout << path;
	cout << endl;
}

void FinishMC( Char_t tgt[2], Char_t cond[32], Int_t eg)
{
	Char_t dir[128];

	sprintf( dir, "histograms/MonteCarlo/%s_%s_%d.root", tgt, cond, eg);

	printf( "End-of-Run macro executing\n");

	TFile f( dir, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	printf( "All histograms saved to %s\n", dir);

	gROOT->ProcessLine( ".L ThetaEff.C");

	ThetaEffBatch( eg, tgt, cond);

	printf( "Theta efficiencies calculated for %s %s %d\n", tgt, cond, eg);
}

void FinishMC2( TString proc, Int_t eg, TString cond)
{
	TString name;

	name = Form( "histograms/MonteCarlo/%s_%d_%s.root", proc.Data(), eg,
			cond.Data());

	printf( "End-of-Run macro executing\n");

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to ");
	cout << name << endl;

//	gROOT->ProcessLine( ".L ThetaEff.C");
//	ThetaEffBatch2( tgt, eg, cond);
//	printf( "Theta efficiencies calculated for %s %d %s\n", tgt, eg, cond);
}

void FinishMC3( TString savedir, TString encl, TString process, TString tgt,
		UInt_t eg)
{
	TString name;

	name = "End-of-Run macro executing";
	cout << name << endl;

	name = Form( "histograms/MonteCarlo/%s/%s/%s_%s_%d.root",
			savedir.Data(), encl.Data(), process.Data(),
			tgt.Data(), eg);

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms save to ");
	cout << name << endl;
}

void FinishMC4( UInt_t chan, TString expt, TString tgt, TString thick,
		TString process, UInt_t effc = 0)
{
	TString dir;

	dir = Form( "histograms/MonteCarlo/%s/%s_%s/%s_chan%d.root",
			expt.Data(), tgt.Data(), thick.Data(), process.Data(), chan);

	printf( "End-of-Run macro executing\n");

	TFile f( dir, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	printf( "All histograms saved to %s\n", dir.Data());

	if ( effc == 1) {

		gROOT->ProcessLine( ".L ThetaEff.C");

		ThetaEff( chan, expt, tgt, thick, process, kFALSE);

		printf( "Theta efficiencies for experiment %s, target %s, thickness %s,"
					" process %s, and tagger channel %d\n", expt.Data(),
					tgt.Data(), thick.Data(), process.Data(), chan);
	}
	else if ( effc == 2) {

		gROOT->ProcessLine( ".L ThetaEff.C");

		CThetaEff( chan, expt, tgt, thick, process, kFALSE);

		printf( "Cos(Theta) efficiencies for experiment %s, target %s,"
					" thickness %s, process %s, and tagger channel %d\n",
					expt.Data(), tgt.Data(), thick.Data(), process.Data(), chan);
	}
	else if ( effc == 3) {

		gROOT->ProcessLine( ".L ThetaEff.C");

		ThetaEff( chan, expt, tgt, thick, process, kFALSE);
		CThetaEff( chan, expt, tgt, thick, process, kFALSE);

		printf( "Both efficiencies for experiment %s, target %s,"
					" thickness %s, process %s, and tagger channel %d\n",
					expt.Data(), tgt.Data(), thick.Data(), process.Data(), chan);
	}
}

void FinishMC5( TString process, Int_t energy, Int_t chan)
{
	TString name;

	name = Form( "histograms/MonteCarlo/%s_e%d_t%d.root",
			process.Data(), energy, chan);

	printf( "End-of-Run macro executing\n");

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to ");
	cout << name << endl;
}

void FinishTE()
{
	Int_t run;
	TString name;

	run = GetRunDave();

	name = Form( "histograms/TaggEff/TaggEff_%d.root", run);

	cout << "End-of-Run macro executing";
	cout << endl;

	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to \"");
	name.Append( "\"");
	cout << name;
	cout << endl;
}
