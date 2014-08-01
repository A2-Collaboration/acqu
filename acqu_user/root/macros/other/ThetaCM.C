gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void ThetaCM( UInt_t chan)
{
	TString name;

	TCanvas* canv;
	if ( !( canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		canv = new TCanvas ( "c1", "Subtraction", 300, 20, 700, 400);

   canv->SetGrid();
   canv->GetFrame()->SetFillColor( 21);
   canv->GetFrame()->SetBorderSize( 12);

	name = Form( "histograms/MonteCarlo/pi0/pi0_p_chan%d.root", chan);
	sim = new TFile( name);

	name = "THR_ThetaCMP";
	tcm = (TH1D*) sim->Get( name);
	tcm->SetMaximum( 1200);
	tcm->Draw();

	name = "THR_ThetaCMCut1P";
	tcm1 = (TH1D*) sim->Get( name);
	tcm1->SetLineStyle( 2);
	tcm1->Draw( "same");

	name = "THR_ThetaCMCut2P";
	tcm2 = (TH1D*) sim->Get( name);
	tcm2->SetLineStyle( 3);
	tcm2->Draw( "same");
}
