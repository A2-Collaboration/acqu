gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void SimCarb( Int_t eg, Char_t hist[32] = "TGG")
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, a, b, fact;
	Char_t name[64], cond[16];

	strcpy( cond, "masscorr");
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	carb = new TFile( name);

	TCanvas *c1 = new TCanvas ( "c1", "Carbon Simulation", 200, 0, 700, 500);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	if ( strcmp( hist, "TGG") == 0) {

		h3c = (TH3D*) carb->Get( "THits_TGG_ThetaCM_P");
		TH1D *h1 = h3c->Project3D("y2");

		sprintf( name, "#Theta^{min}_{#gamma#gamma} for E_{#gamma} = %5.1f MeV",
				eg);
		h1->SetTitle( name);
		x2 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
		col = 4;
		a = x2 - 10;
		b = 180;
	}
	else if ( strcmp( hist, "EmPi0") == 0) {

		h2c = (TH2D*) carb->Get( "THits_EmPi0_P");
		TH1D *h1 = h2c->ProjectionX( "carb");

		sprintf( name, "M_{miss} for E_{#gamma} = %5.1f MeV", eg);
		h1->SetTitle( name);
		x1 = Linear( 145, 932, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		col = 2;
		a = Linear( 145, 900, 400, 600, eg);
		b = Linear( 145, 950, 400, 1000, eg);
	}

	h1->SetLineWidth(2);
	h1->SetFillColor(2);
	h1->GetXaxis()->SetRangeUser(a,b);
	h1->Draw();

	max = h1->GetMaximum();

	max *= 1.10;
	h1->SetMaximum(max);

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineColor(2);
	l1->SetLineStyle(2);
	l1->SetLineWidth(4);
	l1->Draw();

	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineColor(col);
	l2->SetLineStyle(2);
	l2->SetLineWidth(4);
	l2->Draw();

//	sprintf( name, "plots/eps/Sim_%s_%s_%d.eps", hist, cond, eg);
//	c1->Print( name);
}
