gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

TFile p2( "histograms/MonteCarlo/pi0_p_150_2cm.root");
TFile w2( "histograms/MonteCarlo/pi0_w_150_2cm.root");

TFile p5( "histograms/MonteCarlo/pi0_p_150_5cm.root");
TFile w5( "histograms/MonteCarlo/pi0_w_150_5cm.root");

TFile p10( "histograms/MonteCarlo/pi0_p_150_10cm.root");
TFile w10( "histograms/MonteCarlo/pi0_w_150_10cm.root");

void TgtLen( Bool_t save = kFALSE)
{
	Double_t left, right;
	TString name;

	TCanvas *c1 = new TCanvas ( "c1", "Target Length", 200, 20, 1000, 700);
	c1->Divide( 3, 2);

	gStyle->SetOptStat( 0);

	name = "THR_M2g";

	TH1D *m2g_p2 = (TH1D*)p2.Get( name);
	TH1D *m2g_w2 = (TH1D*)w2.Get( name);

	TH1D *m2g_p5 = (TH1D*)p5.Get( name);
	TH1D *m2g_w5 = (TH1D*)w5.Get( name);

	TH1D *m2g_p10 = (TH1D*)p10.Get( name);
	TH1D *m2g_w10 = (TH1D*)w10.Get( name);

	left = 0;
	right = 200;

	name = "#pi^{0} Invariant Mass E_{#gamma} = 150 MeV";

	c1->cd( 1);
	m2g_p2->SetTitle( name);
	m2g_p2->SetLineWidth( 3);
	m2g_p2->SetLineColor( 1);
	m2g_p2->GetXaxis()->SetRangeUser( left, right);
	m2g_p2->Draw();

	m2g_p5->SetLineWidth( 3);
	m2g_p5->SetLineColor( 2);
	m2g_p5->Draw( "same");

	m2g_p10->SetLineWidth( 3);
	m2g_p10->SetLineColor( 4);
	m2g_p10->Draw( "same");

	pt = new TLegend( 0.15, 0.65, 0.35, 0.85);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( m2g_p2, "2 cm", "l");
	pt->AddEntry( m2g_p5, "5 cm", "l");
	pt->AddEntry( m2g_p10, "10 cm", "l");
	pt->Draw();

	name = "THR_MMissPi0P";

	TH1D *mmiss_p2 = (TH1D*)p2.Get( name);
	TH1D *mmiss_w2 = (TH1D*)w2.Get( name);

	TH1D *mmiss_p5 = (TH1D*)p5.Get( name);
	TH1D *mmiss_w5 = (TH1D*)w5.Get( name);

	TH1D *mmiss_p10 = (TH1D*)p10.Get( name);
	TH1D *mmiss_w10 = (TH1D*)w10.Get( name);

	left = 900;
	right = 960;

	c1->cd(2);

	name = "Proton Missing Mass E_{#gamma} = 150 MeV";
	mmiss_p2->SetTitle( name);
	mmiss_p2->SetLineWidth( 3);
	mmiss_p2->SetLineColor( 1);
	mmiss_p2->GetXaxis()->SetRangeUser( left, right);
	mmiss_p2->Draw();
	mmiss_p2->Fit( "gaus");

	mmiss_p5->SetLineWidth( 3);
	mmiss_p5->SetLineColor( 2);
	mmiss_p5->Draw( "same");
	mmiss_p5->Fit( "gaus");

	mmiss_p10->SetLineWidth( 3);
	mmiss_p10->SetLineColor( 4);
	mmiss_p10->Draw( "same");
	mmiss_p10->Fit( "gaus");

	pt = new TLegend( 0.15, 0.65, 0.35, 0.85);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( mmiss_p2, "2 cm", "l");
	pt->AddEntry( mmiss_p5, "5 cm", "l");
	pt->AddEntry( mmiss_p10, "10 cm", "l");
	pt->Draw();

	name = Form( "plots/TgtLen.pdf");
	if ( save == kTRUE) c1->Print( name);
}
