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

//TFile full( "histograms/ARH.root");
//TFile full( "histograms/Compton/Full/Full.root");
TFile full( "histograms/Compton/Full/Full2.root");

Int_t tchan[] = { 286, 282, 277, 273, 269, 265, 261, 257, 253, 249, 245};

void PhotTaggTime( UInt_t ebin = 8, Bool_t print = kFALSE)
{
	UInt_t i, chan_lo, chan_hi, bin_lo, bin_hi;
	UInt_t prompt[2], acc[4];
	Double_t yy, ym, y1, y2;
	Double_t energy, temin, temax;
	TString name;

	TFile* file;

	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

	// Middle channel bin and bin limits
	i = tchan[ebin];
	chan_lo = i-2;
	chan_hi = i+2;

	file = full;
	energy = WeightedPar( file, "energy", chan_lo, chan_hi);
	temin = tcd[chan_hi].energy - tcd[chan_hi].denergy;
	temax = tcd[chan_lo].energy + tcd[chan_lo].denergy;

	gStyle->SetOptStat( 0);

	prompt[0] = 92;
	prompt[1] = 106;
	acc[0] = 20;
	acc[1] = 80;
	acc[2] = 118;
	acc[3] = 178;

	TCanvas *c1 = new TCanvas ( "c1", "Canvas", 200, 0, 1200, 500);
	c1->Divide( 3,1);

	// Phot Tagger Time
	c1->cd( 1);

	name = "COMP_PhotTaggTime_v_TChanHit";
	TH2D *h2D = (TH2D*)full.Get( name);

	bin_lo = h2D->GetXaxis()->FindBin( chan_lo);
	bin_hi = h2D->GetXaxis()->FindBin( chan_hi);
	TH1D *h1 = h2D->ProjectionY( "PhotTaggTime", bin_lo, bin_hi);

	name = Form( " E_{#gamma} = %5.1f (%5.1f - %5.1f) MeV", energy, temin,
			temax);
	h1->SetTitle( name);
	h1->SetLineWidth(2);
	h1->GetXaxis()->SetTitleOffset( 1.1);
	h1->GetYaxis()->SetTitleOffset( 1.3);
	h1->GetXaxis()->SetTitle( "Photon - Tagger Time (ns)");
	h1->GetYaxis()->SetTitle( "Counts");
	h1->GetXaxis()->SetLabelSize( 0.03);
	h1->GetYaxis()->SetLabelSize( 0.03);
	h1->GetXaxis()->CenterTitle();
	h1->GetYaxis()->CenterTitle();
	h1->Draw();

	h1->SetMinimum( 0);
	yy = h1->GetMaximum();
	ym = yy*1.05;
	h1->SetMaximum( ym);
	y1 = yy*0.9;
	y2 = yy;

	l1 = new TLine( prompt[0], 0,prompt[0], ym);
	l1->SetLineWidth(3);
	l1->Draw();
	l2 = new TLine( prompt[1], 0,prompt[1], ym);
	l2->SetLineWidth(3);
	l2->Draw();
	l3 = new TLine( acc[0], 0,acc[0], ym);
	l3->SetLineStyle(2);
	l3->SetLineWidth(3);
	l3->Draw();
	l4 = new TLine( acc[1], 0,acc[1], ym);
	l4->SetLineStyle(2);
	l4->SetLineWidth(3);
	l4->Draw();
	l5 = new TLine( acc[2], 0,acc[2], ym);
	l5->SetLineStyle(2);
	l5->SetLineWidth(3);
	l5->Draw();
	l6 = new TLine( acc[3], 0,acc[3], ym);
	l6->SetLineStyle(2);
	l6->SetLineWidth(3);
	l6->Draw();

	// Photon Missing Mass
	c1->cd( 2);

	name = "COMP_PhotonMmissP_v_TChanPhotP";
	TH2D *h2P = (TH2D*)full.Get( name);
	name = "COMP_PhotonMmissR_v_TChanPhotR";
	TH2D *h2R = (TH2D*)full.Get( name);

	bin_lo = h2P->GetXaxis()->FindBin( chan_lo);
	bin_hi = h2P->GetXaxis()->FindBin( chan_hi);
	TH1D *hP = h2P->ProjectionY( "PhotonMmissP", bin_lo, bin_hi);
	TH1D *hR = h2R->ProjectionY( "PhotonMmissR", bin_lo, bin_hi);
	TH1D *hS = (TH1D*)hP->Clone( "subt");

	// TDC Spectrum
	hP->SetTitle();
	hP->SetLineWidth( 2);
	hP->GetXaxis()->SetTitleOffset( 1.1);
	hP->GetXaxis()->SetTitle( "Proton Missing Mass (MeV)");
	hP->GetXaxis()->SetLabelSize( 0.03);
	hP->GetXaxis()->CenterTitle();
	hP->GetXaxis()->SetRangeUser( 800, 1200);
	hP->Draw();

	hR->SetLineWidth( 2);
	hR->SetLineStyle( 2);
	hR->Scale( 0.0833);
	hR->Draw( "same");

	tl = new TLegend( 0.2, 0.7, 0.4, 0.8);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hP, "prompt");
	tl->AddEntry( hR, "random");
	tl->Draw();

	c1->cd( 3);
	hS->Sumw2();
	hS->Add( hR, -1.0);

	hS->SetTitle();
	hS->GetXaxis()->SetTitleOffset( 1.1);
	hS->GetXaxis()->SetTitle( "Proton Missing Mass (MeV)");
	hS->GetXaxis()->SetLabelSize( 0.03);
	hS->GetXaxis()->CenterTitle();
	hS->SetLineWidth( 2);
	hS->SetMarkerStyle( 20);
	hS->SetMarkerSize( 0.5);
	hS->GetXaxis()->SetRangeUser( 800, 1200);
	hS->Draw();

	tl = new TLegend( 0.2, 0.7, 0.4, 0.8);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hS, "subtracted", "p");
	tl->Draw();

	name = Form( "plots/Compton/random_subt_%d", ebin);
	name.Append( ".pdf");
	if ( print == kTRUE) c1->Print( name);
}

void PlotMmiss( UInt_t ebin = 8, Bool_t print = kFALSE)
{
	UInt_t i, chan_lo, chan_hi, bin_lo, bin_hi;
	Double_t x[2], y[2];
	Double_t energy, temin, temax;
	TString name;

	TFile* file;

	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

	// Middle channel bin and bin limits
	i = tchan[ebin];
	chan_lo = i-2;
	chan_hi = i+2;

	file = full;
	energy = WeightedPar( file, "energy", chan_lo, chan_hi);
	temin = tcd[chan_hi].energy - tcd[chan_hi].denergy;
	temax = tcd[chan_lo].energy + tcd[chan_lo].denergy;

	gStyle->SetOptStat( 0);

	TCanvas *c1 = new TCanvas ( "c1", "Canvas", 200, 0, 1200, 500);
	c1->Divide( 3,1);

	// No cut
	c1->cd( 1);
	name = "COMP_PhotonMmissP_v_TChanPhotP";
	TH2D *h2P = (TH2D*)full.Get( name);
	name = "COMP_PhotonMmissR_v_TChanPhotR";
	TH2D *h2R = (TH2D*)full.Get( name);
	h2P->Sumw2();
	h2P->Add( h2R, -0.0833);
	bin_lo = h2P->GetXaxis()->FindBin( chan_lo);
	bin_hi = h2P->GetXaxis()->FindBin( chan_hi);
	TH1D *hS = h2P->ProjectionY( "PhotonMmissP", bin_lo, bin_hi);

	name = Form( " E_{#gamma} = %5.1f (%5.1f - %5.1f) MeV", energy, temin,
			temax);
	hS->SetTitle( name);
	hS->SetLineWidth( 2);
	hS->GetXaxis()->SetTitleOffset( 1.1);
	hS->GetXaxis()->SetTitle( "Proton Missing Mass (MeV)");
	hS->GetXaxis()->SetLabelSize( 0.03);
	hS->GetYaxis()->SetLabelSize( 0.03);
	hS->GetXaxis()->CenterTitle();
	hS->GetXaxis()->SetRangeUser( 800, 1200);
	hS->Draw();

	x[0] = 0.2;
	x[1] = 0.4;
	y[0] = 0.88;
	y[1] = 0.96;

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hS, "no cut", "p");
	tl->Draw();

	// Require Proton
	c1->cd( 2);
	name = "COMP_PhotonMmissProtP_v_TChanPhotProtP";
	TH2D *h2P = (TH2D*)full.Get( name);
	name = "COMP_PhotonMmissProtR_v_TChanPhotProtR";
	TH2D *h2R = (TH2D*)full.Get( name);
	h2P->Sumw2();
	h2P->Add( h2R, -0.0833);
	bin_lo = h2P->GetXaxis()->FindBin( chan_lo);
	bin_hi = h2P->GetXaxis()->FindBin( chan_hi);
	TH1D *hS = h2P->ProjectionY( "PhotonMmissProtP", bin_lo, bin_hi);

	hS->SetTitle();
	hS->SetLineWidth( 2);
	hS->GetXaxis()->SetTitleOffset( 1.1);
	hS->GetXaxis()->SetTitle( "Proton Missing Mass (MeV)");
	hS->GetXaxis()->SetLabelSize( 0.03);
	hS->GetYaxis()->SetLabelSize( 0.03);
	hS->GetXaxis()->CenterTitle();
	hS->GetXaxis()->SetRangeUser( 800, 1200);
	hS->Draw();

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hS, "require proton hit", "p");
	tl->Draw();

	// Require Proton and OA Cut
	c1->cd( 3);
	name = "COMP_PhotonMmissProtOAP_v_TChanPhotProtOAP";
	TH2D *h2P = (TH2D*)full.Get( name);
	name = "COMP_PhotonMmissProtOAR_v_TChanPhotProtOAR";
	TH2D *h2R = (TH2D*)full.Get( name);
	h2P->Sumw2();
	h2P->Add( h2R, -0.0833);
	bin_lo = h2P->GetXaxis()->FindBin( chan_lo);
	bin_hi = h2P->GetXaxis()->FindBin( chan_hi);
	TH1D *hS = h2P->ProjectionY( "PhotonMmissProtOAP", bin_lo, bin_hi);

	hS->SetTitle();
	hS->SetLineWidth( 2);
	hS->GetXaxis()->SetTitleOffset( 1.1);
	hS->GetXaxis()->SetTitle( "Proton Missing Mass (MeV)");
	hS->GetXaxis()->SetLabelSize( 0.03);
	hS->GetYaxis()->SetLabelSize( 0.03);
	hS->GetXaxis()->CenterTitle();
	hS->GetXaxis()->SetRangeUser( 800, 1200);
	hS->Draw();

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hS, "require proton hit and angle cut", "p");
	tl->Draw();

	name = Form( "plots/Compton/phot_mmis_%d", ebin);
	name.Append( ".pdf");
	if ( print == kTRUE) c1->Print( name);
}

void PlotMmiss2D( UInt_t ebin = 8, Bool_t print = kFALSE)
{
	UInt_t i, chan_lo, chan_hi, bin_lo, bin_hi;
	Double_t x[2], y[2];
	Double_t energy, temin, temax;
	TString name;

	TFile* file;

	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

	// Middle channel bin and bin limits
	i = tchan[ebin];
	chan_lo = i-2;
	chan_hi = i+2;

	file = full;
	energy = WeightedPar( file, "energy", chan_lo, chan_hi);
	temin = tcd[chan_hi].energy - tcd[chan_hi].denergy;
	temax = tcd[chan_lo].energy + tcd[chan_lo].denergy;

	gStyle->SetOptStat( 0);

	TCanvas *c1 = new TCanvas ( "c1", "Canvas", 200, 0, 700, 500);

	name = "PhotProtOAP_v_PhotonMmissProtP_v_TChanPhotProtP";
	TH3D *h3P = (TH3D*)full.Get( name);
	name = "PhotProtOAR_v_PhotonMmissProtR_v_TChanPhotProtR";
	TH3D *h3R = (TH3D*)full.Get( name);
	h3P->Sumw2();
	h3P->Add( h3R, -0.0833);
	bin_lo = h3P->GetXaxis()->FindBin( chan_lo);
	bin_hi = h3P->GetXaxis()->FindBin( chan_hi);
	h3P->GetXaxis()->SetRange( chan_lo, chan_hi);
	TH2D *h2S = h3P->Project3D( "yz");

	name = Form( " E_{#gamma} = %5.1f (%5.1f - %5.1f) MeV", energy, temin,
			temax);
	h2S->SetTitle( name);
	h2S->GetXaxis()->SetTitleOffset( 1.3);
	h2S->GetYaxis()->SetTitleOffset( 1.2);
	h2S->GetXaxis()->SetTitle( "Scattered Photon CM Angle (deg)");
	h2S->GetYaxis()->SetTitle( "Proton Missing Mass (MeV)");
	h2S->GetXaxis()->SetLabelSize( 0.03);
	h2S->GetYaxis()->SetLabelSize( 0.03);
	h2S->GetXaxis()->CenterTitle();
	h2S->GetYaxis()->CenterTitle();
	h2S->RebinY( 20);
	h2S->GetXaxis()->SetRangeUser( 0, 20);
	h2S->GetYaxis()->SetRangeUser( 800, 1200);
//	h2S->Draw( "cont3");
	h2S->Draw( "box");
//	h2S->Draw( "lego");

/*
	x[0] = 0.2;
	x[1] = 0.4;
	y[0] = 0.88;
	y[1] = 0.96;

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hS, "no cut", "p");
	tl->Draw();
*/

//	l = new TLine( 5, 800, 5, 1200);
//	l->SetLineStyle( 2);
//	l->SetLineWidth( 2);
//	l->Draw();

	name = Form( "plots/Compton/phot_mmiss_2D_%d", ebin);
	name.Append( ".pdf");
	if ( print == kTRUE) c1->Print( name);
}
