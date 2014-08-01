gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

//TFile full( "histograms/ARH_full.root");
//TFile full( "histograms/Full/ARH.root");
TFile full( "histograms/Compton/Full/Full.root");

void Plot()
{
	TH1D *h1 = (TH1D*)full.Get( "THR_PhotonTimeDiff");
	TH1D *h2 = (TH1D*)full.Get( "THR_M2g");

	TCanvas *c1 = new TCanvas ( "c1", "Canvas", 200, 0, 900, 500);
   c1->SetFillColor( 38);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(2,1);

	Double_t ym = h1->GetMaximum();
	ym *= 1.2;
	h1->SetMaximum( ym);
	l1 = new TLine(-18, 0,-18, ym);
	l1->SetLineColor(4);
	l1->SetLineStyle(2);
	l1->SetLineWidth(4);
	l2 = new TLine(18, 0,18, ym);
	l2->SetLineColor(4);
	l2->SetLineStyle(2);
	l2->SetLineWidth(4);

	gStyle->SetOptStat(0);
	h1->SetTitle();
	h1->SetFillStyle(3001);
	h1->SetFillColor(2);
	h1->SetLineColor(1);
	h1->SetLineWidth(2);
	h1->GetXaxis()->SetTitle("t_{#gamma#gamma} (ns)");
	h1->GetXaxis()->CenterTitle();
	h1->GetXaxis()->SetLabelSize( 0.035);
	h1->GetYaxis()->SetLabelSize( 0.035);
	h1->GetXaxis()->SetRangeUser(-80,80);

	Double_t ym = h2->GetMaximum();
	ym *= 1.2;
	h2->SetMaximum( ym);
	l3 = new TLine(115, 0,115, ym);
	l3->SetLineColor(4);
	l3->SetLineStyle(2);
	l3->SetLineWidth(4);
	l4 = new TLine(155, 0,155, ym);
	l4->SetLineColor(4);
	l4->SetLineStyle(2);
	l4->SetLineWidth(4);

	h2->SetTitle();
	h2->SetFillStyle(3001);
	h2->SetFillColor(2);
	h2->SetLineColor(1);
	h2->SetLineWidth(2);
	h2->GetXaxis()->SetTitle("M_{#gamma#gamma} (MeV)");
	h2->GetXaxis()->CenterTitle();
	h2->GetXaxis()->SetLabelSize( 0.035);
	h2->GetYaxis()->SetLabelSize( 0.035);
	h2->GetXaxis()->SetRangeUser(0,300);

	c1->cd(1);
	gPad->SetLogy(1);
	h1->Draw();
	l1->Draw();
	l2->Draw();

	c1->cd(2);
	h2->Draw();
	l3->Draw();
	l4->Draw();

	c1->Print( "plots/eps/inv_tgg.eps");
}

void TDC( UInt_t chan1, Int_t chan2)
{
	UInt_t prompt[2], acc[2], bin;
	Double_t yy, ym, y1, y2;
	TString name;

	prompt[0] = 94;
	prompt[1] = 110;
	acc[0] = 120;
	acc[1] = 168;

	TCanvas *c1 = new TCanvas ( "c1", "Canvas", 200, 0, 900, 500);
   c1->SetFillColor( 38);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(2,1);

	// Pi0 or Phot Tagger Time
//	name = "THR_TChanHit_v_Pi0TaggTime";
	name = "COMP_TChanHit_v_PhotTaggTime";
	TH2D *h1 = (TH2D*)full.Get( name);
	bin = h1->GetYaxis()->FindBin( chan1);
	TH1D *p1 = h1->ProjectionX( "chan1", bin, bin);
	bin = h1->GetYaxis()->FindBin( chan2);
	TH1D *p2 = h1->ProjectionX( "chan2", bin, bin);

	// TDC Spectrum
	c1->cd(1);
	p1->Draw();
	sprintf( name, "Channel %d", chan1);
	p1->SetTitle( name);
	p1->SetLineWidth(2);
	p1->GetXaxis()->SetTitleOffset( 1.1);
	p1->GetXaxis()->SetTitle("#pi^{0} - tagger time (ns)");
	p1->GetXaxis()->SetLabelSize( 0.03);
	p1->GetYaxis()->SetLabelSize( 0.03);
	p1->GetXaxis()->CenterTitle();

	yy = p1->GetMaximum();
	ym = yy*1.05;
	p1->SetMaximum(ym);
	y1 = yy*0.9;
	y2 = yy;

	l1 = new TLine( prompt[0], 0,prompt[0], ym);
	l1->SetLineColor(2);
	l1->SetLineStyle(2);
	l1->SetLineWidth(3);
	l1->Draw();
	l2 = new TLine( prompt[1], 0,prompt[1], ym);
	l2->SetLineColor(2);
	l2->SetLineStyle(2);
	l2->SetLineWidth(3);
	l2->Draw();
	l3 = new TLine( acc[0], 0,acc[0], ym);
	l3->SetLineColor(4);
	l3->SetLineStyle(2);
	l3->SetLineWidth(3);
	l3->Draw();
	l4 = new TLine( acc[1], 0,acc[1], ym);
	l4->SetLineColor(4);
	l4->SetLineStyle(2);
	l4->SetLineWidth(3);
	l4->Draw();

	// TDC Spectrum
	c1->cd(2);
	sprintf( name, "Channel %d", chan2);
	p2->SetTitle( name);
	p2->SetLineWidth(2);
	p2->Draw();
	p2->GetXaxis()->SetTitleOffset( 1.1);
	p2->GetXaxis()->SetTitle("#pi^{0} - tagger time (ns)");
	p2->GetXaxis()->SetLabelSize( 0.03);
	p2->GetYaxis()->SetLabelSize( 0.03);
	p2->GetXaxis()->CenterTitle();

	yy = p2->GetMaximum();
	ym = yy*1.05;
	p2->SetMaximum(ym);
	y1 = yy*0.9;
	y2 = yy;

	l1 = new TLine( prompt[0], 0,prompt[0], ym);
	l1->SetLineColor(2);
	l1->SetLineStyle(2);
	l1->SetLineWidth(3);
	l1->Draw();
	l2 = new TLine( prompt[1], 0,prompt[1], ym);
	l2->SetLineColor(2);
	l2->SetLineStyle(2);
	l2->SetLineWidth(3);
	l2->Draw();
	l3 = new TLine( acc[0], 0,acc[0], ym);
	l3->SetLineColor(4);
	l3->SetLineStyle(2);
	l3->SetLineWidth(3);
	l3->Draw();
	l4 = new TLine( acc[1], 0,acc[1], ym);
	l4->SetLineColor(4);
	l4->SetLineStyle(2);
	l4->SetLineWidth(3);
	l4->Draw();

	sprintf( name, "plots/eps/tdc_%d_%d.eps", chan1, chan2);
	c1->Print( name);
}

void PlotAll()
{
	UInt_t prompt[2], acc[4];
	Double_t yy, ym, y1, y2;
	TString name;

	gStyle->SetOptStat( 0);

	prompt[0] = 92;
	prompt[1] = 106;
	acc[0] = 20;
	acc[1] = 80;
	acc[2] = 118;
	acc[3] = 178;

	TCanvas *c1 = new TCanvas ( "c1", "Canvas", 200, 0, 1200, 500);
//	c1->SetFillColor( 38);
//	c1->SetGrid();
//	c1->GetFrame()->SetFillColor( 21);
//	c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 3,1);

	// Phot Tagger Time
	c1->cd( 1);
	name = "COMP_PhotTaggTime";
	TH1D *h1 = (TH1D*)full.Get( name);

	// TDC Spectrum
	h1->SetTitle();
	h1->SetLineWidth(2);
	h1->GetXaxis()->SetTitleOffset( 1.1);
	h1->GetXaxis()->SetTitle("photon - tagger time (ns)");
	h1->GetXaxis()->SetLabelSize( 0.03);
	h1->GetYaxis()->SetLabelSize( 0.03);
	h1->GetXaxis()->CenterTitle();
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

	// Phot Tagger Time
	c1->cd( 2);
	name = "COMP_PhotonMmissP";
	TH1D *hP = (TH1D*)full.Get( name);
	name = "COMP_PhotonMmissR";
	TH1D *hR = (TH1D*)full.Get( name);
	TH1D *hS = (TH1D*)hP->Clone( "subt");

	// TDC Spectrum
	hP->SetTitle();
	hP->SetLineWidth( 2);
	hP->GetXaxis()->SetTitleOffset( 1.1);
	hP->GetXaxis()->SetTitle("Proton Missing Mass (MeV)");
	hP->GetXaxis()->SetLabelSize( 0.03);
	hP->GetYaxis()->SetLabelSize( 0.03);
	hP->GetXaxis()->CenterTitle();

	hR->SetLineWidth( 2);
	hR->SetLineStyle( 2);
	hR->Scale( 0.0833);
	hP->Draw();
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
	hS->GetXaxis()->SetTitle("Proton Missing Mass (MeV)");
	hS->GetXaxis()->SetLabelSize( 0.03);
	hS->GetYaxis()->SetLabelSize( 0.03);
	hS->GetXaxis()->CenterTitle();
	hS->SetLineWidth( 2);
	hS->SetMarkerStyle( 20);
	hS->SetMarkerSize( 0.5);
	hS->Draw();

	tl = new TLegend( 0.2, 0.7, 0.4, 0.8);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( hS, "subtracted", "p");
	tl->Draw();

	name = "plots/Compton/random_subt.pdf";
	c1->Print( name);
}
