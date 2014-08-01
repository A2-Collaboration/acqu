gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void EmPi02D( Int_t pbin, Char_t cond[32]="nmc")
{
	Int_t eg;
	Double_t a, b, c;
	Char_t name[64];

	if ( pbin == 2) eg = 147;
	else if ( pbin == 3) eg = 153;
	else if ( pbin == 4) eg = 157;
	else if ( pbin == 5) eg = 163;
	else if ( pbin == 6) eg = 168;
	else if ( pbin >= 146) eg = pbin;
	else {
		printf( "<E> CarbCont: Invalid Channel Number (2-6)\n");
		exit(-1);
	}

	a = Linear( 145, 918, 400, 933, eg);
	b = Linear( 145, 933, 400, 880, eg);
	c = Linear( 145, 945, 400, 1000, eg);

	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, eg);
	pdat = new TFile( name);
	h3p = (TH3D*) pdat->Get( "EmPi0_TGG_ThetaCM");
	TH2D *h2p = h3p->Project3D("xz1");
	h2p->SetMarkerColor( 4);

	h3pc = (TH3D*) pdat->Get( "EmPi0_TGG_ThetaCM_C");
	TH2D *h2pc = h3pc->Project3D("xz2");
	h2pc->SetMarkerColor( 4);

	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	cdat = new TFile( name);
	h3c = (TH3D*)cdat->Get( "EmPi0_TGG_ThetaCM");
	TH2D *h2c = h3c->Project3D("xz3");

	h3cc = (TH3D*)cdat->Get( "EmPi0_TGG_ThetaCM_C");
	TH2D *h2cc = h3cc->Project3D("xz4");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 800, 500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(2,1);

	c1->cd(1);
	h2p->GetXaxis()->SetRangeUser( 900, 980);
	h2p->SetTitle("Assume Proton Target");
	h2p->GetXaxis()->SetTitle("M_{inv} (MeV)");
	h2p->GetXaxis()->CenterTitle();
	h2p->GetYaxis()->SetTitle("#theta_{CM} (deg)");
	h2p->GetYaxis()->CenterTitle();
	h2p->GetYaxis()->SetTitleOffset(1.4);
	h2p->Draw();
	h2c->Draw("same");

	l = new TLine( a, 0, a, 180);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	l = new TLine( b, 0, b, 180);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	l = new TLine( c, 0, c, 180);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	pt = new TLegend(0.70,0.2,0.85,0.3);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->SetTextSize(0.04);
	h2p->SetFillColor(4);
	h2c->SetFillColor(1);
	pt->AddEntry( h2p, "Proton", "f");
	pt->AddEntry( h2c, "Carbon", "f");
	pt->Draw();

	c1->cd(2);
	h2pc->GetXaxis()->SetRangeUser( 11150, 11220);
	h2pc->SetTitle("Assume Carbon Target");
	h2pc->GetXaxis()->SetTitle("M_{inv} (MeV)");
	h2pc->GetXaxis()->CenterTitle();
	h2pc->GetYaxis()->SetTitle("#theta_{CM} (deg)");
	h2pc->GetYaxis()->CenterTitle();
	h2pc->GetYaxis()->SetTitleOffset(1.4);
	h2pc->Draw();
	h2cc->Draw("same");

	pt->Draw();

//	sprintf( name, "plots/eps/sim_eff_%d_%d.eps", eg, (theta_lo + theta_hi)/2);
	sprintf( name, "plots/pdf/sim_miss_%s_%d.pdf", cond, eg);
//	c1->Print( name);
}

void EmPi02DDat()
{
	TFile full( "histograms/ARH_full.root");
//	TFile empty( "histograms/ARH_empty.root");

	h3 = (TH3D*) full.Get( "EmPi0_TGG_ThetaCM_P");
	TH2D *h2 = h3->Project3D("xz1");
//	h2->SetMarkerColor( 4);

//	h3c = (TH3D*) empty.Get( "EmPi0_TGG_ThetaCM_R");
//	TH2D *h2c = h3c->Project3D("xz2");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	h2->GetXaxis()->SetRangeUser( 880, 950);
	h2->Draw();
//	h2c->Draw("same");

//	sprintf( name, "plots/eps/sim_eff_%d_%d.eps", eg, (theta_lo + theta_hi)/2);
//	c1->Print( name);
}
