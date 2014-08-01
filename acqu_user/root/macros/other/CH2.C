gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void CH2( Int_t eg = 150)
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, x3, a, b, fact;
	Char_t hist[32], cond[32], name[64];

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 200, 0, 700, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
//	c1->Divide(1,2);

	strcpy( hist, "EmPi0");

	c1->cd(1);
	sprintf( name, "histograms/MonteCarlo/pi0_p_%d.root", eg);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/pi0_c_%d.root", eg);
	carb = new TFile( name);

	h3p = (TH3D*) prot->Get( "EmPi0_TGG_ThetaCM");
	TH1D *h1 = h3p->Project3D("z");
	h3c = (TH3D*) carb->Get( "EmPi0_TGG_ThetaCM");
	TH1D *h2 = h3c->Project3D("z2");

	sprintf( name, "E_{#gamma} = %5.1f MeV, %s target", eg, cond);
	h1->SetTitle( name);

	x1 = Linear( 145, 933, 400, 880, eg);
	x2 = Linear( 145, 945, 400, 1000, eg);
	x3 = Linear( 145, 918, 400, 880, eg);
	a = Linear( 145, 900, 400, 600, eg);
	b = Linear( 145, 950, 400, 1000, eg);

	col = 2;

	h1->SetLineWidth(2);
	h1->SetFillColor(2);
	h1->GetXaxis()->SetRangeUser(a,b);
	h1->Draw();

	ch1 = h1->Integral();
	ch2 = h2->Integral();
	fact = (double) ch1/ch2;

	Int_t chan1 = h1->GetXaxis()->FindBin( x1);
	Int_t chan2 = h1->GetXaxis()->FindBin( x2);
	Int_t chan3 = h1->GetXaxis()->FindBin( x3);
	Double_t intp = h1->Integral( chan1, chan2);
	Double_t intc = h2->Integral( chan3, chan1);
	Double_t ratio = intp/intc;
	cout << " 2 cm target" << endl;
	cout << intp << " " << intc << endl;
	cout << "ratio = " << ratio << endl;

	h2->SetFillStyle(3001);
	h2->SetFillColor(4);
	h2->SetLineWidth(2);
	h2->SetLineStyle(1);
	h2->Scale( fact);
	h2->Draw("same");

	TH1D *h3 = (TH1D*)h1->Clone();
	h3->Sumw2();
	h3->Add(h2, 1.0);
	h3->Draw("same");

	max = h3->GetMaximum();

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

	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineColor(col);
	l3->SetLineStyle(2);
	l3->SetLineWidth(4);
	l3->Draw();

	sprintf( name, "plots/pdf/EmPi0_ch2_%d.pdf", eg);
//	c1->Print( name);
}
