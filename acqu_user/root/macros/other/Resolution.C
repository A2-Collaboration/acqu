gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

typedef struct {

	Double_t mid;
	Double_t mean;
	Double_t dmean;
	Double_t sigma;
	Double_t dsigma;

} TData;

TData data[5][9];

void PlotAll( TString type = "theta", Int_t ebin, Bool_t save = kFALSE)
{
	Int_t tbin;

	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *canvas = new TCanvas ( "c1", "Resolutions", 200, 0, 1000, 700);

	canvas->SetFillColor( 42);
	canvas->GetFrame()->SetFillColor( 21);
	canvas->GetFrame()->SetBorderSize( 12);
	canvas->Divide(3,3);

	for ( tbin = 0; tbin < 9; tbin++) {
		canvas->cd( tbin+1);
		if ( type == "theta") TRes( ebin, tbin, kFALSE);
		else ERes( ebin, tbin, kFALSE);
	}

	if ( save == kTRUE) {

		if ( type == "theta") {
			name = Form( "plots/TRes_%d.pdf", 150+15*ebin);
			c1->Print( name);
			name = Form( "xs/TRes_%d.out", 150+15*ebin);
		}
		else {
			name = Form( "plots/ERes_%d.pdf", 150+15*ebin);
			c1->Print( name);
			name = Form( "xs/ERes_%d.out", 150+15*ebin);
		}


		ofstream outFile( name);
		if ( !outFile.is_open()) {
			cout << "Error opening file ";
			cout << name;
			cout << endl;
			break;
		}
		for ( tbin = 0; tbin < 9; tbin++) {
			name = Form( "%3d  %5.1f  %6.2f +/- %4.2f  %5.2f +/- %4.2f",
					10+20*tbin, data[ebin][tbin].mid, data[ebin][tbin].mean,
					data[ebin][tbin].dmean, data[ebin][tbin].sigma,
					data[ebin][tbin].dsigma);
			outFile << name << endl;
		}
		outFile.close();
	}
}

void TRes( Int_t ebin, Int_t tbin, Bool_t save = kFALSE)
{
	Int_t eg, theta;
	Double_t left, right, max;
	TString name;

	eg = 150 + 15*ebin;
	theta = 10 + 20*tbin;

	TF1* gaus;

	name = Form( "histograms/MonteCarlo/pi0_e%d_t%d.root", eg, theta);
	TFile sim( name);

	name = "EFF_Pi0Theta";
	hist = (TH1D*) sim.Get( name);

	name = Form( "tbin%d", tbin);
	hist->SetName( name);

	TCanvas* canv;
	if ( !( canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		TCanvas* canv = new TCanvas ( "c1", "CH2 Comparison", 200, 0, 700, 500);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.035, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	left = theta - 20;
	right = theta + 20;

	if ( tbin == 0) name = Form( "E_{#gamma} = %d MeV", eg);
	else name = "";
	hist->SetTitle( name);
	hist->SetLineWidth( 3);
	hist->SetLineColor( 2);
	max = hist->GetMaximum();
	max *= 1.10;
	hist->SetMaximum( max);
	hist->Fit( "gaus", "Q");
	hist->GetXaxis()->SetRangeUser( left, right);
	hist->GetXaxis()->CenterTitle();
	hist->GetXaxis()->SetTitle( "#theta_{#pi} (deg)");
	hist->DrawCopy();

	gaus = hist->GetFunction( "gaus");

	data[ebin][tbin].mid = theta;
	data[ebin][tbin].mean = gaus->GetParameter(1);
	data[ebin][tbin].dmean = gaus->GetParError(1);
	data[ebin][tbin].sigma = gaus->GetParameter(2);
	data[ebin][tbin].dsigma = gaus->GetParError(2);

	Double_t x[2], y[2];

	y[0] = 0.65;
	y[1] = 0.85;
	x[0] = 0.60;
	x[1] = 0.80;
	if ( theta >= 170) {
		x[0] = 0.15;
		x[1] = 0.35;
	}
	pl = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pl->SetBorderSize(0);
	pl->SetFillStyle(0);
	pl->SetTextAlign(12);
	pl->SetTextSize(0.04);
	name = Form( "#theta_{#pi} = %d deg", theta);
	pl->AddText( name);
	name = Form( "#theta_{0} = %6.2f #pm %4.2f deg\n", data[ebin][tbin].mean,
			data[ebin][tbin].dmean);
	pl->AddText( name);
	name = Form( "#sigma = %5.2f #pm %4.2f deg\n", data[ebin][tbin].sigma,
			data[ebin][tbin].dsigma);
	pl->AddText( name);
	pl->Draw();

	l1 = new TLine( theta, 0, theta, max);
	l1->SetLineColor( 1);
	l1->SetLineStyle( 2);
	l1->SetLineWidth( 4);
	l1->Draw();

	name = Form( "plots/TRes.pdf");
	if ( save == kTRUE) c1->Print( name);
}

void ERes( Int_t ebin, Int_t tbin, Bool_t save = kFALSE)
{
	Int_t eg, theta;
	Double_t left, right, max;
	Double_t width;
	Double_t mom, KE;
	TString name;

	eg = 150 + 15*ebin;
	theta = 10 + 20*tbin;

	width = 10;
	if ( ebin > 1) width = 40;

	mom = qp( (double) eg, kMP_MEV, (double) theta*kD2R, kMPI0_MEV);
	KE = Energy( mom, kMPI0_MEV) - kMPI0_MEV;

	TF1* gaus;

	name = Form( "histograms/MonteCarlo/pi0_e%d_t%d.root", eg, theta);
	TFile sim( name);

	name = "EFF_Pi0KE";
	hist = (TH1D*) sim.Get( name);

	TCanvas* canv;
	if ( !( canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		TCanvas* canv = new TCanvas ( "c1", "CH2 Comparison", 200, 0, 700, 500);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.035, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	left = KE - width;
	right = KE + width;

	if ( tbin == 0) name = Form( "E_{#gamma} = %d MeV", eg);
	else name = "";
	hist->SetTitle( name);
	hist->SetLineWidth( 3);
	hist->SetLineColor( 2);
	max = hist->GetMaximum();
	max *= 1.10;
	hist->SetMaximum( max);
	hist->Fit( "gaus", "Q");
	hist->GetXaxis()->SetRangeUser( left, right);
	hist->GetXaxis()->CenterTitle();
	hist->GetXaxis()->SetTitle( "T_{#pi} (MeV)");
	hist->DrawCopy();

	gaus = hist->GetFunction( "gaus");

	data[ebin][tbin].mid = KE;
	data[ebin][tbin].mean = gaus->GetParameter(1);
	data[ebin][tbin].dmean = gaus->GetParError(1);
	data[ebin][tbin].sigma = gaus->GetParameter(2);
	data[ebin][tbin].dsigma = gaus->GetParError(2);

	Double_t x[2], y[2];

	y[0] = 0.60;
	y[1] = 0.90;
	x[0] = 0.60;
	x[1] = 0.80;
	pl = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pl->SetBorderSize(0);
	pl->SetFillStyle(0);
	pl->SetTextAlign(12);
	pl->SetTextSize(0.04);
	name = Form( "#theta_{#pi} = %d deg", theta);
	pl->AddText( name);
	name = Form( "T_{#pi} = %5.1f MeV", KE);
	pl->AddText( name);
	name = Form( "T_{0} = %6.2f #pm %4.2f MeV\n", data[ebin][tbin].mean,
			data[ebin][tbin].dmean);
	pl->AddText( name);
	name = Form( "#sigma = %5.2f #pm %4.2f MeV\n", data[ebin][tbin].sigma,
			data[ebin][tbin].dsigma);
	pl->AddText( name);
	pl->Draw();

	l1 = new TLine( KE, 0, KE, max);
	l1->SetLineColor( 1);
	l1->SetLineStyle( 2);
	l1->SetLineWidth( 4);
	l1->Draw();

	name = Form( "plots/ERes.pdf");
	if ( save == kTRUE) c1->Print( name);
}
