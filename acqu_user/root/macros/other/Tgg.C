gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

typedef struct {

	Int_t egamma;
	Double_t energy;
	Double_t denergy;
	Double_t etag;
	Double_t d_etag;
	Double_t edet[10];
	Double_t d_edet[10];
	Double_t pa_ratio;
	Double_t xs;
	Double_t dxs;

} TChanData;
TChanData tcd[352];

void Init()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");
}


void PlotTggCut( Double_t Eg = 150)
{
	Double_t min, offset;
	TString tgt, name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 800, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	c1->cd( 1);
	TF1 *f1 = new TF1( "f1", TggOut, 0, 180, 3);
	offset = 0;
	f1->SetParameters( Eg, kM_C12_MEV, offset);
	TF1 *f2 = new TF1( "f2", TggOut, 0, 180, 3);
	offset = 3;
	f2->SetParameters( Eg, kMP_MEV, offset);
	f2->SetLineStyle( 2);

	name = Form( "histograms/MonteCarlo/pi0/pi0_p_%d.root", (int) Eg);
	simp = new TFile( name);
	name = "THR_TGGP_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)simp->Get( name);

	min = f2->Eval( 0) - 10;
	h2sim->GetYaxis()->SetRangeUser( min, 180);
	h2sim->Draw( "zcol");

	f1->Draw( "same");
	f2->Draw( "same");

	// Cut Efficiency - Proton
	Int_t cts, cts_cut1;
	Double_t eff;
	name = "THR_ThetaCMP";
	TH1D *h1 = (TH1D*)simp->Get( name);
	cts = h1->Integral();
	name = "THR_ThetaCMCut2P";
	TH1D *h1_cut1 = (TH1D*)simp->Get( name);
	cts_cut1 = h1_cut1->Integral();

	eff = (double) cts_cut1/cts;
	cout << cts_cut1;
	cout << " " << cts;
	cout << " " << eff;
	cout << endl;

	c1->cd( 2);

	name = Form( "histograms/MonteCarlo/pi0/pi0_c_%d.root", (int) Eg);
	simc = new TFile( name);
	name = "THR_TGGP_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)simc->Get( name);

	min = f2->Eval( 0) - 10;
	h2sim->GetYaxis()->SetRangeUser( min, 180);
	h2sim->Draw( "zcol");

	f1->Draw( "same");
	f2->Draw( "same");

	// Cut Efficiency - C12
	name = "THR_ThetaCMP";
	TH1D *h1c = (TH1D*)simc->Get( name);
	cts = h1c->Integral();
	name = "THR_ThetaCMCut2P";
	TH1D *h1c_cut1 = (TH1D*)simc->Get( name);
	cts_cut1 = h1c_cut1->Integral();

	eff = (double) cts_cut1/cts;
	cout << cts_cut1;
	cout << " " << cts;
	cout << " " << eff;
	cout << endl;

}

Double_t TggOut( Double_t *x, Double_t *par)
{
	Double_t thcm, Eg, mtgt, Tgg, offset;
	Double_t q_pi, T_pi;

	thcm = x[0];

	Eg = par[0];
	mtgt = par[1];
	offset = par[2];

	q_pi = qp_thcm( Eg, mtgt, thcm, kMPI0_MEV);
	T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
	Tgg = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

	return( Tgg - offset);
}

void PlotTggCutChan( UInt_t chan = 300)
{
	UInt_t eg;
	Double_t min, offset, Eg;
	TString tgt, name;

	Eg = tcd[chan].energy;
	eg = tcd[chan].egamma;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 800, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	gStyle->SetOptStat( 0);

	c1->cd( 1);
	TF1 *f1 = new TF1( "f1", TggOut, 0, 180, 3);
	offset = 0;
	f1->SetParameters( Eg, kM_C12_MEV, offset);
	TF1 *f2 = new TF1( "f2", TggOut, 0, 180, 3);
	offset = 3;
	f2->SetParameters( Eg, kMP_MEV, offset);
	f2->SetLineStyle( 2);

	name = Form( "histograms/MonteCarlo/pi0/LH2_10cm/pi0_p_chan%d.root", chan);
	simp = new TFile( name);
	name = "THR_TGGCut1P_v_ThetaCMCut1P";
	TH2D *h2sim = (TH2D*)simp->Get( name);

	min = f2->Eval( 0) - 10;
	h2sim->GetYaxis()->SetRangeUser( min, 180);
	h2sim->SetTitle( "proton");
	h2sim->GetXaxis()->SetTitle( "#theta^{CM}_{#pi} (deg)");
	h2sim->GetXaxis()->CenterTitle();
	h2sim->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h2sim->GetYaxis()->CenterTitle();
	h2sim->GetYaxis()->SetTitleOffset( 1.5);
	h2sim->Draw( "zcol");

//	f1->Draw( "same");
	f2->Draw( "same");

	// Cut Efficiency - Proton
	Int_t cts, cts_cut1;
	Double_t eff;
	name = "THR_ThetaCMCut1P";
	TH1D *h1 = (TH1D*)simp->Get( name);
	cts = h1->Integral();
	name = "THR_ThetaCMCut2P";
	TH1D *h1_cut1 = (TH1D*)simp->Get( name);
	cts_cut1 = h1_cut1->Integral();

	eff = (double) cts_cut1/cts;
	cout << cts_cut1;
	cout << " " << cts;
	cout << " " << eff;
	cout << endl;

	name = Form( "Efficiency = %6.4f\n", eff);
	pl = new TPaveLabel( 0.5, 0.1, 0.7, 0.2, name, "NDC");
	pl->SetBorderSize( 0);
	pl->SetFillColor( 0);
	pl->SetTextSize( 0.5);
	pl->Draw();

	c1->cd( 2);

	name = Form( "histograms/MonteCarlo/pi0/LH2_10cm/pi0_w_chan%d.root", chan);
	simc = new TFile( name);
	name = "THR_TGGCut1P_v_ThetaCMCut1P";
	TH2D *h2sim = (TH2D*)simc->Get( name);

	min = f2->Eval( 0) - 10;
	h2sim->GetYaxis()->SetRangeUser( min, 180);
	h2sim->SetTitle( "carbon");
	h2sim->GetXaxis()->SetTitle( "#theta^{CM}_{#pi} (deg)");
	h2sim->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h2sim->GetYaxis()->SetTitleOffset( 1.5);
	h2sim->GetXaxis()->CenterTitle();
	h2sim->GetYaxis()->CenterTitle();
	h2sim->Draw( "zcol");

//	f1->Draw( "same");
	f2->Draw( "same");

	// Cut Efficiency - C12
	name = "THR_ThetaCMCut1P";
	TH1D *h1c = (TH1D*)simc->Get( name);
	cts = h1c->Integral();
	name = "THR_ThetaCMCut2P";
	TH1D *h1c_cut1 = (TH1D*)simc->Get( name);
	cts_cut1 = h1c_cut1->Integral();

	eff = (double) cts_cut1/cts;
	cout << cts_cut1;
	cout << " " << cts;
	cout << " " << eff;
	cout << endl;

	name = Form( "Efficiency = %6.4f\n", eff);
	pl = new TPaveLabel( 0.5, 0.1, 0.7, 0.2, name, "NDC");
	pl->SetBorderSize( 0);
	pl->SetFillColor( 0);
	pl->SetTextSize( 0.5);
	pl->Draw();

	name = Form( "plots/Pi0/TggCut_eff_chan%d.pdf", chan);
	c1->Print( name);
}

void TggChanSim( TString tgt = "subt", UInt_t i = 300, Double_t theta = 90,
		Double_t thetabin = 1)
{
	UInt_t rebin, tbin_lo, tbin_hi;
	Double_t q_pi, T_pi, x1, x2, x3, left, right, eg, deg, max;
	Double_t thetalo, thetahi;
	Char_t tchar;
	TString name;

	thetalo = theta - thetabin/2;
	thetahi = theta + thetabin/2;

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 600, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	rebin = 1;

	if ( tgt == "subt") tchar = 'p';
	else if ( tgt == "empty") tchar = 'w';
	name = Form( "histograms/MonteCarlo/pi0/LH2_10cm/pi0_%c_chan%d.root",
			tchar, i);

	sim = new TFile( name);
	name = "THR_TGGP_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)sim->Get( name);

	// Proton Opening Angle
	q_pi = qp_thcm( eg, kMP_MEV, thetalo, kMPI0_MEV);
	T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
	x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;
	cout << "tgg_min = " << x1;
//	x1 -= TggCutOffset;
	cout << " tgg_cut = " << x1;

	// 12-C Opening Angle
	q_pi = qp_thcm( eg, kM_C12_MEV, thetalo, kMPI0_MEV);
	T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
	x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

	cout << " tgg_c12 = " << x2;
	cout << endl;

	x3 = 180;

//	left = x2 - 10;
	left = 90;
	right = 180;

	tbin_lo = h2sim->GetXaxis()->FindBin( thetalo);
	tbin_hi = h2sim->GetXaxis()->FindBin( thetahi)-1;
	h2sim->GetXaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *h1sim = h2sim->ProjectionY( "simhist");
	h1sim->Rebin( rebin);

	h1sim->SetTitle();

	name = "#theta_{#gamma#gamma} (deg)";
	h1sim->GetXaxis()->SetTitle( name);
	h1sim->GetXaxis()->SetRangeUser( left, right);

	h1sim->SetMarkerStyle( 21);
	h1sim->SetMarkerSize( 0.5);
	h1sim->SetLineWidth( 2);
	h1sim->Draw( "E");
	c1->Update();
	max = gPad->GetUymax();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineWidth( 2);
	l1->SetLineStyle( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineWidth( 2);
	l2->SetLineStyle( 2);
	l2->Draw();

	name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.5);
	pl->SetTextFont( 132);
	pl->Draw();

	name = Form ( "#theta^{#pi}_{CM} = %5.1f #pm %4.1f deg", theta, thetabin/2);
	pl2 = new TPaveLabel( 0.55, 0.90, 0.85, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.5);
	pl2->SetTextFont( 132);
	pl2->Draw();

	name = Form( "plots/Pi0/Tgg_sim_chan%d_t%d.pdf", i, (int) theta);
	c1->Print( name);
}
