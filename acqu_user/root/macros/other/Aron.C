gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

Int_t eg[] = {150, 165, 180, 195, 210};
Double_t fact[] = {60, 21, 16, 12, 8};

void MmissPlot( UInt_t ebin, Bool_t save = kFALSE, Bool_t cut = kTRUE,
		Double_t TggOffset = 0)
{
	UInt_t tbin;

	c1 = new TCanvas ( "c1", "CH2 Comparison", 200, 0, 1000, 700);

	c1->SetFillColor( 42);
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 3,3);

	for ( tbin = 0; tbin < 9; tbin++) {
		c1->cd( tbin+1);
//		MmissLoop( ebin, tbin, kFALSE, cut, TggOffset);
		MmissTheta( ebin, tbin, kFALSE, cut, TggOffset);
	}

	TString name = Form( "plots/Aron/Mmiss_%d.pdf", eg[ebin]);
	if ( save == kTRUE) c1->Print( name);
}

void MmissLoop( UInt_t ebin = 0, UInt_t tbin = 0, Bool_t save = kFALSE,
		Bool_t cut = kTRUE, Double_t TggOffset = 0)
{
	Double_t theta_lo;
	Double_t mPi0, mProt;
	Double_t q_pi0, T_pi0;
	Double_t Tgg, TggCut, msep;
	Double_t a, b;
	Double_t max;
	TString name;

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.030, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	TCanvas* canv;
	if ( !( canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		TCanvas* canv = new TCanvas ( "c1", "CH2 Comparison", 200, 0, 700, 500);

	theta_lo = tbin*20;

	mPi0 = kMPI0_MEV;
	mProt = kMP_MEV;

	q_pi0 = qp_thcm( eg[ebin], mProt, theta_lo, mPi0);
	T_pi0 = Energy( q_pi0, mPi0) - mPi0;
	Tgg = Tgg_Min( T_pi0, mPi0)/kD2R;
	if ( theta_lo < 90) TggCut = Tgg - TggOffset;
	else TggCut = Tgg;

	if ( ebin == 0) msep = 934;
	else if ( ebin == 1) msep = 934;
	else if ( ebin == 2) msep = 932;
	else if ( ( ebin == 3) || ( ebin == 4)) msep = 930;

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_p_%d.root", eg[ebin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_c_%d.root", eg[ebin]);
	carb = new TFile( name);

	name = "MMissP_v_TGGP_v_ThetaCMP";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

	if ( cut == kTRUE) {
		h3p->GetYaxis()->SetRangeUser( TggCut, 180);
		h3c->GetYaxis()->SetRangeUser( TggCut, 180);
	}
	h3p->GetXaxis()->SetRangeUser( theta_lo, 180);
	h3c->GetXaxis()->SetRangeUser( theta_lo, 180);

	TH1D *h1p = h3p->Project3D( "z");
	TH1D *h1c = h3c->Project3D( "z2");

	name = Form( "prot%d", tbin);
	TH1D *hp = (TH1D*)h1p->Clone( name);
	name = Form( "carb%d", tbin);
	TH1D *hc = (TH1D*)h1c->Clone( name);

	hc->SetFillStyle( 3001);
	hc->SetFillColor( 4);
	hc->SetLineWidth( 2);
	hc->SetLineStyle( 1);

	Double_t totc = hc->Integral();
	Double_t effc = totc/1e5;

	hp->SetLineWidth( 2);
	hp->SetFillColor( 2);
	Double_t totp = hp->Integral();
	Double_t effp = totp/1e5;
	hp->Scale( 1/fact[ebin]);

	a = Linear( 145, 900, 400, 600, eg[ebin]);
	b = Linear( 145, 950, 400, 1000, eg[ebin]);

	if ( tbin == 0)
		name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	else
		name = "";

	TH1D *hs = (TH1D*)hp->Clone();
	hs->SetTitle( name);
	hs->Sumw2();
	hs->Add( hc, 1.0);
	hs->GetXaxis()->SetRangeUser( a, b);
	hs->GetXaxis()->SetTitle( "M_{miss} (MeV)");

	max = hs->GetMaximum();

	hs->Draw();
	hp->Draw( "same");
	hc->Draw( "same");

	l1 = new TLine( msep, 0, msep, max);
	l1->SetLineColor( 1);
	l1->SetLineStyle( 2);
	l1->SetLineWidth( 4);
	l1->Draw();

	pl = new TPaveText( 0.15, 0.60, 0.45, 0.85, "NDC");
	pl->SetBorderSize( 0);
	pl->SetFillStyle( 0);
	pl->SetTextAlign( 12);
	pl->SetTextSize( 0.04);
	name = Form( "#theta_{#gamma#gamma} = %d - 180 deg\n", (int) TggCut);
	pl->AddText( name);
	name = Form( "#theta^{#pi}_{cm} = %d - 180 deg\n", (int) theta_lo);
	pl->AddText( name);
	name = Form( "proton eff = %5.3f\n", effp);
	pl->AddText( name);
	name = Form( "carbon eff = %5.3f\n", effc);
	pl->AddText( name);
	pl->Draw();

	if ( tbin == 0) {
		tl = new TLegend( 0.15, 0.2, 0.45, 0.4);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.04);

		tl->AddEntry( hc, "^{12}C");
		tl->AddEntry( hp, "p");
		tl->AddEntry( hs, "sum", "P");

		tl->Draw();
	}

	name = Form( "plots/Aron/Mmiss_%d_%d.pdf", eg[ebin], (int) theta_lo);
	if ( save == kTRUE) canv->Print( name);
}

void MmissTheta( Int_t ebin = 0, Int_t tbin = 0, Bool_t save = kFALSE,
		Bool_t cut = kTRUE, Double_t TggOffset = 0)
{
	Int_t th_lo, th_mid, th_hi;
	Int_t bin_lo, bin_hi;
	Double_t a, b;
	Double_t max;
	Double_t mPi0, mProt;
	Double_t q_pi0, T_pi0;
	Double_t Tgg, TggCut, msep;
	TString name;

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_p_%d.root", eg[ebin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_c_%d.root", eg[ebin]);
	carb = new TFile( name);

	th_lo = tbin*20;
	th_mid = 10 + tbin*20;
	th_hi = 20 + tbin*20;

	mPi0 = kMPI0_MEV;
	mProt = kMP_MEV;

	q_pi0 = qp_thcm( eg[ebin], mProt, th_lo, mPi0);
	T_pi0 = Energy( q_pi0, mPi0) - mPi0;
	Tgg = Tgg_Min( T_pi0, mPi0)/kD2R;
	if ( tbin < 5) TggCut = Tgg - TggOffset;
	else TggCut = Tgg;

	if ( ebin == 0) msep = 934;
	else if ( ebin == 1) msep = 934;
	else if ( ebin == 2) msep = 932;
	else if ( ( ebin == 3) || ( ebin == 4)) msep = 930;

	name = "MMissP_v_TGGP_v_ThetaCMP";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

	bin_lo = h3p->GetXaxis()->FindBin( th_lo);
	bin_hi = h3p->GetXaxis()->FindBin( th_hi)-1;
	h3p->GetXaxis()->SetRange( bin_lo, bin_hi);
	h3c->GetXaxis()->SetRange( bin_lo, bin_hi);

	bin_lo = h3p->GetYaxis()->FindBin( TggCut);
	bin_hi = h3p->GetYaxis()->FindBin( 180)-1;
	if ( cut == kTRUE) {
		h3p->GetYaxis()->SetRangeUser( bin_lo, bin_hi);
		h3c->GetYaxis()->SetRangeUser( bin_lo, bin_hi);
	}

	TH1D *h1p = h3p->Project3D( "z");
	TH1D *h1c = h3c->Project3D( "z2");

	name = Form( "prot%d", tbin);
	TH1D *hp = (TH1D*)h1p->Clone( name);
	name = Form( "carb%d", tbin);
	TH1D *hc = (TH1D*)h1c->Clone( name);

	hc->SetFillStyle( 3001);
	hc->SetFillColor( 4);
	hc->SetLineWidth( 2);
	hc->SetLineStyle( 1);

	hp->SetLineWidth( 2);
	hp->SetFillColor( 2);
	hp->Scale( 1/fact[ebin]);

	bin_lo = hp->GetXaxis()->FindBin( msep);
	bin_hi = hp->GetXaxis()->FindBin( 1000);
	Double_t cutp = hp->Integral( bin_lo, bin_hi);
	Double_t cutc = hc->Integral( bin_lo, bin_hi);

	name = Form( "evgen/2cm/pi0_p_%d_in.root", eg[ebin]);
	evgen = new TFile( name);
	evhist = (TH1D*)evgen->Get( "h5");
	name = Form( "egp%d", eg[ebin]);
	TH1D *evcl = (TH1D*)evhist->Clone( name);
	evcl->Scale( 1/fact[ebin]);

	bin_lo = evcl->GetXaxis()->FindBin( th_lo);
	bin_hi = evcl->GetXaxis()->FindBin( th_hi)-1;
	Double_t totp = evcl->Integral( bin_lo, bin_hi);

	Double_t effp = cutp/totp;
	Double_t cont = cutc/(cutp + cutc);

	TH1D *hs = (TH1D*)hp->Clone();
	hs->Sumw2();
	hs->Add( hc, 1.0);
	max = hs->GetMaximum();
	max *= 1.4;
	hs->SetMaximum( max);

	if ( tbin == 0) name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	else name = "";
	hs->SetTitle( name);
	hs->GetXaxis()->SetTitle( "M_{miss} (MeV)");
	a = Linear( 145, 900, 400, 600, eg[ebin]);
	b = Linear( 145, 950, 400, 1000, eg[ebin]);

	hs->Draw();
	hs->GetXaxis()->SetRangeUser( a, b);
	hp->Draw( "same");
	hc->Draw( "same");

	l1 = new TLine( msep, 0, msep, max);
	l1->SetLineColor( 1);
	l1->SetLineStyle( 2);
	l1->SetLineWidth( 4);
	l1->Draw();

	pl = new TPaveText( 0.15, 0.60, 0.45, 0.85, "NDC");
	pl->SetBorderSize( 0);
	pl->SetFillStyle( 0);
	pl->SetTextAlign( 12);
	pl->SetTextSize( 0.04);
	name = Form( "#theta_{#gamma#gamma} = %d - 180 deg\n", (int) TggCut);
	pl->AddText( name);
	name = Form( "#theta^{#pi}_{cm} = %d - %d deg\n", th_lo, th_hi);
	pl->AddText( name);
	name = Form( "proton eff = %5.3f\n", effp);
	pl->AddText( name);
	name = Form( "carbon ratio = %5.3f\n", cont);
	pl->AddText( name);
	pl->Draw();

	name = Form( " egamma = %3d  theta = %3d  effp = %5.3f  cont = %5.3f",
			eg[ebin], th_mid, effp, cont);
	cout << name << endl;

	if ( tbin == 0) {
		tl = new TLegend( 0.15, 0.2, 0.45, 0.4);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.04);

		tl->AddEntry( hc, "^{12}C");
		tl->AddEntry( hp, "p");
		tl->AddEntry( hs, "sum", "P");

		tl->Draw();
	}

	TString name = Form( "plots/Aron/MmissTbin_%d_%d.pdf", eg[ebin], th_mid);
	if ( save == kTRUE) c1->Print( name);
}

void Mmiss2D( Int_t bin = 0, Bool_t save = kFALSE)
{
	Double_t Tgg, Tgg_c12, tgg_cut, theta_cut, msep;
	TString name, proj;

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 1050,
			700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 3,2);

	Tgg = Tgg_Min( qT_max(eg[ebin], kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	Tgg_c12 = Tgg_Min( qT_max(eg[ebin], kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;

	if ( bin == 0) {
		tgg_cut = Tgg + 10;
		theta_cut = 90;
		msep = 934;
	}
	else if ( bin == 1) {
		tgg_cut = Tgg + 10;
		theta_cut = 85;
		msep = 934;
	}
	else if ( bin == 2) {
		tgg_cut = Tgg + 8;
		theta_cut = 80;
		msep = 932;
	}
	else if ( bin == 3) {
		tgg_cut = Tgg + 8;
		theta_cut = 80;
		msep = 930;
	}
	else if ( bin == 4) {
		tgg_cut = Tgg + 8;
		theta_cut = 75;
		msep = 930;
	}

	cout << "CUTS:" << endl;
	cout << "  tgg: " << tgg_cut;
	cout << "  theta: " << theta_cut;
	cout << "  msep: " << msep;
	cout << endl;

	l1 = new TLine( 900, tgg_cut, 960, tgg_cut);
	l1->SetLineColor( 2);
	l1->SetLineStyle( 2);
	l1->SetLineWidth( 4);
	l2 = new TLine( 900, Tgg_c12, 960, Tgg_c12);
	l2->SetLineColor( 4);
	l2->SetLineStyle( 2);
	l2->SetLineWidth( 4);
	l3 = new TLine( 900, theta_cut, 960, theta_cut);
	l3->SetLineColor( 4);
	l3->SetLineStyle( 2);
	l3->SetLineWidth( 4);
	l4 = new TLine( msep, 0, msep, 180);
	l4->SetLineColor( 1);
	l4->SetLineStyle( 2);
	l4->SetLineWidth( 4);

// Gamma-Gamma Opening Angle vs. Missing Mass

	// Proton
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = "MMissP_v_TGGP_v_ThetaCMP";
	h3p = (TH3D*) prot->Get( name);
	h3p->SetTitle( "Proton");
	h3p->GetZaxis()->SetRangeUser( 900, 960);
	h3p->GetZaxis()->SetTitle( "M_{miss} (MeV)");
	proj = "yz";
	TH2D *h2p_tgg = h3p->Project3D( proj);
	h2p_tgg->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2p_tgg->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	name = Form( "E_{#gamma} = %3d MeV, Opening Angle Proton", eg[bin]);
	h2p_tgg->SetTitle( name);
	c1->cd( 1);
	h2p_tgg->Draw( "col");
	l1->Draw();

	// Carbon
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);
	name = "MMissP_v_TGGP_v_ThetaCMP";
	h3c = (TH3D*) carb->Get( name);
	h3c->SetTitle( "Carbon");
	h3c->GetZaxis()->SetRangeUser( 900, 960);
	h3c->GetZaxis()->SetTitle( "M_{miss} (MeV)");
	proj = "yz2";
	TH2D *h2c_tgg = h3c->Project3D( proj);
	h2c_tgg->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2c_tgg->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h2c_tgg->SetTitle( "Opening Angle Carbon");
	c1->cd( 2);
	h2c_tgg->Draw( "col");
//	l2->Draw();
	l1->Draw();

	// Sum
	TH2D *h2tgg_sum = (TH2D*)h2c_tgg->Clone();
	h2tgg_sum->Add(h2p_tgg, 1/fact[bin]);
	c1->cd( 3);
	h2tgg_sum->Draw( "col");
	h2tgg_sum->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2tgg_sum->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h2tgg_sum->SetTitle( "Opening Angle Sum");
	l1->Draw();
//	l2->Draw();

// Pion Angle vs. Missing Mass
	proj = "xz";

	// Proton
	TH2D *h2p_th = h3p->Project3D( proj);
	c1->cd( 4);
	h2p_th->Draw( "col");
	h2p_th->GetYaxis()->SetTitle( "#theta^{cm}_{#pi} (deg)");
	h2p_th->SetTitle( "Pion Angle Proton");
	l3->Draw();
	l4->Draw();

	// Carbon
	proj.Append( "2");
	TH2D *h2c_th = h3c->Project3D( proj);
	c1->cd(5);
	h2c_th->Draw( "col");
	h2c_th->GetYaxis()->SetTitle( "#theta^{cm}_{#pi} (deg)");
	h2c_th->SetTitle( "Pion Angle Carbon");
	l3->Draw();
	l4->Draw();

	// Sum
	TH2D *h2th_sum = (TH2D*)h2c_th->Clone();
	h2th_sum->Add(h2p_th, 1/fact[bin]);
	c1->cd(6);
	h2th_sum->Draw( "col");
	h2th_sum->GetYaxis()->SetTitle( "#theta^{cm}_{#pi} (deg)");
	h2th_sum->SetTitle( "Pion Angle Sum");
	l3->Draw();
	l4->Draw();

	name = Form( "plots/Aron/Mmiss2D_%d.pdf", eg[bin]);
	if ( save == kTRUE) c1->Print( name);
}

void TGG2D( Int_t bin = 0)
{
	Double_t Tgg, Tgg_c12, tgg_cut, theta_cut, msep;
	Int_t a[3], b[3];

	TString name, proj;

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 800, 0, 1000,
			500);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	name = "MMiss_TGG_ThetaCM_NC_P";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

// Gamma-Gamma Opening Angle vs. CM Theta
	proj = "yx";

	Tgg = Tgg_Min( qT_max(eg[bin], kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	Tgg_c12 = Tgg_Min( qT_max(eg[bin], kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;

	a[0] = 0;
	a[2] = 180;
	if ( bin == 0) {
		b[0] = 134;
		a[1] = 92.5;
		b[1] = 144;
		b[2] = 164;
	}
	else if ( bin == 1) {
		b[0] = 110;
		a[1] = 67.5;
		b[1] = 114;
		b[2] = 150;
	}
	else if ( bin == 2) {
		b[0] = 96;
		a[1] = 60;
		b[1] = 102;
		b[2] = 130;
	}

	l1 = new TLine( 20, Tgg_c12-10, 20, 180);
	l1->SetLineColor( 4);
	l1->SetLineStyle( 2);
	l1->SetLineWidth( 4);
	l2 = new TLine( 40, Tgg_c12-10, 40, 180);
	l2->SetLineColor( 4);
	l2->SetLineStyle( 2);
	l2->SetLineWidth( 4);
	l3 = new TLine( 60, Tgg_c12-10, 60, 180);
	l3->SetLineColor( 4);
	l3->SetLineStyle( 2);
	l3->SetLineWidth( 4);
	l4 = new TLine( 80, Tgg_c12-10, 80, 180);
	l4->SetLineColor( 4);
	l4->SetLineStyle( 2);
	l4->SetLineWidth( 4);
	l5 = new TLine( 100, Tgg_c12-10, 100, 180);
	l5->SetLineColor( 4);
	l5->SetLineStyle( 2);
	l5->SetLineWidth( 4);
	l6 = new TLine( 120, Tgg_c12-10, 120, 180);
	l6->SetLineColor( 4);
	l6->SetLineStyle( 2);
	l6->SetLineWidth( 4);
	l7 = new TLine( 140, Tgg_c12-10, 140, 180);
	l7->SetLineColor( 4);
	l7->SetLineStyle( 2);
	l7->SetLineWidth( 4);
	l8 = new TLine( 160, Tgg_c12-10, 160, 180);
	l8->SetLineColor( 4);
	l8->SetLineStyle( 2);
	l8->SetLineWidth( 4);


/*
	l1 = new TLine( a[0], b[0], a[1], b[1]);
	l1->SetLineColor( 4);
	l1->SetLineStyle( 1);
	l1->SetLineWidth( 4);
	l2 = new TLine( a[1], b[1], a[2], b[2]);
	l2->SetLineColor( 4);
	l2->SetLineStyle( 1);
	l2->SetLineWidth( 4);
*/

	// Proton
	TH2D *h2p_tgg = h3p->Project3D( proj);
	c1->cd( 1);
	h2p_tgg->Draw( "col");
	h2p_tgg->GetXaxis()->SetTitle( "#theta_{#pi}^{CM} (deg)");
	h2p_tgg->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2p_tgg->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	name = Form( "E_{#gamma} = %3d MeV, Proton", eg[bin]);
	h2p_tgg->SetTitle( name);
	l1->Draw();
	l2->Draw();
	l3->Draw();
	l4->Draw();
	l5->Draw();
	l6->Draw();
	l7->Draw();
	l8->Draw();

	// Carbon
	proj.Append( "2");
	TH2D *h2c_tgg = h3c->Project3D( proj);
	c1->cd( 2);
	h2c_tgg->Draw( "col");
	h2c_tgg->GetXaxis()->SetTitle( "#theta_{#pi}^{CM} (deg)");
	h2c_tgg->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2c_tgg->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h2c_tgg->SetTitle( "Carbon");
	l1->Draw();
	l2->Draw();
	l3->Draw();
	l4->Draw();
	l5->Draw();
	l6->Draw();
	l7->Draw();
	l8->Draw();

	name = Form( "plots/Aron/TGG_2D_%d.pdf", eg[bin]);
	c1->Print( name);
}

void TGG2DTheta( Int_t bin = 0, Int_t tbin = 0)
{
	Double_t Tgg, Tgg_c12, tgg_cut, theta_cut, msep;
	Int_t th_lo, th_hi;
	Int_t a[3], b[3];

	TString name, proj;

	th_lo = tbin*20;
	th_hi = 20 + tbin*20;

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.35, "xyz");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 800, 0, 1000,
			500);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	name = "MMiss_TGG_ThetaCM_NC_P";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

// Gamma-Gamma Opening Angle vs. CM Theta
	proj = "yx";

	Tgg = Tgg_Min( qT_max(eg[bin], kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	Tgg_c12 = Tgg_Min( qT_max(eg[bin], kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;

	a[0] = 0;
	a[2] = 180;
	if ( bin == 0) {
		b[0] = 134;
		a[1] = 92.5;
		b[1] = 144;
		b[2] = 164;
	}
	else if ( bin == 1) {
		b[0] = 110;
		a[1] = 67.5;
		b[1] = 114;
		b[2] = 150;
	}
	else if ( bin == 2) {
		b[0] = 96;
		a[1] = 60;
		b[1] = 102;
		b[2] = 130;
	}

	l1 = new TLine( a[0], b[0], a[1], b[1]);
	l1->SetLineColor( 4);
	l1->SetLineStyle( 1);
	l1->SetLineWidth( 4);
	l2 = new TLine( a[1], b[1], a[2], b[2]);
	l2->SetLineColor( 4);
	l2->SetLineStyle( 1);
	l2->SetLineWidth( 4);

	// Proton
	TH2D *h2p_tgg = h3p->Project3D( proj);
	c1->cd( 1);
	h2p_tgg->Draw( "col");
	h2p_tgg->GetXaxis()->SetTitle( "#theta_{#pi}^{CM} (deg)");
	h2p_tgg->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2p_tgg->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	name = Form( "E_{#gamma} = %3d MeV, Proton", eg[bin]);
	h2p_tgg->SetTitle( name);
//	l1->Draw();
//	l2->Draw();

	// Carbon
	proj.Append( "2");
	TH2D *h2c_tgg = h3c->Project3D( proj);
	c1->cd( 2);
	h2c_tgg->Draw( "col");
	h2c_tgg->GetXaxis()->SetTitle( "#theta_{#pi}^{CM} (deg)");
	h2c_tgg->GetYaxis()->SetRangeUser(Tgg_c12-10, 180);
	h2c_tgg->GetYaxis()->SetTitle( "#theta_{#gamma#gamma} (deg)");
	h2c_tgg->SetTitle( "Carbon");
//	l1->Draw();
//	l2->Draw();

//	name = Form( "plots/Aron/TGG_2D_%d.pdf", eg[bin]);
//	c1->Print( name);
}

void TGGDiff( Int_t bin = 0)
{
	Double_t a, b;
	TString name;

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 800, 0, 700, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.25, "xyz");

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/p_%d_2cm.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/c_%d_2cm.root", eg[bin]);
	carb = new TFile( name);

	hp = (TH1D*) prot->Get( "TGGDiff");
	hc = (TH1D*) carb->Get( "TGGDiff");

	name = Form( "E_{#gamma} = %3d MeV, 2-cm CH_{2} target", eg[bin]);
	hp->SetTitle( name);
	hp->GetXaxis()->SetTitle( "#theta_{#gamma#gamma}^{sim}-"
		"#theta_{#gamma#gamma}^{calc} (deg)");

	cout << " 2 cm target" << endl;
	cout << "scaling factor = " << 1/fact[bin] << endl;

	hp->SetLineWidth( 2);
	hp->SetFillColor( 2);
//	hp->Scale( 1/fact[bin]);
	hp->Draw();

	hc->SetFillStyle( 3001);
	hc->SetFillColor( 4);
	hc->SetLineWidth( 2);
	hc->SetLineStyle( 1);
	hc->Draw( "same");

	tl = new TLegend( 0.15, 0.6, 0.45, 0.8);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.04);
	tl->AddEntry( hp, "p");
	tl->AddEntry( hc, "^{12}C");
	tl->Draw();

//	TH1D *hsum = (TH1D*)hp->Clone();
//	hsum->Sumw2();
//	hsum->Add(hc, 1.0);
//	hsum->Draw( "same");

	name = Form( "plots/Aron/TGGDiff_%d.pdf", eg[bin]);
	c1->Print( name);
}

void qTDiff( Int_t bin = 0)
{
	TString name;

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 800, 0, 700, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.25, "xyz");

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/p_%d_2cm.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/c_%d_2cm.root", eg[bin]);
	carb = new TFile( name);

	hp = (TH1D*) prot->Get( "qTDiff");
	hc = (TH1D*) carb->Get( "qTDiff");

	name = Form( "E_{#gamma} = %3d MeV, 2-cm CH_{2} target", eg[bin]);
	hp->SetTitle( name);
	hp->GetXaxis()->SetTitle( "T_{#pi}^{sim}-T_{#pi}^{calc} (MeV)");

	cout << " 2 cm target" << endl;
	cout << "scaling factor = " << 1/fact[bin] << endl;

	hp->SetLineWidth( 2);
	hp->SetFillColor( 2);
//	hp->Scale( 1/fact[bin]);
	hp->Draw();

	hc->SetFillStyle( 3001);
	hc->SetFillColor( 4);
	hc->SetLineWidth( 2);
	hc->SetLineStyle( 1);
	hc->Draw( "same");

	tl = new TLegend( 0.15, 0.6, 0.45, 0.8);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.04);
	tl->AddEntry( hp, "p");
	tl->AddEntry( hc, "^{12}C");
	tl->Draw();

//	TH1D *hsum = (TH1D*)hp->Clone();
//	hsum->Sumw2();
//	hsum->Add(hc, 1.0);
//	hsum->Draw( "same");

	name = Form( "plots/Aron/qTDiff_%d.pdf", eg[bin]);
	c1->Print( name);
}

/*
void TGG2D( Int_t bin = 0)
{
	Double_t a, b;
	TString name;

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 800, 0, 700, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.25, "xyz");

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/p_%d_2cm.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/c_%d_2cm.root", eg[bin]);
	carb = new TFile( name);

	hp = (TH2D*) prot->Get( "TGGDiff_LTheta");
	hc = (TH2D*) carb->Get( "TGGDiff_LTheta");

	name = Form( "E_{#gamma} = %3d MeV, 2-cm CH_{2} target", eg[bin]);
	hp->SetTitle( name);
	hp->GetXaxis()->SetTitle( "#theta^{lab}_{#pi} (deg)");
	hp->GetYaxis()->SetTitle( "#theta_{#gamma#gamma}^{sim}-"
		"#theta_{#gamma#gamma}^{calc} (deg)");

	cout << " 2 cm target" << endl;
	cout << "scaling factor = " << 1/fact[bin] << endl;

//	hp->SetLineWidth( 2);
	hp->SetMarkerColor( 2);
//	hp->SetFillColor( 2);
//	hp->Scale( 1/fact[bin]);
	hp->Draw();

//	hc->SetFillStyle( 3001);
//	hc->SetFillColor( 4);
	hc->SetMarkerColor( 4);
//	hc->SetLineWidth( 2);
//	hc->SetLineStyle( 1);
	hc->Draw( "same");

//	tl = new TLegend( 0.75, 0.7, 0.95, 0.9);
//	tl->SetFillStyle( 0);
//	tl->SetBorderSize( 0);
//	tl->SetTextSize( 0.04);
//	tl->AddEntry( hp, "p", "f");
//	tl->AddEntry( hc, "^{12}C", "f");
//	tl->Draw();

//	TH1D *hsum = (TH1D*)hp->Clone();
//	hsum->Sumw2();
//	hsum->Add(hc, 1.0);
//	hsum->Draw( "same");

	name = Form( "plots/Aron/TGGDiff2D_%d.pdf", eg[bin]);
	c1->Print( name);
}
*/

void qT2D( Int_t bin = 0)
{
	Double_t a, b;
	TString name;

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 800, 0, 700, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	gStyle->SetOptStat( 0);
	gStyle->SetLabelSize( 0.025, "xyz");
	gStyle->SetTitleOffset( 1.25, "xyz");

	name = Form( "histograms/MonteCarlo/pi0/PolTarg/p_%d_2cm.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MonteCarlo/pi0/PolTarg/c_%d_2cm.root", eg[bin]);
	carb = new TFile( name);

	hp = (TH2D*) prot->Get( "qTDiff_LTheta");
	hc = (TH2D*) carb->Get( "qTDiff_LTheta");

	name = Form( "E_{#gamma} = %3d MeV, 2-cm CH_{2} target", eg[bin]);
	hp->SetTitle( name);
	hp->GetXaxis()->SetTitle( "#theta^{lab}_{#pi} (deg)");
	hp->GetYaxis()->SetTitle( "T_{#pi}^{sim}-T_{#pi}^{calc} (MeV)");

	cout << " 2 cm target" << endl;
	cout << "scaling factor = " << 1/fact[bin] << endl;

//	hp->SetLineWidth( 2);
	hp->SetMarkerColor( 2);
//	hp->SetFillColor( 2);
//	hp->Scale( 1/fact[bin]);
	hp->Draw();

//	hc->SetFillStyle( 3001);
//	hc->SetFillColor( 4);
	hc->SetMarkerColor( 4);
//	hc->SetLineWidth( 2);
//	hc->SetLineStyle( 1);
	hc->Draw( "same");

//	tl = new TLegend( 0.75, 0.7, 0.95, 0.9);
//	tl->SetFillStyle( 0);
//	tl->SetBorderSize( 0);
//	tl->SetTextSize( 0.04);
//	tl->AddEntry( hp, "p", "f");
//	tl->AddEntry( hc, "^{12}C", "f");
//	tl->Draw();

//	TH1D *hsum = (TH1D*)hp->Clone();
//	hsum->Sumw2();
//	hsum->Add(hc, 1.0);
//	hsum->Draw( "same");

	name = Form( "plots/Aron/qTDiff2D_%d.pdf", eg[bin]);
	c1->Print( name);
}
