/*********************************************************************
*                                                                    *
*                       ThetaEff.C                                   *
*                                                                    *
*  This ROOT macro calculates the CM theta efficiency from the       *
*  event generator input distribution and the GEANT output run       *
*  through the AcquRoot code.                                        *
*                                                                    *
*  Revision date                    Reason                           *
*  ----------------------------------------------------------------- *
*   16.09.05    Dave Hornidge       first version                    *
*   10.09.06    Dave Hornidge       updated version                  *
*   15.09.06    Dave Hornidge       including the new 3D plots       *
*   23.01.11    Dave Hornidge       adding cos(theta)						*
*		               	                                             *
**********************************************************************/

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void EffRange( UInt_t start, UInt_t end, TString experiment = "pi0",
		TString encl = "LH2_10cm", TString proc = "pi0_p", UInt_t rebin = 20)
{
	UInt_t i;

	for ( i = start; i <= end; ThetaEff( i++, experiment, encl, proc, kFALSE,
				rebin));
}

void ThetaEff( UInt_t chan = 300, TString expt = "pi0", TString tgt = "LH2",
		TString thick = "10cm", TString proc = "pi0_p", Bool_t plot = kTRUE,
		UInt_t rebin = 20)
{
	UInt_t i, end;
	Int_t yld_ev, yld_ge;
	Double_t eff, deff;
	TString encl;
	TString base, name;

	encl = tgt + "_" + thick;

	base = Form( "%s_chan%d", proc.Data(), chan);

	name = "evgen/" + thick + "/" + base + "_in.root";
	evgen = new TFile( name);

	name = "histograms/MonteCarlo/" + expt + "/" + encl + "/" + base +
		".root";
	geant = new TFile( name);

	name = "h5";
	ev_hist = (TH1D*) evgen->Get( name);
	TH1D *ev_clone = (TH1D*)ev_hist->Clone( "evgen");

//	name = "THR_ThetaCMP";
	name = "THR_ThetaCMCut1P";
//	name = "THR_ThetaCMCut2P";
	h1 = (TH1D*) geant->Get( name);
	TH1D *geant_clone = (TH1D*)h1->Clone( "geant");

	yld_ge = h1->Integral();
	yld_ev = ev_hist->Integral();
	eff = (double) yld_ge/yld_ev;
	deff = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));

	ev_hist->Rebin( rebin);
	h1->Rebin( rebin);
	h1->Sumw2();
	h1->Divide( ev_hist);

	name = Form( "xs/pi0/eff/DetEff_chan_%d.out", chan);
	ofstream outFile( name);

	name = Form( "  0 %6.4f %6.4f", eff, deff);
	outFile << name << endl;
	cout << name << endl;

	end = 180/rebin;

	for ( i = 1; i <= end; i++) {
		name = Form( "%3.0f %6.4f %6.4f", h1->GetBinCenter( i),
				h1->GetBinContent( i), h1->GetBinError( i));
		cout << name << endl;
		outFile << name << endl;
	}
	outFile.close();

	if ( plot == kTRUE) {
		c1 = new TCanvas( "canvas", "Set-Up Efficiency", 700, 700);
		c1->Divide(1,2);

		c1->cd(1);
		ev_clone->Draw();
		ev_clone->SetFillColor( 46);
		ev_clone->SetLineWidth( 2);
		ev_clone->SetMinimum( 0);
		geant_clone->Draw( "same");
		geant_clone->SetFillColor( 38);
		geant_clone->SetLineWidth( 2);

		c1->cd(2);
		h1->Draw( "E");
		h1->SetLineWidth( 2);
		h1->SetMinimum( 0);
		h1->SetMaximum( 1);
	}
}

//
// Detection Efficiency as a function of Cos(Theta)
//
// Note: currently set up for 20 bins.
//
void CThetaEff( UInt_t chan = 300, TString expt = "pi0", TString tgt = "LH2",
		TString thick = "10cm", TString proc = "pi0_p", Bool_t plot = kTRUE)
{
	UInt_t i, end;
	Int_t yld_ev, yld_ge;
	Double_t eff, deff;
	TString encl;
	TString base, name;

	encl = tgt + "_" + thick;

	base = Form( "%s_chan%d", proc.Data(), chan);

	name = "evgen/" + thick + "/" + base + "_in.root";
	evgen = new TFile( name);

	name = "histograms/MonteCarlo/" + expt + "/" + encl + "/" + base +
		".root";
	geant = new TFile( name);

	name = "h16";
	ev_hist = (TH1D*) evgen->Get( name);
	TH1D *ev_clone = (TH1D*)ev_hist->Clone( "evgen");

	name = "THR_CosThetaCMCut3P";
	h1 = (TH1D*) geant->Get( name);
	TH1D *geant_clone = (TH1D*)h1->Clone( "geant");

	yld_ge = h1->Integral();
	yld_ev = ev_hist->Integral();
	eff = (double) yld_ge/yld_ev;
	deff = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));

	h1->Sumw2();
	h1->Divide( ev_hist);

	name = Form( "xs/pi0/eff/DetEff_cth_chan_%d.out", chan);
	ofstream outFile( name);

	name = Form( "total %6.4f %6.4f", eff, deff);
	outFile << name << endl;
	cout << name << endl;

	for ( i = 1; i <= 20; i++) {
		name = Form( "%5.2f %6.4f %6.4f", h1->GetBinCenter( i),
				h1->GetBinContent( i), h1->GetBinError( i));
		cout << name << endl;
		outFile << name << endl;
	}
	outFile.close();

	if ( plot == kTRUE) {
		c1 = new TCanvas( "canvas", "Set-Up Efficiency", 700, 700);
		c1->Divide(1,2);

		c1->cd(1);
		ev_clone->Draw();
		ev_clone->SetFillColor( 46);
		ev_clone->SetLineWidth( 2);
		ev_clone->SetMinimum( 0);
		geant_clone->Draw( "same");
		geant_clone->SetFillColor( 38);
		geant_clone->SetLineWidth( 2);

		c1->cd(2);
		h1->Draw( "E");
		h1->SetLineWidth( 2);
		h1->SetMinimum( 0);
		h1->SetMaximum( 1);
	}
}

void ThetaEffCompAll( Int_t e_lo = 200, Int_t e_hi = 300)
{
	for ( Int_t eg = e_lo; eg <= e_hi; eg += 10) ThetaEffComp( eg);
}

void ThetaEffComp( Int_t eg, TString tgt = "p")
{
	Int_t i, bins, yld_ev, yld_ge, lo, hi, theta_lo, theta_hi, binsize;
	Double_t max, t_eff, dt_eff;
	Double_t th[9], dth[9], eff[9], deff[9];
	TString name;

	c1 = new TCanvas( "canvas", "Set-Up Efficiency", 200, 0, 1000, 800);
	c1->Divide(2,2);
	gStyle->SetOptStat(0);

	bins = 9;
	binsize = 180/bins;

	name = Form( "evgen/comp_%s_%d_in.root", tgt.Data(), eg);
	evgen = new TFile( name);
	name = Form( "histograms/MCH_comp_%s_%d.root", tgt.Data(), eg);
	geant = new TFile( name);

// Event Generator

	name = "h2";
	ev_h = (TH1D*) evgen->Get( name);
	ev_h->SetFillStyle(3001);
	ev_h->SetFillColor(18);
	ev_h->SetLineWidth(2);

// Gamma ONLY

	name = "MMissPhot_Theta_TChan_P";
	ge_h3 = (TH3D*) geant->Get( name);
	TH1D *ge_h = ge_h3->Project3D("y");
	ge_h->SetFillColor(9);
	ge_h->SetLineWidth(2);

	name = "Photon Only";
	c1->cd(1);
	if ( ev_h->GetMaximum() > ge_h->GetMaximum()) {
		max = ev_h->GetMaximum();
		max *= 1.30;
		ev_h->SetMaximum(max);
		ev_h->SetTitle( name);
		ev_h->Draw();
		ev_h->GetXaxis()->SetTitle("#theta_{CM} (deg)");

		ge_h->Draw("same");
	}
	else {
		max = ge_h->GetMaximum();
		max *= 1.30;
		ge_h->SetMaximum(max);
		ge_h->SetTitle( name);
		ge_h->Draw();
		ge_h->GetXaxis()->SetTitle("#theta_{CM} (deg)");

		ev_h->Draw("same");
	}
	tl = new TLegend( 0.6, 0.7, 0.8, 0.85);
	tl->SetFillStyle(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.04);
	tl->AddEntry( ev_h, "Event Generator");
	tl->AddEntry( ge_h, "Analysed");
	tl->Draw();

	name = Form( "E_{#gamma} = %3d MeV\n", eg);
	pl = new TPaveLabel( 0.4, 0.9, 0.6, 0.99, name, "NDC");
	pl->SetBorderSize(0);
	pl->SetFillStyle(0);
//	pl->SetTextAlign(12);
	pl->SetTextSize(0.4);
	pl->Draw();

	yld_ge = ge_h->Integral();
	yld_ev = ev_h->Integral();
	t_eff = (double) yld_ge/yld_ev;
	dt_eff = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));

	cout << " tot ";
	cout << yld_ge;
	cout << " ";
	cout << yld_ev;
	cout << " ";
	cout << t_eff;
	cout << endl;

	Double_t evc, gec;
	evc = 0;
	gec = 0;

	for ( i = 0; i < bins; i++) {
		theta_lo = i*binsize;
		theta_hi = (i+1)*binsize;

		lo = ev_h->GetXaxis()->FindBin( theta_lo);
		hi = ev_h->GetXaxis()->FindBin( theta_hi)-1;
		yld_ev = ev_h->Integral( lo, hi);

		evc += yld_ev;

		lo = ge_h->GetXaxis()->FindBin( theta_lo);
		hi = ge_h->GetXaxis()->FindBin( theta_hi)-1;
		yld_ge = ge_h->Integral( lo, hi);

		gec += yld_ge;

		th[i] = (theta_lo+theta_hi)/2;
		dth[i] = 0;
		eff[i] = (double) yld_ge/yld_ev;
		deff[i] = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));
		cout << th[i];
		cout << "  ";
		cout << yld_ev;
		cout << "  ";
		cout << yld_ge;
		cout << "  ";
		cout << eff[i];
		cout << endl;
	}
//	cout << " check ";
//	cout << gec;
//	cout << " ";
//	cout << evc;
//	cout << " ";
//	cout << gec/evc;
//	cout << endl;

/*
	// Output g only efficiencies
	name = Form( "efficiency/Eff_comp_g_%s_%d.out", tgt.Data(), eg);
	ofstream outFile( name);
	outFile << " 0 ";
	outFile.width(7);
	outFile << yld_ev << " ";
	outFile.width(7);
	outFile << yld_ge << " ";
	outFile.width(10);
	outFile.precision(4);
	outFile << eff << " ";
	outFile.precision(2);
	outFile.width(10);
	outFile << deff << endl;

	for ( i = 1; i <= bins; i++) {
		outFile.width(2);
		outFile << i << " ";
		outFile.width(7);
		outFile << yld_ev << " ";
		outFile.width(7);
		outFile << yld_ge << " ";
		outFile.width(10);
		outFile.precision(4);
		outFile << eff[i] << " ";
		outFile.precision(2);
		outFile.width(10);
		outFile << deff[i] << endl;
	}
	outFile.close();
*/

// Gamma AND Proton

	c1->cd(2);
	gr = new TGraphErrors( bins, th, eff, dth, deff);
	gr->SetTitle();
	gr->SetMarkerStyle( 20);
	gr->SetMarkerSize( 1.0);
	gr->SetMarkerColor( 4);
	gr->GetXaxis()->SetTitle("#theta_{CM} (deg)");
	gr->GetYaxis()->SetTitle("#varepsilon_{#gamma'}");
	gr->GetXaxis()->SetRangeUser(0,180);
	gr->GetYaxis()->SetRangeUser(0,1.4);
	gr->Draw("AP");

	l1 = new TLine(0,1,180,1);
	l1->SetLineStyle(2);
	l1->SetLineWidth(1);
	l1->Draw();

	name = "MMissPhot_Theta_TChan_CutP";
	ge_h3cut = (TH3D*) geant->Get( name);
	TH1D *ge_h2 = ge_h3cut->Project3D("y");
	ge_h2->SetFillColor(9);
	ge_h2->SetLineWidth(2);

	name = "Photon and Proton";
	c1->cd(3);
	if ( ev_h->GetMaximum() > ge_h2->GetMaximum()) {
		max = ev_h->GetMaximum();
		max *= 1.30;
		ev_h->SetMaximum(max);
		ev_h->SetTitle(name);
		ev_h->GetXaxis()->SetTitle("#theta_{CM} (deg)");
		ev_h->Draw();

		ge_h2->Draw("same");
	}
	else {
		max = ge_h2->GetMaximum();
		max *= 1.30;
		ge_h2->SetMaximum(max);
		ge_h2->SetTitle(name);
		ge_h2->GetXaxis()->SetTitle("#theta_{CM} (deg)");
		ge_h2->Draw();

		ev_h->Draw("same");
	}

	yld_ge = ge_h2->Integral();
	yld_ev = ev_h->Integral();
	t_eff = (double) yld_ge/yld_ev;
	dt_eff = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));

	cout << endl;
	cout << " tot ";
	cout << yld_ge;
	cout << " ";
	cout << yld_ev;
	cout << " ";
	cout << t_eff;
	cout << endl;

	evc = 0;
	gec = 0;
	for ( i = 0; i < bins; i++) {
		theta_lo = i*binsize;
		theta_hi = (i+1)*binsize;

		lo = ev_h->GetXaxis()->FindBin( theta_lo);
		hi = ev_h->GetXaxis()->FindBin( theta_hi)-1;
		yld_ev = ev_h->Integral( lo, hi);

		evc += yld_ev;

		lo = ge_h2->GetXaxis()->FindBin( theta_lo);
		hi = ge_h2->GetXaxis()->FindBin( theta_hi)-1;
		yld_ge = ge_h2->Integral( lo, hi);

		gec += yld_ge;

		th[i] = (theta_lo+theta_hi)/2;
		dth[i] = 0;
		eff[i] = (double) yld_ge/yld_ev;
		deff[i] = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));
		cout << th[i];
		cout << "  ";
		cout << yld_ev;
		cout << "  ";
		cout << yld_ge;
		cout << "  ";
		cout << eff[i];
		cout << endl;
	}
//	cout << " check ";
//	cout << gec;
//	cout << " ";
//	cout << evc;
//	cout << " ";
//	cout << gec/evc;
//	cout << endl;

	c1->cd(4);
	gr2 = new TGraphErrors( bins, th, eff, dth, deff);
	gr2->SetTitle();
	gr2->SetMarkerStyle( 20);
	gr2->SetMarkerSize( 1.0);
	gr2->SetMarkerColor( 4);
	gr2->GetXaxis()->SetTitle("#theta_{CM} (deg)");
	gr2->GetYaxis()->SetTitle("#varepsilon_{#gamma'}");
	gr2->GetXaxis()->SetRangeUser(0,180);
	gr2->GetYaxis()->SetRangeUser(0,1.4);
	gr2->Draw("AP");

	l1->Draw();

	// Output g + p efficiencies
//	name = Form( "efficiency/Eff_comp_gp_%s_%d.out", tgt.Data(), eg);
	name = Form( "efficiency/eff_%d.dat", eg);
	ofstream outFile( name);

	for ( i = 0; i < bins; i++) {
		outFile.width(2);
		outFile << 10+i*20;
		outFile << " ";
		outFile.precision(4);
		outFile << eff[i];
		outFile << " ";
//		outFile.precision(2);
//		outFile.width(10);
//		outFile << deff[i];
		outFile << endl;
	}
	outFile.close();

	name = Form( "efficiency/comp_teff_%d.pdf", eg);
	c1->Print( name);
}

void ThetaEffEff( Bool_t plot = kTRUE, TString tgt = "p", Int_t eg = 100)
{
	Int_t rebin;
	Int_t yld_ev, yld_ge;
	Double_t eff, deff;
	TString base, name;

	rebin = 10;

	base = Form( "eff_%s_%d", tgt.Data(), eg);

	name = "evgen/" + base + "_in.root";
	evgen = new TFile( name);

	name = "histograms/MonteCarlo/" + base + ".root";
	geant = new TFile( name);

	name = "h2";
	ev_hist = (TH1D*) evgen->Get( name);
	TH1D *ev_clone = (TH1D*)ev_hist->Clone( "evgen_clone");

	if ( tgt == "p") name = "EFF_ProtTheta";
	else if ( tgt == "g") name = "EFF_PhotTheta";
	h1 = (TH1D*) geant->Get( name);
	TH1D *geant_clone = (TH1D*)h1->Clone( "geant_clone");

	yld_ge = h1->Integral();
	yld_ev = ev_hist->Integral();
	eff = (double) yld_ge/yld_ev;
	deff = sqrt( yld_ge/Sqr( yld_ev) + Sqr( yld_ge)/pow( yld_ev, 3));

	ev_hist->Rebin( rebin);
	h1->Rebin( rebin);
	h1->Sumw2();
	h1->Divide( ev_hist);

	name = Form( "xs/eff/det_eff_new/Eff_%s_%d.out", tgt.Data(), eg);
	ofstream outFile( name);

	name = Form( "  0 %6.4f %6.4f", eff, deff);
	outFile << name << endl;
	cout << name << endl;

	if ( plot == kTRUE) {
		c1 = new TCanvas( "canvas", "Set-Up Efficiency", 700, 700);
		c1->Divide(1,2);

		c1->cd(1);
		ev_clone->Draw();
		ev_clone->SetFillColor( 46);
		ev_clone->SetLineWidth( 2);
		geant_clone->Draw( "same");
		geant_clone->SetFillColor( 38);
		geant_clone->SetLineWidth( 2);

		c1->cd(2);
		h1->Draw( "E");
		h1->SetLineWidth( 2);
	}
}
