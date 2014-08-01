/*********************************************************************
*																							*
*						CheckComp.C														*
*																							*
*	This ROOT macro...		 														*
*																							*
*	Revision date                    Reason                           *
*  ----------------------------------------------------------------- *
*	2009.05.07		Dave Hornidge		first version							*
*																							*
**********************************************************************/

gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

Int_t eg[] = { 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300};
Int_t th[] = { 0, 10, 30, 50, 70, 90, 110, 130, 150, 170};

Double_t effp, dil_eff;

// The standard parameters are:
// 	LH2:
// 		240 MeV:		4	945	10		kTRUE
// 		280 MeV:		8	940	10		kTRUE
// 	Butanol:
// 		240 MeV:		4	945	5		kFALSE
// 		280 MeV:		8	940	5		kFALSE

void CheckEbin( Int_t ebin = 4, Double_t msep = 945, Double_t open_cut = 10,
		Bool_t lin = kTRUE, Bool_t plot = kFALSE)
{
	Int_t tbin;
	TString name;

	if ( lin == kTRUE) name = Form( "efficiency/eff_LH2_%d_%d_%d.out",
			(int) open_cut, (int) msep, eg[ebin]);
	else name = Form( "efficiency/eff_butanol_%d_%d_%d.out",
			(int) open_cut, (int) msep, eg[ebin]);

	ofstream outFile( name);
	if ( !outFile.is_open())
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit( -1);
	}

	outFile << " 10 0.000 0.000";
	outFile << endl;
	outFile << " 30 0.000 0.000";
	outFile << endl;
	outFile << " 50 0.000 0.000";
	outFile << endl;
	for ( tbin = 4; tbin < 9; tbin++) {
//		CheckComp( ebin, tbin, "pi0_p", open_cut);
		AddBG( ebin, tbin, msep, open_cut, lin, plot);
		name = Form( "%3d %5.3f %5.3f", th[tbin], effp, dil_eff);
		outFile << name;
		outFile << endl;
	}
	outFile << "170 0.000 0.000";
	outFile << endl;

	outFile.close();

}

void CheckComp( Int_t ebin = 8, Int_t tbin = 0, TString reaction = "pi0_p",
		Double_t open_cut = 10)
{
	Int_t i, bin_lo, bin_hi, pbin;
	Double_t fact[100][100];
	Double_t left, right, open_cut;
	Double_t x[2], y[2];
	Double_t max, msep;
	Double_t pc, bc, pctot, evp, cont;
	Int_t junk;
	TString name, title;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 500, 0, 700, 500);

	gStyle->SetOptStat( 0);

	name = Form( "histograms/MonteCarlo/compton/LH2_10cm/compton_p_%d.root",
			eg[ebin]);
	comp = new TFile( name);

	name = Form( "histograms/MonteCarlo/compton/LH2_10cm/%s_%d.root",
			reaction.Data(), eg[ebin]);
	back = new TFile( name);

	i = 0;
	name = "includes/txs_ratio.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) 
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit(1);
	}
	for ( i = 0; i < 11; i++)
		inFile >> junk >> fact[i][0] >> fact[i][1] >> fact[i][2]
			>> fact[i][3] >> fact[i][4] >> fact[i][5] >> fact[i][6];
	inFile.close();

	title = "#pi^{0} p";
	pbin = 0;
	if ( reaction == "comp_qfc12") {
		title = "QF #gamma from ^{12}C";
		pbin = 1;
	}
	else if ( reaction == "pi0_qfc12") {
		title = "QF #pi^{0} from ^{12}C";
		pbin = 2;
	}
	else if ( reaction == "comp_qfhe3") {
		title = "QF #gamma from ^{3}He";
		pbin = 3;
	}
	else if ( reaction == "comp_qfhe4") {
		title = "QF #gamma from ^{4}He";
		pbin = 4;
	}
	else if ( reaction == "pi0_qfhe3") {
		title = "QF #pi^{0} from ^{3}He";
		pbin = 5;
	}
	else if ( reaction == "pi0_qfhe4") {
		title = "QF #pi^{0} from ^{4}He";
		pbin = 6;
	}

	// Default background process is pi0_p
	cout << fact[ebin][pbin] << endl;

	left = 900;
	right = 1050;

	msep = 940;
	if ( eg[ebin] == 240) msep = 945;

	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";

	// Compton
	h3comp = (TH3D*) comp->Get( name);
	h3comp->GetXaxis()->SetRangeUser( 0, open_cut);
	if ( tbin == 0) {
		bin_lo = h3comp->GetYaxis()->FindBin( 0);
		bin_hi = h3comp->GetYaxis()->FindBin( 180);
	}
	else {
		bin_lo = h3comp->GetYaxis()->FindBin( th[tbin] - 10);
		bin_hi = h3comp->GetYaxis()->FindBin( th[tbin] + 10) - 1;
	}
	h3comp->GetYaxis()->SetRange( bin_lo, bin_hi);
	h3comp->GetZaxis()->SetRangeUser( left, right);
	hcomp = (TH1D*) h3comp->Project3D( "z");

	// Background
	h3back = (TH3D*) back->Get( name);
	h3back->GetXaxis()->SetRangeUser( 0, open_cut);
	h3back->GetYaxis()->SetRange( bin_lo, bin_hi);
	h3back->GetZaxis()->SetRangeUser( left, right);
	hback = (TH1D*) h3back->Project3D( "z2");

	// Event Generator Information
	name = Form( "evgen/10cm/compton_p_%d_in.root", eg[ebin]);
	evgen_p = new TFile( name);
	name = "h5";
	ev_p = (TH1D*) evgen_p->Get( name);
	if ( tbin == 0) {
		bin_lo = ev_p->GetXaxis()->FindBin( 0);
		bin_hi = ev_p->GetXaxis()->FindBin( 180);
	}
	else
	{
		bin_lo = ev_p->GetXaxis()->FindBin( th[tbin] - 10);
		bin_hi = ev_p->GetXaxis()->FindBin( th[tbin] + 10) - 1;
	}
	cout << bin_lo << " " << bin_hi << endl;
	evp = ev_p->Integral( bin_lo, bin_hi);

	hback->SetFillStyle(3001);
	hback->SetLineWidth(2);
	hback->SetFillColor(17);
//	hback->SetFillColor(4);
//	hback->SetLineColor(4);
	hback->Scale( fact[ebin][pbin]);

	hcomp->SetFillStyle(3008);
	hcomp->SetFillColor(13);
	hcomp->SetLineWidth(2);
//	hcomp->SetFillStyle(3001);
//	hcomp->SetFillColor(2);
//	hcomp->SetLineColor(2);
//	hcomp->Scale( fact[ebin][pbin]);

	name = "";
	if ( hcomp->GetMaximum() > hback->GetMaximum()) {
		max = hcomp->GetMaximum();
		max *= 1.10;
		hcomp->SetMaximum( max);
		hcomp->SetTitle( name);
		hcomp->GetXaxis()->SetTitle( "M_{miss} (MeV)");
		hcomp->GetXaxis()->CenterTitle();
		hcomp->Draw();
		hback->Draw( "same");
	}
	else {
		max = hback->GetMaximum();
		max *= 1.10;
		hback->SetMaximum( max);
		hback->SetTitle( name);
		hback->GetXaxis()->SetTitle( "M_{miss} (MeV)");
		hback->GetXaxis()->CenterTitle();
		hback->Draw();
		hcomp->Draw( "same");
	}

	x[0] = 0.15;
	x[1] = 0.30;
	if ( eg[ebin] == 240) {
		x[0] = 0.5;
		x[1] = 0.65;
	}
	y[0] = 0.7;
	y[1] = 0.9;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.04);
	tl->SetTextFont(132);
	tl->AddEntry( hcomp, "#gamma p");
	tl->AddEntry( hback, title);
	tl->Draw();

//	name = Form( "E_{#gamma} = %3d #pm 1 MeV", eg[ebin]);
	name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	pl = new TPaveLabel( 0.15, 0.87, 0.40, 1.00, name, "NDC");
	pl->SetFillStyle(0);
	pl->SetBorderSize(0);
	pl->SetTextSize(0.4);
	pl->SetTextFont(132);
	pl->Draw();

	y[0] = 0.75;
	if ( tbin != 0) {
		y[0] = 0.60;
	}
	y[1] = 0.87;
	x[0] = 0.6;
	x[1] = 0.9;
	pt = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextSize(0.04);
	pt->SetTextFont(132);
	if ( tbin != 0) {
		name = Form( "#theta^{CM}_{#gamma} = %d #pm 10 deg", th[tbin]);
		pt->AddText( name);
	}
	name = Form( "#theta_{open} < %d deg", (int) open_cut);
	pt->AddText( name);
	name = Form( "M_{miss} < %d MeV", (int) msep);
	pt->AddText( name);
	pt->Draw();

	bin_lo = hcomp->FindBin( 900);
	bin_hi = hcomp->FindBin( msep);
	pc = hcomp->Integral( bin_lo, bin_hi);
	bc = hback->Integral( bin_lo, bin_hi);
	pctot = hcomp->Integral();

	effp = pc/evp;
	if ( pc != 0) cont = bc/pc;
	else cont = 0;

	cout << evp;
	cout << " ";
	cout << pctot;
	cout << " ";
	cout << pc;
	cout << " ";
	cout << bc;
	cout << " ";
	cout << cont;
	cout << " ";
	cout << effp;
	cout << endl;

	l1 = new TLine( msep, 0, msep, max);
	l1->SetLineStyle(3);
	l1->SetLineWidth(4);
	l1->Draw();

	name = Form( "efficiency/MMissCut%d_%s_%d_%d", (int) open_cut,
			reaction.Data(), eg[ebin], th[tbin]);
//	name.Append( ".eps");
	name.Append( ".pdf");
//	c1->Print( name);
}

void Plot2D( Int_t ebin = 8, TString proc = "compton", TString tgt = "p")
{
	Double_t left, right;
	TString name, title, reaction;
	TString s1, s2, s3;

	reaction = proc + "_" + tgt;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);

	gStyle->SetOptStat( 0);

//	name = Form( "histograms/MonteCarlo/Compton/LH2_10cm/%s_%d.root",
//			reaction.Data(), eg[ebin]);
	name = Form( "histograms/MonteCarlo/Compton/Poltarg/%s_%d.root",
			reaction.Data(), eg[ebin]);
	comp = new TFile( name);
	name = "COMP_PhotonMmissProtP_v_PhotProtOAP";
	hcomp = (TH2D*) comp->Get( name);

	Int_t counts = hcomp->Integral();
	cout << counts;
	cout << endl;

	title = "#gamma p #rightarrow ";
	if ( proc == "compton") s1 = "#gamma ";
	else if ( proc == "pi0") s1 = "#pi^{0} ";
	s3 = "QF ";
	if ( tgt == "p") {
		s2 = "p";
		s3 = "";
	}
	else if ( tgt = "qfc12") s2 = "^{12}C";
	else if ( tgt = "qfhe3") s2 = "^{3}He";
	else if ( tgt = "qfhe4") s2 = "^{4}He";
	name = s3 + s1 + s2;
	title.Append( name);

	hcomp->SetTitle( title);
	hcomp->SetFillStyle( 3001);
	hcomp->SetFillColor( 2);
	hcomp->SetLineWidth( 2);
	hcomp->SetLineColor( 2);
	hcomp->GetXaxis()->SetTitle( "#theta_{open} (deg)");
	hcomp->GetXaxis()->CenterTitle();
	hcomp->GetXaxis()->SetTitleSize( 0.04);
	hcomp->GetXaxis()->SetTitleOffset( 1.2);
	hcomp->GetXaxis()->SetLabelSize( 0.03);
	hcomp->GetYaxis()->SetTitle( "M_{miss} (MeV/c^{2})");
	hcomp->GetYaxis()->CenterTitle();
	hcomp->GetYaxis()->SetTitleSize( 0.04);
	hcomp->GetYaxis()->SetTitleOffset( 1.2);
	hcomp->GetYaxis()->SetLabelSize( 0.03);
	hcomp->Draw( "colz");

	name = Form( "E_{#gamma} = %3d MeV and Proton Detected",
		eg[ebin]);
	pl = new TPaveLabel( 0.5, 0.9, 0.9, 0.99, name, "NDC");
	pl->SetBorderSize(0);
	pl->SetFillStyle(0);
	pl->SetTextSize(0.4);
	pl->Draw();

	name = Form( "plots/Compton/Mmiss_2D_%s_%d", reaction.Data(), eg[ebin]);

//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);
}

void CheckComp2D( Int_t ebin = 8, TString proc = "comp_p",
		Double_t open_cut = 10)
{
	Int_t bin;
	TString name, title;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 400, 0, 700, 500);

	gStyle->SetOptStat( 0);

	name = Form( "histograms/MonteCarlo/%s_%d.root", proc.Data(), eg[ebin]);
	comp = new TFile( name);

	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";
	h3comp = (TH3D*) comp->Get( name);
	h3comp->GetXaxis()->SetRangeUser( 0, open_cut);
	h3comp->GetYaxis()->SetRangeUser( 0, 180);
	h3comp->GetZaxis()->SetRangeUser( 900, 1000);

	h2comp = (TH2D*) h3comp->Project3D( "zy");

	name = "";
	h2comp->SetTitle( name);
	h2comp->GetXaxis()->SetTitle( "#theta_{#gamma}^{CM} (deg)");
	h2comp->GetYaxis()->SetTitle( "M_{miss} (MeV)");
	h2comp->Draw( "zcol");

	Int_t counts = h2comp->Integral();
	cout << counts;
	cout << endl;

//	name = Form( "plots/pdf/comp_%s_%d", proc.Data(), eg[ebin]);
//	name.Append( "_cut");
//	name.Append( ".eps");
//	name.Append( ".pdf");
//	c1->Print( name);
}

void AddBGAll()
{
	Int_t ebin, tbin;
	Double_t msep, open_cut;
	Bool_t lin, plot;

	tbin = 0;
	plot = kTRUE;

// 240 MeV
	ebin = 4;
	msep = 945;

	// Linear Polarization and Hydrogen Target
	open_cut = 10;
	lin = kTRUE;
	AddBG( ebin, tbin, msep, open_cut, lin, plot);

	// Circular Polarization and Frozen-Spin Target
	open_cut = 5;
	lin = kFALSE;
	AddBG( ebin, tbin, msep, open_cut, lin, plot);

// 280 MeV
	ebin = 8;
	msep = 940;

	// Linear Polarization and Hydrogen Target
	open_cut = 10;
	lin = kTRUE;
	AddBG( ebin, tbin, msep, open_cut, lin, plot);

	// Circular Polarization and Frozen-Spin Target
	open_cut = 5;
	lin = kFALSE;
	AddBG( ebin, tbin, msep, open_cut, lin, plot);
}

void AddBG( Int_t ebin = 8, Int_t tbin = 0, Double_t msep = 945,
		Double_t open_cut = 10, Bool_t lin = kTRUE, Bool_t plot = kFALSE,
		Bool_t save = kFALSE)
{
	Int_t i, j;
	Int_t rebin, pbin;
	Int_t ebin_lo, ebin_hi, tbin_lo, tbin_hi;
	Double_t fact[100][100];
	Double_t left, right, open_cut;
	Double_t x[2], y[2];
	Double_t max;
	Double_t pc, bc, pctot, evp;
	Double_t cont;
	Int_t junk;
	TString name, reaction;
	Bool_t CH2;

	CH2 = kTRUE;
//	CH2 = kFALSE;

	left = 850;
	right = 1050;

	rebin = 1;

	i = 0;
	name = "includes/txs_ratio.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) 
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit(1);
	}
	for ( i = 0; i < 11; i++)
	{
		inFile >> junk >> fact[i][0] >> fact[i][1] >> fact[i][2]
			>> fact[i][3] >> fact[i][4] >> fact[i][5] >> fact[i][6];
//		for ( j = 0; j <= 6; fact[i][j++] /= 10);
	}
	inFile.close();

	if ( plot == kTRUE) {

		if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
		TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 500, 0, 700,
				500);

		gStyle->SetOptStat( 0);
	}

	name = Form( "histograms/MonteCarlo/comp_p_%d.root", eg[ebin]);
	comp = new TFile( name);

	reaction = "pi0_p";
	name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
			eg[ebin]);
	back1 = new TFile( name);

	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";

	h3comp = (TH3D*) comp->Get( name);
	h3comp->GetXaxis()->SetRangeUser( 0, open_cut);
	if ( tbin == 0) {
		tbin_lo = h3comp->GetYaxis()->FindBin( 0);
		tbin_hi = h3comp->GetYaxis()->FindBin( 180);
	}
	else {
		tbin_lo = h3comp->GetYaxis()->FindBin( th[tbin] - 10);
		tbin_hi = h3comp->GetYaxis()->FindBin( th[tbin] + 10)-1;
	}
	ebin_lo = h3comp->GetZaxis()->FindBin( left);
	ebin_hi = h3comp->GetZaxis()->FindBin( right)-1;

	// Compton
	h3comp->GetYaxis()->SetRange( tbin_lo, tbin_hi);
	h3comp->GetZaxis()->SetRange( ebin_lo, ebin_hi);
	hcomp = (TH1D*) h3comp->Project3D( "z");

	// Pi0P Background - pbin 0
	pbin = 0;
	h3back1 = (TH3D*) back1->Get( name);
	h3back1->GetXaxis()->SetRangeUser( 0, open_cut);
	h3back1->GetYaxis()->SetRange( tbin_lo, tbin_hi);
	h3back1->GetZaxis()->SetRange( ebin_lo, ebin_hi);
	hback1 = (TH1D*) h3back1->Project3D( "z1");
	hback1->Scale( fact[ebin][pbin]);

	// Add them all together
	hback = (TH1D*) hback1->Clone( "all");

	if ( lin == kFALSE)
	{
		reaction = "comp_qfc12";
		name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
				eg[ebin]);
		back2 = new TFile( name);

		reaction = "pi0_qfc12";
		name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
				eg[ebin]);
		back3 = new TFile( name);

		if ( CH2 == kFALSE) {

			reaction = "comp_qfhe3";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back4 = new TFile( name);

			reaction = "comp_qfhe4";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back5 = new TFile( name);

			reaction = "pi0_qfhe3";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back6 = new TFile( name);

			reaction = "pi0_qfhe4";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back7 = new TFile( name);
		}

		name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";

		// Comp QF C12 Background - pbin 1
		pbin = 1;
		h3back2 = (TH3D*) back2->Get( name);
		h3back2->GetXaxis()->SetRangeUser( 0, open_cut);
		h3back2->GetYaxis()->SetRange( tbin_lo, tbin_hi);
		h3back2->GetZaxis()->SetRange( ebin_lo, ebin_hi);
		hback2 = (TH1D*) h3back2->Project3D( "z2");
		hback2->Scale( fact[ebin][pbin]);

		// Pi0 QF C12 Background - pbin 2
		pbin = 2;
		h3back3 = (TH3D*) back3->Get( name);
		h3back3->GetXaxis()->SetRangeUser( 0, open_cut);
		h3back3->GetYaxis()->SetRange( tbin_lo, tbin_hi);
		h3back3->GetZaxis()->SetRange( ebin_lo, ebin_hi);
		hback3 = (TH1D*) h3back3->Project3D( "z3");
		hback3->Scale( fact[ebin][pbin]);

		hback->Add( hback2, 1.0);
		hback->Add( hback3, 1.0);

		if ( CH2 == kFALSE) {

			// Comp QF He3 Background - pbin 3
			pbin = 3;
			h3back4 = (TH3D*) back4->Get( name);
			h3back4->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back4->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back4->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback4 = (TH1D*) h3back4->Project3D( "z4");
			hback4->Scale( fact[ebin][pbin]);

			// Comp QF He4 Background - pbin 4
			pbin = 4;
			h3back5 = (TH3D*) back5->Get( name);
			h3back5->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back5->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back5->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback5 = (TH1D*) h3back5->Project3D( "z5");
			hback5->Scale( fact[ebin][pbin]);

			// Pi0 QF He3 Background - pbin 5
			pbin = 5;
			h3back6 = (TH3D*) back6->Get( name);
			h3back6->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back6->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back6->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback6 = (TH1D*) h3back6->Project3D( "z6");
			hback6->Scale( fact[ebin][pbin]);

			// Pi0 QF He4 Background - pbin 6
			pbin = 6;
			h3back7 = (TH3D*) back7->Get( name);
			h3back7->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back7->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back7->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback7 = (TH1D*) h3back7->Project3D( "z7");
			hback7->Scale( fact[ebin][pbin]);

			hback->Add( hback4, 1.0);
			hback->Add( hback5, 1.0);
			hback->Add( hback6, 1.0);
			hback->Add( hback7, 1.0);
		}

	}

	// Event Generator Information
	name = Form( "evgen/comp_p_%d_in.root", eg[ebin]);
	evgen_p = new TFile( name);
	name = "h2";
	ev_p = (TH1D*) evgen_p->Get( name);
	evp = ev_p->Integral( tbin_lo, tbin_hi);

	hcomp->SetFillStyle(3004);
	hcomp->SetFillColor(16);
	hcomp->SetLineWidth(2);
	hcomp->SetLineStyle(1);

	hback->SetFillStyle(3005);
	hback->SetFillColor(16);
	hback->SetLineWidth(3);
	hback->SetLineStyle(2);

	hback->Rebin( rebin);

	name = "";
	if ( hcomp->GetMaximum() > hback->GetMaximum()) {
		max = hcomp->GetMaximum();
		max *= 1.20;
		hcomp->SetMaximum( max);
		hcomp->SetTitle( name);
		hcomp->GetXaxis()->SetTitle( "M_{miss} (MeV)");
		hcomp->GetXaxis()->CenterTitle();
		if ( plot == kTRUE) {
			hcomp->Draw();
			hback->Draw( "same");
		}
	}
	else {
		max = hback->GetMaximum();
		max *= 1.20;
		hback->SetMaximum( max);
		hback->SetTitle( name);
		hback->GetXaxis()->SetTitle( "M_{miss} (MeV)");
		hback->GetXaxis()->CenterTitle();
		if ( plot == kTRUE) {
			hback->Draw();
			hcomp->Draw( "same");
		}
	}

	x[0] = 0.15;
	x[1] = 0.4;
	y[0] = 0.68;
	y[1] = 0.88;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.04);
	tl->SetTextFont(132);
	tl->AddEntry( hcomp, "Compton");
	tl->AddEntry( hback, "Background");
	if ( plot == kTRUE) tl->Draw();

//	name = Form( "E_{#gamma} = %3d #pm 1 MeV", eg[ebin]);
	name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	pl = new TPaveLabel( 0.15, 0.87, 0.40, 1.00, name, "NDC");
	pl->SetFillStyle(0);
	pl->SetBorderSize(0);
	pl->SetTextSize(0.4);
	pl->SetTextFont(132);
	if ( plot == kTRUE) pl->Draw();
	if ( tbin != 0) {
		name = Form( "#theta^{CM}_{#gamma} = %d #pm 10 deg", th[tbin]);
		pl2 = new TPaveLabel( 0.60, 0.87, 0.85, 1.00, name, "NDC");
		pl2->SetFillStyle(0);
		pl2->SetBorderSize(0);
		pl2->SetTextSize(0.4);
		pl2->SetTextFont(132);
	}

	x[0] = 0.07;
	x[1] = 0.37;
	y[0] = 0.55;
	y[1] = 0.65;
	pt = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextSize(0.04);
	pt->SetTextFont(132);
	name = Form( "#theta_{open} < %d deg", (int) open_cut);
	pt->AddText( name);
	name = Form( "M_{miss} < %d MeV", (int) msep);
	pt->AddText( name);
	if ( plot == kTRUE) pt->Draw();

	ebin_lo = hcomp->FindBin( 900);
	ebin_hi = hcomp->FindBin( msep);
	pc = hcomp->Integral( ebin_lo, ebin_hi);
	bc = hback->Integral( ebin_lo, ebin_hi);
	pctot = hcomp->Integral();

	effp = pc/evp;
	if ( pc != 0) cont = bc/pc;
	else cont = 0;

	if ( (pc + bc) != 0) dil_eff = pc/(pc + bc);
	else dil_eff = 0;

	name = Form( "%3d %5.0f %5.0f %5.0f %5.0f %5.3f %5.3f %5.3f", th[tbin], evp,
			pctot, pc, bc, cont, effp, dil_eff);
	cout << name;
	cout << endl;

	l1 = new TLine( msep, 0, msep, max);
	l1->SetLineStyle(3);
	l1->SetLineWidth(4);

	if ( plot == kTRUE) {

		l1->Draw();

		if ( lin == kTRUE) 
			name = Form( "xs/efficiency/MMissCut%d_LH2_%d_%d", (int) open_cut,
					eg[ebin], th[tbin]);
		else if ( CH2 == kFALSE)
			name = Form( "xs/efficiency/MMissCut%d_butanol_%d_%d", (int) open_cut,
					eg[ebin], th[tbin]);
		else
			name = Form( "xs/efficiency/MMissCut%d_CH2_%d_%d", (int) open_cut,
					eg[ebin], th[tbin]);

		name.Append( ".pdf");
		if ( save == kTRUE) c1->Print( name);

	}
}

void AddBG_EMiss( Int_t ebin = 8, Int_t tbin = 0, Double_t msep = 20,
		Double_t open_cut = 10, Bool_t lin = kTRUE, Bool_t plot = kFALSE)
{
	Int_t i, j;
	Int_t rebin, pbin;
	Int_t ebin_lo, ebin_hi, tbin_lo, tbin_hi;
	Double_t fact[100][100];
	Double_t left, right, open_cut;
	Double_t x[2], y[2];
	Double_t max;
	Double_t pc, bc, pctot, evp;
	Double_t cont;
	Int_t junk;
	TString name, reaction;
	Bool_t CH2;

	CH2 = kTRUE;
//	CH2 = kFALSE;

	left = -50;
	right = 150;

	rebin = 1;

	i = 0;
	name = "includes/txs_ratio.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) 
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit(1);
	}
	for ( i = 0; i < 11; i++)
	{
		inFile >> junk >> fact[i][0] >> fact[i][1] >> fact[i][2]
			>> fact[i][3] >> fact[i][4] >> fact[i][5] >> fact[i][6];
//		for ( j = 0; j <= 6; fact[i][j++] /= 10);
	}
	inFile.close();

	if ( plot == kTRUE) {

		if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
		TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 500, 0, 700,
				500);

		gStyle->SetOptStat( 0);
	}

	name = Form( "histograms/MonteCarlo/comp_p_%d.root", eg[ebin]);
	comp = new TFile( name);

	reaction = "pi0_p";
	name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
			eg[ebin]);
	back1 = new TFile( name);

	name = "EMissPhot_Theta_OpenA_P";

	h3comp = (TH3D*) comp->Get( name);
	h3comp->GetXaxis()->SetRangeUser( 0, open_cut);
	if ( tbin == 0) {
		tbin_lo = h3comp->GetYaxis()->FindBin(0);
		tbin_hi = h3comp->GetYaxis()->FindBin(180);
	}
	else {
		tbin_lo = h3comp->GetYaxis()->FindBin(th[tbin] - 10);
		tbin_hi = h3comp->GetYaxis()->FindBin(th[tbin] + 10)-1;
	}
	ebin_lo = h3comp->GetZaxis()->FindBin( left);
	ebin_hi = h3comp->GetZaxis()->FindBin( right)-1;

	// Compton
	h3comp->GetYaxis()->SetRange( tbin_lo, tbin_hi);
	h3comp->GetZaxis()->SetRange( ebin_lo, ebin_hi);
	hcomp = (TH1D*) h3comp->Project3D( "z");

	// Pi0P Background - pbin 0
	pbin = 0;
	h3back1 = (TH3D*) back1->Get( name);
	h3back1->GetXaxis()->SetRangeUser( 0, open_cut);
	h3back1->GetYaxis()->SetRange( tbin_lo, tbin_hi);
	h3back1->GetZaxis()->SetRange( ebin_lo, ebin_hi);
	hback1 = (TH1D*) h3back1->Project3D( "z1");
	hback1->Scale( fact[ebin][pbin]);

	// Add them all together
	hback = (TH1D*) hback1->Clone( "all");

	if ( lin == kFALSE)
	{
		reaction = "comp_qfc12";
		name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
				eg[ebin]);
		back2 = new TFile( name);

		reaction = "pi0_qfc12";
		name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
				eg[ebin]);
		back3 = new TFile( name);

		if ( CH2 == kFALSE) {

			reaction = "comp_qfhe3";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back4 = new TFile( name);

			reaction = "comp_qfhe4";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back5 = new TFile( name);

			reaction = "pi0_qfhe3";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back6 = new TFile( name);

			reaction = "pi0_qfhe4";
			name = Form( "histograms/MonteCarlo/%s_%d.root",
					reaction.Data(), eg[ebin]);
			back7 = new TFile( name);
		}

		name = "EMissPhot_Theta_OpenA_P";

		// Comp QF C12 Background - pbin 1
		pbin = 1;
		h3back2 = (TH3D*) back2->Get( name);
		h3back2->GetXaxis()->SetRangeUser( 0, open_cut);
		h3back2->GetYaxis()->SetRange( tbin_lo, tbin_hi);
		h3back2->GetZaxis()->SetRange( ebin_lo, ebin_hi);
		hback2 = (TH1D*) h3back2->Project3D( "z2");
		hback2->Scale( fact[ebin][pbin]);

		// Pi0 QF C12 Background - pbin 2
		pbin = 2;
		h3back3 = (TH3D*) back3->Get( name);
		h3back3->GetXaxis()->SetRangeUser( 0, open_cut);
		h3back3->GetYaxis()->SetRange( tbin_lo, tbin_hi);
		h3back3->GetZaxis()->SetRange( ebin_lo, ebin_hi);
		hback3 = (TH1D*) h3back3->Project3D( "z3");
		hback3->Scale( fact[ebin][pbin]);

		hback->Add( hback2, 1.0);
		hback->Add( hback3, 1.0);

		if ( CH2 == kFALSE) {

			// Comp QF He3 Background - pbin 3
			pbin = 3;
			h3back4 = (TH3D*) back4->Get( name);
			h3back4->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back4->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back4->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback4 = (TH1D*) h3back4->Project3D( "z4");
			hback4->Scale( fact[ebin][pbin]);

			// Comp QF He4 Background - pbin 4
			pbin = 4;
			h3back5 = (TH3D*) back5->Get( name);
			h3back5->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back5->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back5->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback5 = (TH1D*) h3back5->Project3D( "z5");
			hback5->Scale( fact[ebin][pbin]);

			// Pi0 QF He3 Background - pbin 5
			pbin = 5;
			h3back6 = (TH3D*) back6->Get( name);
			h3back6->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back6->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back6->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback6 = (TH1D*) h3back6->Project3D( "z6");
			hback6->Scale( fact[ebin][pbin]);

			// Pi0 QF He4 Background - pbin 6
			pbin = 6;
			h3back7 = (TH3D*) back7->Get( name);
			h3back7->GetXaxis()->SetRangeUser( 0, open_cut);
			h3back7->GetYaxis()->SetRange( tbin_lo, tbin_hi);
			h3back7->GetZaxis()->SetRange( ebin_lo, ebin_hi);
			hback7 = (TH1D*) h3back7->Project3D( "z7");
			hback7->Scale( fact[ebin][pbin]);

			hback->Add( hback4, 1.0);
			hback->Add( hback5, 1.0);
			hback->Add( hback6, 1.0);
			hback->Add( hback7, 1.0);
		}

	}

	// Event Generator Information
	name = Form( "evgen/comp_p_%d_in.root", eg[ebin]);
	evgen_p = new TFile( name);
	name = "h2";
	ev_p = (TH1D*) evgen_p->Get( name);
	evp = ev_p->Integral( tbin_lo, tbin_hi);

	hcomp->SetFillStyle(3004);
	hcomp->SetFillColor(16);
	hcomp->SetLineWidth(2);
	hcomp->SetLineStyle(1);

	hback->SetFillStyle(3005);
	hback->SetFillColor(16);
	hback->SetLineWidth(3);
	hback->SetLineStyle(2);

	hback->Rebin( rebin);

	name = "";
	if ( hcomp->GetMaximum() > hback->GetMaximum()) {
		max = hcomp->GetMaximum();
		max *= 1.20;
		hcomp->SetMaximum( max);
		hcomp->SetTitle( name);
		hcomp->GetXaxis()->SetTitle( "E_{miss} (MeV)");
		hcomp->GetXaxis()->CenterTitle();
		if ( plot == kTRUE) {
			hcomp->Draw();
			hback->Draw( "same");
		}
	}
	else {
		max = hback->GetMaximum();
		max *= 1.20;
		hback->SetMaximum( max);
		hback->SetTitle( name);
		hback->GetXaxis()->SetTitle( "E_{miss} (MeV)");
		hback->GetXaxis()->CenterTitle();
		if ( plot == kTRUE) {
			hback->Draw();
			hcomp->Draw( "same");
		}
	}

	x[0] = 0.15;
	x[1] = 0.4;
	y[0] = 0.68;
	y[1] = 0.88;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.04);
	tl->SetTextFont(132);
	tl->AddEntry( hcomp, "Compton");
	tl->AddEntry( hback, "Background");
	if ( plot == kTRUE) tl->Draw();

//	name = Form( "E_{#gamma} = %3d #pm 1 MeV", eg[ebin]);
	name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	pl = new TPaveLabel( 0.15, 0.87, 0.40, 1.00, name, "NDC");
	pl->SetFillStyle(0);
	pl->SetBorderSize(0);
	pl->SetTextSize(0.4);
	pl->SetTextFont(132);
	if ( plot == kTRUE) pl->Draw();
	if ( tbin != 0) {
		name = Form( "#theta^{CM}_{#gamma} = %d #pm 10 deg", th[tbin]);
		pl2 = new TPaveLabel( 0.60, 0.87, 0.85, 1.00, name, "NDC");
		pl2->SetFillStyle(0);
		pl2->SetBorderSize(0);
		pl2->SetTextSize(0.4);
		pl2->SetTextFont(132);
	}

	x[0] = 0.07;
	x[1] = 0.37;
	y[0] = 0.55;
	y[1] = 0.65;
	pt = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextSize(0.04);
	pt->SetTextFont(132);
	name = Form( "#theta_{open} < %d deg", (int) open_cut);
	pt->AddText( name);
	name = Form( "E_{miss} < %d MeV", (int) msep);
	pt->AddText( name);
	if ( plot == kTRUE) pt->Draw();

	ebin_lo = hcomp->FindBin( 0);
	ebin_hi = hcomp->FindBin( msep);
	pc = hcomp->Integral( ebin_lo, ebin_hi);
	bc = hback->Integral( ebin_lo, ebin_hi);
	pctot = hcomp->Integral();

	effp = pc/evp;
	if ( pc != 0) cont = bc/pc;
	else cont = 0;

	if ( (pc + bc) != 0) dil_eff = pc/(pc + bc);
	else dil_eff = 0;

	name = Form( "%3d %5.0f %5.0f %5.0f %5.0f %5.3f %5.3f %5.3f", th[tbin], evp,
			pctot, pc, bc, cont, effp, dil_eff);
	cout << name;
	cout << endl;

	l1 = new TLine( msep, 0, msep, max);
	l1->SetLineStyle(3);
	l1->SetLineWidth(4);

	if ( plot == kTRUE) {

		l1->Draw();

		if ( lin == kTRUE) 
			name = Form( "efficiency/EMissCut%d_LH2_%d_%d", (int) open_cut,
					eg[ebin], th[tbin]);
		else if ( CH2 == kFALSE)
			name = Form( "efficiency/EMissCut%d_butanol_%d_%d", (int) open_cut,
					eg[ebin], th[tbin]);
		else
			name = Form( "efficiency/EMissCut%d_CH2_%d_%d", (int) open_cut,
					eg[ebin], th[tbin]);

		name.Append( ".eps");
	//	name.Append( ".pdf");
		c1->Print( name);

	}
}

void CheckCompSum( Int_t ebin = 8, Int_t tbin = 0, TString reaction = "pi0_p",
		Double_t open_cut = 10, Bool_t save = kFALSE)
{
	Int_t i, bin_lo, bin_hi, pbin;
	Double_t fact[100][100];
	Double_t left, right, open_cut;
	Double_t x[2], y[2];
	Double_t max, msep;
	Double_t pc, bc, pctot, evp, cont;
	Int_t junk;
	TString name, title;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 500, 0, 700, 500);

	gStyle->SetOptStat( 0);

	name = Form( "histograms/MonteCarlo/comp_p_%d.root", eg[ebin]);
	comp = new TFile( name);

	name = Form( "histograms/MonteCarlo/%s_%d.root", reaction.Data(),
			eg[ebin]);
	back = new TFile( name);

	i = 0;
	name = "includes/txs_ratio.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) 
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit(1);
	}
	for ( i = 0; i < 11; i++)
		inFile >> junk >> fact[i][0] >> fact[i][1] >> fact[i][2]
			>> fact[i][3] >> fact[i][4] >> fact[i][5] >> fact[i][6];
	inFile.close();

	title = "#pi^{0} p";
	pbin = 0;
	if ( reaction == "comp_qfc12") {
		title = "QF #gamma from ^{12}C";
		pbin = 1;
	}
	else if ( reaction == "pi0_qfc12") {
		title = "QF #pi^{0} from ^{12}C";
		pbin = 2;
	}
	else if ( reaction == "comp_qfhe3") {
		title = "QF #gamma from ^{3}He";
		pbin = 3;
	}
	else if ( reaction == "comp_qfhe4") {
		title = "QF #gamma from ^{4}He";
		pbin = 4;
	}
	else if ( reaction == "pi0_qfhe3") {
		title = "QF #pi^{0} from ^{3}He";
		pbin = 5;
	}
	else if ( reaction == "pi0_qfhe4") {
		title = "QF #pi^{0} from ^{4}He";
		pbin = 6;
	}

	// Default background process is pi0_p
//	if ( ( ebin == 4) || ( ebin == 8)) fact[ebin][pbin] /= 10; // 1e6 events
	cout << fact[ebin][pbin] << endl;


//	left = 900;
//	right = 1050;
	left = 800;
	right = 1100;

	msep = 940;
	if ( eg[ebin] == 240) msep = 945;

	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";

	// Compton
	h3comp = (TH3D*) comp->Get( name);
	h3comp->GetXaxis()->SetRangeUser( 0, open_cut);
	if ( tbin == 0) {
		bin_lo = h3comp->GetYaxis()->FindBin(0);
		bin_hi = h3comp->GetYaxis()->FindBin(180);
	}
	else {
		bin_lo = h3comp->GetYaxis()->FindBin(th[tbin] - 10);
		bin_hi = h3comp->GetYaxis()->FindBin(th[tbin] + 10)-1;
	}
	h3comp->GetYaxis()->SetRange( bin_lo, bin_hi);
	h3comp->GetZaxis()->SetRangeUser( left, right);
	hcomp = (TH1D*) h3comp->Project3D( "z");

	// Background
	h3back = (TH3D*) back->Get( name);
	h3back->GetXaxis()->SetRangeUser( 0, open_cut);
	h3back->GetYaxis()->SetRange( bin_lo, bin_hi);
	h3back->GetZaxis()->SetRangeUser( left, right);
	hback = (TH1D*) h3back->Project3D( "z2");

	// Event Generator Information
	name = Form( "evgen/comp_p_%d_in.root", eg[ebin]);
	evgen_p = new TFile( name);
	name = "h2";
	ev_p = (TH1D*) evgen_p->Get( name);
	evp = ev_p->Integral( bin_lo, bin_hi);

	hback->SetFillStyle(3001);
	hback->SetLineWidth(2);
	hback->SetFillColor(17);
//	hback->SetFillColor(4);
//	hback->SetLineColor(4);
	hback->Scale( fact[ebin][pbin]);

	hcomp->SetFillStyle( 3008);
	hcomp->SetFillColor( 13);
	hcomp->SetLineWidth( 2);
//	hcomp->SetFillStyle( 3001);
//	hcomp->SetFillColor( 2);
//	hcomp->SetLineColor( 2);
//	hcomp->Scale( fact[ebin][pbin]);

	hsum = (TH1D*) hcomp->Clone( "sum");
	hsum->Add( hback, 1.0);

	name = "";
	max = hsum->GetMaximum();
	max *= 1.4;
	hsum->SetFillStyle( 0);
	hsum->SetLineWidth( 4);
	hsum->SetLineStyle( 2);
	hsum->SetMaximum( max);
	hsum->SetTitle( name);
	hsum->GetXaxis()->SetTitle( "M_{miss} (MeV)");
	hsum->GetXaxis()->CenterTitle();

	hsum->Draw();
	hcomp->Draw( "same");
	hback->Draw( "same");

	x[0] = 0.15;
	x[1] = 0.30;
	if ( eg[ebin] == 240) {
		x[0] = 0.5;
		x[1] = 0.65;
	}
	y[0] = 0.7;
	y[1] = 0.9;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.04);
	tl->SetTextFont( 132);
	tl->AddEntry( hcomp, "#gamma p");
	tl->AddEntry( hback, title);
	tl->AddEntry( hsum, "Sum");
	tl->Draw();

	name = Form( "E_{#gamma} = %3d #pm 1 MeV", eg[ebin]);
//	name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	pl = new TPaveLabel( 0.15, 0.87, 0.40, 1.00, name, "NDC");
	pl->SetFillStyle(0);
	pl->SetBorderSize(0);
	pl->SetTextSize(0.4);
	pl->SetTextFont(132);
	pl->Draw();

	y[0] = 0.75;
//	if ( tbin != 0) {
//		y[0] = 0.60;
//	}
	y[1] = 0.87;
	x[0] = 0.6;
	x[1] = 0.9;
	pt = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pt->SetBorderSize( 0);
	pt->SetFillStyle( 0);
	pt->SetTextSize( 0.04);
	pt->SetTextFont( 132);
	pt->SetTextAlign( 12);
	if ( tbin != 0) {
		name = Form( "#theta^{CM}_{#gamma} = %d #pm 10 deg", th[tbin]);
		pt->AddText( name);
	}
	name = Form( "#theta_{open} < %d deg", (int) open_cut);
	pt->AddText( name);
//	name = Form( "M_{miss} < %d MeV", (int) msep);
//	pt->AddText( name);
	pt->Draw();

	bin_lo = hcomp->FindBin( 900);
	bin_hi = hcomp->FindBin( msep);
	pc = hcomp->Integral( bin_lo, bin_hi);
	bc = hback->Integral( bin_lo, bin_hi);
	pctot = hcomp->Integral();

	effp = pc/evp;
	if ( pc != 0) cont = bc/pc;
	else cont = 0;

	cout << evp;
	cout << " ";
	cout << pctot;
	cout << " ";
	cout << pc;
	cout << " ";
	cout << bc;
	cout << " ";
	cout << cont;
	cout << " ";
	cout << effp;
	cout << endl;

//	l1 = new TLine( msep, 0, msep, max);
//	l1->SetLineStyle(3);
//	l1->SetLineWidth(4);
//	l1->Draw();

	name = Form( "plots/PhotonMmissSim_%d_%d_%d", (int) open_cut, eg[ebin],
			th[tbin]);
	name.Append( ".pdf");
	if ( save == kTRUE) c1->Print( name);
}

void Plot2D_new()
{
	UInt_t ebin;
	Double_t left, right;
	TString proc, tgt;
	TString name, title, reaction;
	TString s1, s2, s3;

	ebin = 10;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 700);
	c1->Divide( 2, 2);

	gStyle->SetOptStat( 0);

	// Compton from proton
	c1->cd( 1);
	proc = "compton";
	tgt = "p";

	reaction = proc + "_" + tgt;

	name = Form( "histograms/MonteCarlo/Compton/Poltarg/%s_%d.root",
			reaction.Data(), eg[ebin]);
	comp = new TFile( name);
	name = "COMP_PhotonMmissProtP_v_PhotProtOAP";
	hcomp = (TH2D*) comp->Get( name);

	Int_t counts = hcomp->Integral();
	cout << counts;
	cout << endl;

	title = "#gamma p #rightarrow ";
	if ( proc == "compton") s1 = "#gamma ";
	else if ( proc == "pi0") s1 = "#pi^{0} ";
	s3 = "QF ";
	if ( tgt == "p") {
		s2 = "p";
		s3 = "";
	}
	else if ( tgt = "qfc12") s2 = "^{12}C";
	else if ( tgt = "qfhe3") s2 = "^{3}He";
	else if ( tgt = "qfhe4") s2 = "^{4}He";
	name = s3 + s1 + s2;
	title.Append( name);

	hcomp->SetTitle( title);
	hcomp->SetFillStyle( 3001);
	hcomp->SetFillColor( 2);
	hcomp->SetLineWidth( 2);
	hcomp->SetLineColor( 2);
	hcomp->GetXaxis()->SetTitle( "#theta_{open} (deg)");
	hcomp->GetXaxis()->CenterTitle();
	hcomp->GetXaxis()->SetTitleSize( 0.04);
	hcomp->GetXaxis()->SetTitleOffset( 1.2);
	hcomp->GetXaxis()->SetLabelSize( 0.03);
	hcomp->GetYaxis()->SetTitle( "M_{miss} (MeV/c^{2})");
	hcomp->GetYaxis()->CenterTitle();
	hcomp->GetYaxis()->SetTitleSize( 0.04);
	hcomp->GetYaxis()->SetTitleOffset( 1.2);
	hcomp->GetYaxis()->SetLabelSize( 0.03);
	hcomp->Draw( "colz");
	hcomp->GetXaxis()->SetRangeUser( 0, 40);

	l1 = new TLine( 10, 800, 10, 1100);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 4);
	l1->Draw();

	// Pi0 from proton
	c1->cd( 2);
	proc = "pi0";
	tgt = "p";

	reaction = proc + "_" + tgt;

	name = Form( "histograms/MonteCarlo/Compton/Poltarg/%s_%d.root",
			reaction.Data(), eg[ebin]);
	comp = new TFile( name);
	name = "COMP_PhotonMmissProtP_v_PhotProtOAP";
	hcomp = (TH2D*) comp->Get( name);

	Int_t counts = hcomp->Integral();
	cout << counts;
	cout << endl;

	title = "#gamma p #rightarrow ";
	if ( proc == "compton") s1 = "#gamma ";
	else if ( proc == "pi0") s1 = "#pi^{0} ";
	s3 = "QF ";
	if ( tgt == "p") {
		s2 = "p";
		s3 = "";
	}
	else if ( tgt = "qfc12") s2 = "^{12}C";
	else if ( tgt = "qfhe3") s2 = "^{3}He";
	else if ( tgt = "qfhe4") s2 = "^{4}He";
	name = s3 + s1 + s2;
	title.Append( name);

	hcomp->SetTitle( title);
	hcomp->SetFillStyle( 3001);
	hcomp->SetFillColor( 2);
	hcomp->SetLineWidth( 2);
	hcomp->SetLineColor( 2);
	hcomp->GetXaxis()->SetTitle( "#theta_{open} (deg)");
	hcomp->GetXaxis()->CenterTitle();
	hcomp->GetXaxis()->SetTitleSize( 0.04);
	hcomp->GetXaxis()->SetTitleOffset( 1.2);
	hcomp->GetXaxis()->SetLabelSize( 0.03);
	hcomp->GetYaxis()->SetTitle( "M_{miss} (MeV/c^{2})");
	hcomp->GetYaxis()->CenterTitle();
	hcomp->GetYaxis()->SetTitleSize( 0.04);
	hcomp->GetYaxis()->SetTitleOffset( 1.2);
	hcomp->GetYaxis()->SetLabelSize( 0.03);
	hcomp->Draw( "colz");
	hcomp->GetXaxis()->SetRangeUser( 0, 40);

	l1->Draw();

	// QF Compton from carbon
	c1->cd( 3);
	proc = "compton";
	tgt = "qfc12";

	reaction = proc + "_" + tgt;

	name = Form( "histograms/MonteCarlo/Compton/Poltarg/%s_%d.root",
			reaction.Data(), eg[ebin]);
	comp = new TFile( name);
	name = "COMP_PhotonMmissProtP_v_PhotProtOAP";
	hcomp = (TH2D*) comp->Get( name);

	Int_t counts = hcomp->Integral();
	cout << counts;
	cout << endl;

	title = "#gamma p #rightarrow ";
	if ( proc == "compton") s1 = "#gamma ";
	else if ( proc == "pi0") s1 = "#pi^{0} ";
	s3 = "QF ";
	if ( tgt == "p") {
		s2 = "p";
		s3 = "";
	}
	else if ( tgt = "qfc12") s2 = "^{12}C";
	else if ( tgt = "qfhe3") s2 = "^{3}He";
	else if ( tgt = "qfhe4") s2 = "^{4}He";
	name = s3 + s1 + s2;
	title.Append( name);

	hcomp->SetTitle( title);
	hcomp->SetFillStyle( 3001);
	hcomp->SetFillColor( 2);
	hcomp->SetLineWidth( 2);
	hcomp->SetLineColor( 2);
	hcomp->GetXaxis()->SetTitle( "#theta_{open} (deg)");
	hcomp->GetXaxis()->CenterTitle();
	hcomp->GetXaxis()->SetTitleSize( 0.04);
	hcomp->GetXaxis()->SetTitleOffset( 1.2);
	hcomp->GetXaxis()->SetLabelSize( 0.03);
	hcomp->GetYaxis()->SetTitle( "M_{miss} (MeV/c^{2})");
	hcomp->GetYaxis()->CenterTitle();
	hcomp->GetYaxis()->SetTitleSize( 0.04);
	hcomp->GetYaxis()->SetTitleOffset( 1.2);
	hcomp->GetYaxis()->SetLabelSize( 0.03);
	hcomp->Draw( "colz");
	hcomp->GetXaxis()->SetRangeUser( 0, 40);

	l1->Draw();

	// QF Pi0 from carbon
	c1->cd( 4);
	proc = "pi0";
	tgt = "qfc12";

	reaction = proc + "_" + tgt;

	name = Form( "histograms/MonteCarlo/Compton/Poltarg/%s_%d.root",
			reaction.Data(), eg[ebin]);
	comp = new TFile( name);
	name = "COMP_PhotonMmissProtP_v_PhotProtOAP";
	hcomp = (TH2D*) comp->Get( name);

	Int_t counts = hcomp->Integral();
	cout << counts;
	cout << endl;

	title = "#gamma p #rightarrow ";
	if ( proc == "compton") s1 = "#gamma ";
	else if ( proc == "pi0") s1 = "#pi^{0} ";
	s3 = "QF ";
	if ( tgt == "p") {
		s2 = "p";
		s3 = "";
	}
	else if ( tgt = "qfc12") s2 = "^{12}C";
	else if ( tgt = "qfhe3") s2 = "^{3}He";
	else if ( tgt = "qfhe4") s2 = "^{4}He";
	name = s3 + s1 + s2;
	title.Append( name);

	hcomp->SetTitle( title);
	hcomp->SetFillStyle( 3001);
	hcomp->SetFillColor( 2);
	hcomp->SetLineWidth( 2);
	hcomp->SetLineColor( 2);
	hcomp->GetXaxis()->SetTitle( "#theta_{open} (deg)");
	hcomp->GetXaxis()->CenterTitle();
	hcomp->GetXaxis()->SetTitleSize( 0.04);
	hcomp->GetXaxis()->SetTitleOffset( 1.2);
	hcomp->GetXaxis()->SetLabelSize( 0.03);
	hcomp->GetYaxis()->SetTitle( "M_{miss} (MeV/c^{2})");
	hcomp->GetYaxis()->CenterTitle();
	hcomp->GetYaxis()->SetTitleSize( 0.04);
	hcomp->GetYaxis()->SetTitleOffset( 1.2);
	hcomp->GetYaxis()->SetLabelSize( 0.03);
	hcomp->Draw( "colz");
	hcomp->GetXaxis()->SetRangeUser( 0, 40);

	l1->Draw();

	name = Form( "plots/OpenAngleCut");
//	name.Append( ".eps");
	name.Append( ".pdf");
	c1->Print( name);

}
