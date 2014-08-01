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

TCanvas *c1 = new TCanvas( "c1", "Random Subtraction", 200, 0, 1000, 700);
c1->SetFillColor( 18);

// The standard parameters are:
// 	LH2:
// 		240 MeV:		4	945	10		kTRUE
// 		280 MeV:		8	940	10		kTRUE
// 	Butanol:
// 		240 MeV:		4	945	5		kFALSE
// 		280 MeV:		8	940	5		kFALSE


void PMPlot( UInt_t tbin)
{
	UInt_t ebin;

	TString name;

	Double_t open_cut = 5;

//	c1->Divide( 4, 3);
	c1->Divide( 3, 2);

//	for ( ebin = 0; ebin <= 10; ebin++)
	for ( ebin = 5; ebin <= 10; ebin++)
	{
//		c1->cd( ebin+1);
		c1->cd( ebin-4);
		CheckCompSum( ebin, tbin, "pi0_p", 5, kFALSE);
	}
	name = Form( "plots/PhotonMmissSim_oa%d_t%d.pdf", (int) open_cut,
			(int) th[tbin]);
//	c1->Print( name);
}

void CheckCompSum( Int_t ebin = 8, Int_t tbin = 0, TString reaction = "pi0_p",
		Double_t open_cut = 5, Bool_t save = kFALSE)
{
	Int_t i, bin_lo, bin_hi, pbin;
	Double_t fact[100][100];
	Double_t left, right, open_cut;
	Double_t x[2], y[2];
	Double_t max, msep;
	Double_t pc, bc, pctot, evp, cont;
	Int_t junk;
	TString name, title;

	gStyle->SetOptStat( 0);

	name = Form( "histograms/MonteCarlo/compton_p_%d.root", eg[ebin]);
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
	if ( reaction == "compton_qfc12") {
		title = "QF #gamma from ^{12}C";
		pbin = 1;
	}
	else if ( reaction == "pi0_qfc12") {
		title = "QF #pi^{0} from ^{12}C";
		pbin = 2;
	}
	else if ( reaction == "compton_qfhe3") {
		title = "QF #gamma from ^{3}He";
		pbin = 3;
	}
	else if ( reaction == "compton_qfhe4") {
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
	name = Form( "evgen/compton_p_%d_in.root", eg[ebin]);
	evgen_p = new TFile( name);
	name = "h4";
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

	name = Form( "#theta^{CM}_{#gamma} = %d #pm 10 deg", th[tbin]);
//	name = Form( "E_{#gamma} = %3d MeV", eg[ebin]);
	pl = new TPaveLabel( 0.15, 0.87, 0.40, 1.00, name, "NDC");
	pl->SetFillStyle(0);
	pl->SetBorderSize(0);
	pl->SetTextSize(0.4);
	pl->SetTextFont(132);
//	if ( ebin == 0) pl->Draw();
	if ( ebin == 5) pl->Draw();

	y[0] = 0.75;
	y[1] = 0.87;
	x[0] = 0.55;
	x[1] = 0.85;
	pt = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pt->SetBorderSize( 0);
	pt->SetFillStyle( 0);
	pt->SetTextSize( 0.04);
	pt->SetTextFont( 132);
	pt->SetTextAlign( 12);
	name = Form( "E_{#gamma} = %3d #pm 1 MeV", eg[ebin]);
	pt->AddText( name);
//	name = Form( "#theta_{open} < %d deg", (int) open_cut);
//	pt->AddText( name);
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
