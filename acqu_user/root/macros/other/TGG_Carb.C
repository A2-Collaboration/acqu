gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

TFile full( "histograms/ARH_full.root");
//TFile empty( "histograms/ARH_empty.root");

TH3D *h_p = (TH3D*)full.Get( "THits_TGG_ThetaCM_C_P");
TH3D *h_r = (TH3D*)full.Get( "THits_TGG_ThetaCM_C_R");

void TGG_Carb( Int_t chan)
{
	Int_t i, j, bin;
	Int_t egamma;
	Double_t pa_ratio[284], pa, yield;
	Double_t emin, emax, y1, ym, yy;
	Double_t junk, eg, energy[284];
	Double_t Tgg, Tgg_c12;
	Int_t xmin, xmax;
	Char_t filename[64], name[64], hist[64], string[64];

	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 20, 20, 900, 500);
   c1->SetFillColor( 38);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	strcpy( filename, "xs/chan_subt/chan_subt.out");
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> pa ;
		pa_ratio[i] = pa;
	}
	inFile.close();
	pa = pa_ratio[chan];

	strcpy( filename, "xs/tageng450.out");
	ifstream inFile( filename);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> junk >> junk >> junk >> eg ;
		energy[i] = eg;
	}
	inFile.close();

	eg = energy[chan];
	egamma = (int)(eg + 0.5);

	if ( eg >= E_thr( kMP_MEV, kMPI0_MEV))
		Tgg = Tgg_Min( qT_max(eg, kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg = 180.;

	if ( eg >= E_thr( kM_C12_MEV, kMPI0_MEV))
		Tgg_c12 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
	else Tgg_c12 = 180.;

	// Full
	//  Prompt
	h_p->GetXaxis()->SetRangeUser(0,180);	 // The entire ThetaCM range
	h_p->GetZaxis()->SetRangeUser(chan,chan);// The appropriate tagger channel
	TH1D *p_p = h_p->Project3D("y");

	//  Random
	h_r->GetXaxis()->SetRangeUser(0,180);	 // The entire ThetaCM range
	h_r->GetZaxis()->SetRangeUser(chan,chan);// The appropriate tagger channel
	TH1D *p_r = h_r->Project3D("y");

	cout << "Subtraction Ratio = " << pa << endl;

	p_p->Sumw2();
	p_p->Add( p_r, -pa);

	p_p->SetMarkerColor(4);
	p_p->SetLineColor(4);
	p_p->SetLineWidth(4);
	sprintf( string, "#theta_{#gamma#gamma} for E_{#gamma} = %5.1f MeV", eg);
	p_p->SetTitle(string);
	p_p->Draw();
	p_p->GetXaxis()->SetRangeUser(Tgg_c12-10,180);

	yy = p_p->GetMaximum();
	ym = yy*1.2;
	p_p->SetMaximum(ym);

	l5 = new TLine(Tgg_c12-10,0,180,0);
	l5->SetLineStyle(1);
	l5->SetLineWidth(1);
	l5->Draw();

	l6 = new TLine(Tgg,0,Tgg,ym);
	l6->SetLineStyle(2);
	l6->SetLineWidth(4);
	l6->Draw();

	l7 = new TLine(Tgg_c12,0,Tgg_c12,ym);
	l7->SetLineStyle(2);
	l7->SetLineWidth(4);
	l7->Draw();

	gStyle->SetOptStat(0);

//	sprintf( name, "plots/eps/TGG_C_%d.eps", egamma);
//	c1->Print( name);
}
