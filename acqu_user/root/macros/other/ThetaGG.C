gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

TFile full( "histograms/ARH_full.root");
TFile empty( "histograms/ARH_empty.root");

TH3D *h_f_p = (TH3D*)full.Get( "THits_TGG_ThetaCM_P");
TH3D *h_f_r = (TH3D*)full.Get( "THits_TGG_ThetaCM_R");
TH1D *h_s_f = (TH1D*)full.Get( "SumScalers12to363");

TH3D *h_e_p = (TH3D*)empty.Get( "THits_TGG_ThetaCM_P");
TH3D *h_e_r = (TH3D*)empty.Get( "THits_TGG_ThetaCM_R");
TH1D *h_s_e = (TH1D*)empty.Get( "SumScalers12to363");

void ThetaGG( Int_t chan)
{
	Int_t i, j, bin;
	Int_t egamma;
	Double_t pa_ratio[284], pa, yield;
	Double_t emin, emax, y1, ym, yy;
	Double_t junk, eg, energy[284];
	Double_t Tgg, Tgg_c12;
	Double_t r_fe, fscalers, escalers;
	Int_t xmin, xmax;
	Char_t filename[64], name[64], hist[64], string[64];

	// Scalers
	fscalers = h_s_f->GetBinContent( chan);
	escalers = h_s_e->GetBinContent( chan);
	r_fe = fscalers/escalers;
	cout << " full-empty ratio = " << r_fe << endl;

	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 20, 20, 900, 500);
   c1->SetFillColor( 38);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(2,1);

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
	h_f_p->GetXaxis()->SetRangeUser(0,180);	 // The entire ThetaCM range
	h_f_p->GetZaxis()->SetRangeUser(chan,chan);// The appropriate tagger channel
	TH1D *p_f_p = h_f_p->Project3D("y");
	TH1D *c_f_p = (TH1D*)p_f_p->Clone( "fpnew");

	//  Random
	h_f_r->GetXaxis()->SetRangeUser(0,180);	 // The entire ThetaCM range
	h_f_r->GetZaxis()->SetRangeUser(chan,chan);// The appropriate tagger channel
	TH1D *p_f_r = h_f_r->Project3D("y");

	cout << "Subtraction Ratio = " << pa << endl;

	c_f_p->Sumw2();
	c_f_p->Add( p_f_r, -pa);

	c1->cd(1);
	c_f_p->SetMarkerColor(4);
	c_f_p->SetLineColor(4);
	c_f_p->SetLineWidth(4);
	sprintf( string, "#theta_{#gamma#gamma} for E_{#gamma} = %5.1f MeV", eg);
	c_f_p->SetTitle(string);
	c_f_p->Draw();
	c_f_p->GetXaxis()->SetRangeUser(Tgg_c12-10,180);
	TH1D *sub = (TH1D*)c_f_p->Clone( "sub");
	c_f_p->GetXaxis()->SetTitleOffset( 1.1);
	c_f_p->GetXaxis()->SetTitle("#theta_{#gamma#gamma} (deg)");
	c_f_p->GetXaxis()->CenterTitle();

	yy = c_f_p->GetMaximum();
	ym = yy*1.2;
	c_f_p->SetMaximum(ym);

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

	// Empty
	//  Prompt
	h_e_p->GetXaxis()->SetRangeUser(0,180);	 // The entire ThetaCM range
	h_e_p->GetZaxis()->SetRangeUser(chan,chan);// The appropriate tagger channel
	TH1D *p_e_p = h_e_p->Project3D("y");
	TH1D *c_e_p = (TH1D*)p_e_p->Clone( "epnew");

	//  Random
	h_e_r->GetXaxis()->SetRangeUser(0,180);	 // The entire ThetaCM range
	h_e_r->GetZaxis()->SetRangeUser(chan,chan);// The appropriate tagger channel
	TH1D *p_e_r = h_e_r->Project3D("y");

	c_e_p->Sumw2();
	c_e_p->Add( p_e_r, -pa);
	c_e_p->SetLineColor(2);
	c_e_p->SetLineWidth(4);
	c_e_p->Scale( r_fe);
	c_e_p->Draw("same");

	gStyle->SetOptStat(0);

	c1->cd(2);
	sub->Add( c_e_p, -1.0);
	sub->SetLineColor(1);
	sub->SetLineWidth(4);
	sub->SetTitle();
	sub->Draw();
	sub->GetXaxis()->SetRangeUser(Tgg_c12-10,180);
	sub->GetXaxis()->SetTitleOffset( 1.1);
	sub->GetXaxis()->SetTitle("#theta_{#gamma#gamma} (deg)");
	sub->GetXaxis()->CenterTitle();

	yy = sub->GetMaximum();
	ym = yy*1.2;
	sub->SetMaximum(ym);

	l1 = new TLine(Tgg_c12-10,0,180,0);
	l1->SetLineStyle(1);
	l1->SetLineWidth(1);
	l1->Draw();

	l2 = new TLine(Tgg,0,Tgg,ym);
	l2->SetLineStyle(2);
	l2->SetLineWidth(4);
	l2->SetLineColor(4);
	l2->Draw();

	l3 = new TLine(Tgg_c12,0,Tgg_c12,ym);
	l3->SetLineStyle(2);
	l3->SetLineWidth(4);
	l3->SetLineColor(4);
	l3->Draw();


	sprintf( name, "plots/eps/TGG_%d.eps", egamma);
	c1->Print( name);
}
