gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

Int_t eg[] = {150, 158, 165, 180, 195, 210, 225};
Int_t th[] = {10, 30, 50, 70, 90, 110, 130, 150, 170};
Double_t fact[] = {60, 40, 21, 16, 12, 8.1, 4.8};

Double_t msep[7][9];

msep[0][0] = 937.75;
msep[0][1] = 937.75;
msep[0][2] = 937.75;
msep[0][3] = 937;
msep[0][4] = 935.75;
msep[0][5] = 934;
msep[0][6] = 933.5;
msep[0][7] = 930.5;
msep[0][8] = 930;
msep[1][0] = 937.75;
msep[1][1] = 937.75;
msep[1][2] = 937.75;
msep[1][3] = 937;
msep[1][4] = 935.75;
msep[1][5] = 934;
msep[1][6] = 933.5;
msep[1][7] = 930.5;
msep[1][8] = 930;
msep[2][0] = 939.5;
msep[2][1] = 938;
msep[2][2] = 938;
msep[2][3] = 937.5;
msep[2][4] = 934.75;
msep[2][5] = 931.75;
msep[2][6] = 929;
msep[2][7] = 926;
msep[2][8] = 924;
msep[3][0] = 939.5;
msep[3][1] = 938;
msep[3][2] = 938;
msep[3][3] = 938.25;
msep[3][4] = 933.0;
msep[3][5] = 928;
msep[3][6] = 924;
msep[3][7] = 918;
msep[3][8] = 915;

msep[4][0] = 939;
msep[4][1] = 939;
msep[4][2] = 939;
msep[4][3] = 938;
msep[4][4] = 932.5;
msep[4][5] = 925;
msep[4][6] = 918;
msep[4][7] = 910;
msep[4][8] = 910;

msep[5][0] = 940;
msep[5][1] = 940;
msep[5][2] = 940;
msep[5][3] = 938;
msep[5][4] = 932;
msep[5][5] = 922;
msep[5][6] = 912;
msep[5][7] = 910;
msep[5][8] = 910;

msep[6][0] = 940;
msep[6][1] = 940;
msep[6][2] = 940;
msep[6][3] = 935;
msep[6][4] = 927;
msep[6][5] = 916;
msep[6][6] = 910;
msep[6][7] = 910;
msep[6][8] = 910;

void TGGCut( Int_t bin = 0, Int_t tbin = 0)
{
	Int_t i;
	Int_t bin_lo, bin_hi;
	Int_t th_lo, th_hi;
	Double_t pc, cc, max;
	Double_t left, right;
	Double_t qm, pm, q_pi, T_pi, Tgg, Eg, qth, tgg_cut, tgg_inc;
	TString name, proj;

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.025,"xyz");
	gStyle->SetTitleOffset(1.35,"xyz");

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 100, 0, 800, 800);
	c1->Divide(3,3);

	th_lo = th[tbin]-10;
	th_hi = th[tbin]+10;

	Eg = (double) eg[bin];
	qth = (double) th[tbin];
	qm = kMPI0_MEV;

	// For 12C
//	pm = kM_C12_MEV; 
	// For proton
	pm = kMP_MEV; 
	q_pi = qp_thcm( Eg, pm, qth, qm);
	T_pi = Energy( q_pi, qm) - qm;
	Tgg = Tgg_Min( T_pi, qm)/kD2R-7;

	cout << qth;
	cout << " ";
	cout << q_pi;
	cout << " ";
	cout << T_pi;
	cout << " ";
	cout << Tgg;
	cout << endl;

	tgg_inc = (180 - Tgg)/8;

	name = Form( "histograms/MCH_pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MCH_pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);

	name = "MMissPi0_TGG_ThetaCM_P";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

	bin_lo = h3p->GetXaxis()->FindBin(th_lo);
	bin_hi = h3p->GetXaxis()->FindBin(th_hi)-1;
	h3p->GetXaxis()->SetRange( bin_lo, bin_hi);
	h3c->GetXaxis()->SetRange( bin_lo, bin_hi);

//	left = msep[bin][tbin] - 2;
//	right = msep[bin][tbin] + 5;
	left = msep[bin][tbin] - 10;
//	right = msep[bin][tbin] + 5;
//	left = 900;
	right = 960;

	h3p->GetZaxis()->SetRangeUser( left, right);
	h3c->GetZaxis()->SetRangeUser( left, right);

	name = Form( "efficiency/efficiency_%d_%d.out", eg[bin], th[tbin]);
	ofstream outFile( name);
	if ( !outFile.is_open())
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit( -1);
	}

// Event Generator Information
	name = Form( "simfiles/pi0/pi0_p_%d_out.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "simfiles/pi0/pi0_c_%d_out.root", eg[bin]);
	carb = new TFile( name);

	name = "h5";
	ev_p = (TH1D*) prot->Get( name);
	ev_c = (TH1D*) carb->Get( name);

	ev_p->Scale( 1/fact[bin]);

	bin_lo = ev_p->GetXaxis()->FindBin(th_lo);
	bin_hi = ev_p->GetXaxis()->FindBin(th_hi)-1;

	Double_t evp = ev_p->Integral( bin_lo, bin_hi);
	Double_t evc = ev_c->Integral( bin_lo, bin_hi);

	outFile << "total counts in evgen ";
	outFile << evp;
	outFile << " ";
	outFile << evc;
	outFile << endl;

	i = 0;
	for ( tgg_cut = Tgg; tgg_cut < 180; tgg_cut += tgg_inc)
	{
		c1->cd(++i)->SetLogy();
//		c1->cd(++i);

		h3p->GetYaxis()->SetRangeUser( tgg_cut, 180);
		h3c->GetYaxis()->SetRangeUser( tgg_cut, 180);

		// Plot Missing Mass
		proj = "z";

		TH1D *h1p = h3p->Project3D( proj);
		proj.Append("2");
		TH1D *h1c = h3c->Project3D( proj);

		if ( i == 1)
		{
			name = Form( "E_{#gamma} = %3d MeV, #theta_{#pi}^{CM} = %d-%d deg,",
					eg[bin], th_lo, th_hi);
			tl = new TLegend( 0.75, 0.6, 0.90, 0.8);
			tl->SetFillStyle(0);
			tl->SetBorderSize(0);
			tl->SetTextSize(0.04);

			tl->AddEntry( h1c, "^{12}C");
			tl->AddEntry( h1p, "p");
//			tl->AddEntry( h3, "sum", "P");

			tl->Draw();
		}
		else name = "";

		// Proton
		h1p->SetFillStyle(3001);
		h1p->SetLineWidth(2);
		h1p->SetFillColor(2);
		h1p->Scale( 1/fact[bin]);

		// Carbon
		h1c->SetFillStyle(3001);
		h1c->SetFillColor(4);
		h1c->SetLineWidth(2);

		pc = h1p->Integral();
		cc = h1c->Integral();

		outFile << tgg_cut;
		outFile << " ";
		outFile << pc;
		outFile << " ";
		outFile << cc;

		if ( h1c->GetMaximum() > h1p->GetMaximum())
		{
			h1c->SetTitle( name);
			max = h1c->GetMaximum();
			max *= 1.10;
			h1c->SetMaximum( max);
			h1c->DrawClone();
			h1p->DrawClone("same");
		}
		else
		{
			h1p->SetTitle( name);
			max = h1p->GetMaximum();
			max *= 1.10;
			h1p->SetMaximum( max);
			h1p->DrawClone();
			h1c->DrawClone("same");
		}
//		cout << "max = " << max << endl;

		l1 = new TLine( msep[bin][tbin], 0, msep[bin][tbin], max);
		l1->SetLineColor(1);
		l1->SetLineStyle(2);
		l1->SetLineWidth(4);
		l1->Draw();

		pl = new TPaveText( 0.5, 0.7, 0.8, 0.9, "NDC");
		pl->SetBorderSize(0);
		pl->SetFillStyle(0);
		pl->SetTextAlign(12);
		pl->SetTextSize(0.04);
		name = Form( "#theta_{#gamma#gamma} = %d - 180 deg\n", (int) tgg_cut);
		pl->AddText( name);
		pl->Draw();

		Int_t midc, midp, cbins, pbins;
		Double_t yp, yc;

		midc = h1c->GetXaxis()->FindBin( msep[bin][tbin]);
		midp = h1p->GetXaxis()->FindBin( msep[bin][tbin]);
		
		cbins = h1c->GetXaxis()->GetNbins();
		pbins = h1c->GetXaxis()->GetNbins();

		yc = h1c->Integral(midc, cbins);
		yp = h1p->Integral(midp, pbins);

		outFile << " ";
		outFile << yp;
		outFile << " ";
		outFile << yc;
		outFile << " ";
		if ( yp != 0) outFile << yc/yp;
		else outFile << " 0.00";
		outFile << " ";
		outFile << yp/evp;
		outFile << endl;

	}
	outFile.close();

// Output pdf
	name = Form( "efficiency/eff/MmissCut_%d_%d.pdf", eg[bin], th[tbin]);
	c1->Print( name);
}

void TGGAll( Int_t bin)
{
	Int_t tbin;

	for ( tbin = 0; tbin <= 8; tbin++)
	{
		TGGCut( bin, tbin);
	}
}

void TGG2DCut( Int_t bin = 0, TString tgt = "p")
{
	Int_t i, tbin;
	Int_t th_lo, th_hi;
	Double_t Tgg_c12;
	TString name, proj;

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.025,"xyz");
	gStyle->SetTitleOffset(1.35,"xyz");

	Tgg_c12 = Tgg_Min( qT_max(eg[bin], kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 100, 0, 800, 800);
	c1->Divide(3,3);

	name = Form( "histograms/MCH_pi0_%s_%d.root", tgt.Data(), eg[bin]);
	sim = new TFile( name);
	name = "MMissPi0_TGG_ThetaCM_NC_P";
	h3 = (TH3D*) sim->Get( name);

	h3->GetYaxis()->SetRangeUser( Tgg_c12-10, 180);
	h3->GetZaxis()->SetRangeUser( 900, 960);

	name = Form( "evgen/pi0_%s_%d_in.root", tgt.Data(), eg[bin]);
	ev = new TFile( name);
	name = "h5";
	hev = (TH1D*) ev->Get( name);

	Int_t bin_lo, bin_hi;

	i = 0;
	Int_t totc = 0;
	Int_t tot_evc = 0;
	cout << h3->Integral();
	cout << " ";
	cout << hev->Integral();
	cout << endl;
	for ( tbin = 0; tbin < 9; tbin++)
	{
		c1->cd(++i);

		th_lo = th[tbin]-10;
		th_hi = th[tbin]+10;

		bin_lo = h3->GetXaxis()->FindBin( th_lo);
		bin_hi = h3->GetXaxis()->FindBin( th_hi)-1;
		h3->GetXaxis()->SetRange( bin_lo, bin_hi);

		proj = "yz";
		TH2D *h2 = h3->Project3D( proj);

		Int_t counts = h2->Integral();
		totc += counts;
		cout << counts;

		bin_lo = hev->FindBin(th_lo);
		bin_hi = hev->FindBin(th_hi)-1;
		Int_t evcounts = hev->Integral( bin_lo, bin_hi);
		tot_evc += evcounts;
		cout << " ";
		cout << evcounts;

		Double_t eff = (double) counts/evcounts;
		cout << " ";
		cout << eff;
		cout << endl;

		if ( i == 1)
		{
			name = Form( "E_{#gamma} = %3d MeV, Target = %s",
					eg[bin], tgt.Data(), th_lo, th_hi);
		}
		else name = "";
		h2->SetTitle( name);

		h2->GetXaxis()->SetTitle("M_{miss} (MeV)");
		h2->GetYaxis()->SetTitle("#theta_{#gamma#gamma} (deg)");
		h2->DrawClone("col");

		pl = new TPaveText( 0.2, 0.15, 0.5, 0.35, "NDC");
		pl->SetBorderSize(0);
		pl->SetFillStyle(0);
		pl->SetTextAlign(12);
		pl->SetTextSize(0.04);
		name = Form( "#theta_{#pi}^{CM} = %d - %d deg\n", th_lo, th_hi);
		pl->AddText( name);
		name = Form( "Counts = %d\n", counts);
		pl->AddText( name);
		pl->Draw();

	}
	cout << " tot counts = " << totc;
	cout << " " << tot_evc;
	cout << " " << (double) totc/tot_evc;
	cout << endl;

	name = Form( "plots/pdf/TGG_v_Mmiss_%s_%d.pdf", tgt.Data(), eg[bin]);
	c1->Print( name);
}

void TGG2DCutSingle( Int_t tbin = 2)
{
	Int_t bin;
	Int_t bin_lo, bin_hi;
	Int_t th_lo, th_hi;
	Int_t left, right;
	Double_t pm, qm, Eg, qth, T_pi, q_pi;
	Double_t Tgg, Tgg_c12;
	TString tgt, name, proj;

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.045,"xyz");
	gStyle->SetTitleSize(0.045,"xyz");
	gStyle->SetTitleOffset(1.0,"xyz");

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 100, 0, 800, 1000);
	c1->Divide(2,7);

	th_lo = th[tbin]-10;
	th_hi = th[tbin]+10;

	left = 900;
	right = 960;
//	left = 920;
//	right = 945;

	qth = (double) th[tbin];
	qm = kMPI0_MEV;
	pm = kMP_MEV; 

	// Proton
	tgt = "p";
	for ( bin = 0; bin < 7; bin++)
	{
		c1->cd(2*bin+1);

		Eg = (double) eg[bin];
		q_pi = qp_thcm( Eg, pm, qth, qm);
		T_pi = Energy( q_pi, qm) - qm;
		Tgg = Tgg_Min( T_pi, qm)/kD2R-7;

		name = Form( "histograms/MCH_%s_%d_pi0.root",
				tgt.Data(), eg[bin]);
		sim = new TFile( name);
		name = "MMissPi0_TGG_ThetaCM_P";
		h3p = (TH3D*) sim->Get( name);

		bin_lo = h3p->GetXaxis()->FindBin( th_lo);
		bin_hi = h3p->GetXaxis()->FindBin( th_hi)-1;
		h3p->GetXaxis()->SetRange( bin_lo, bin_hi);

		Tgg_c12 = Tgg_Min( qT_max(eg[bin], kM_C12_MEV, kMPI0_MEV),
				kMPI0_MEV)/kD2R;
		h3p->GetYaxis()->SetRangeUser( Tgg_c12-10, 180);
		h3p->GetZaxis()->SetRangeUser( left, right);

		proj = "yz";
		TH2D *h2p = h3p->Project3D( proj);

		if ( bin == 0)
			name = Form( "  #theta_{#pi}^{CM} = %d - %d deg \n", th_lo, th_hi);
		else name = "";
		h2p->SetTitle( name);

		h2p->GetXaxis()->SetTitle("M_{miss} (MeV)");
		h2p->GetXaxis()->CenterTitle();
		h2p->GetYaxis()->SetTitle("#theta_{#gamma#gamma} (deg)");
		h2p->GetYaxis()->SetTitleOffset( 0.9);
		h2p->GetYaxis()->CenterTitle();
		h2p->DrawClone("col");

		pt = new TPaveText( 0.15, 0.7, 0.35, 0.9, "NDC");
		pt->SetBorderSize(0);
		pt->SetFillStyle(0);
		pt->SetTextAlign(12);
		pt->SetTextSize(0.05);
		name = Form( "E_{#gamma} = %d MeV\n", eg[bin]);
		pt->AddText( name);
		pt->Draw();

		if ( bin == 0) {
			pl = new TPaveLabel( 0.65, 0.9, 0.85, 0.98, "Proton", "NDC");
			pl->SetBorderSize(0);
			pl->SetFillStyle(0);
			pl->SetTextSize(0.7);
			pl->Draw();
		}

		l1 = new TLine( msep[bin][tbin], Tgg_c12 - 10, msep[bin][tbin], 180);
		l2 = new TLine( left, Tgg, right, Tgg);
		l1->Draw();
		l2->Draw();
	}

	// Carbon
	tgt = "c";
	for ( bin = 0; bin < 7; bin++)
	{
		c1->cd(2*bin+2);

		Eg = (double) eg[bin];
		q_pi = qp_thcm( Eg, pm, qth, qm);
		T_pi = Energy( q_pi, qm) - qm;
		Tgg = Tgg_Min( T_pi, qm)/kD2R-7;

		name = Form( "histograms/MCH_%s_%d_pi0.root",
				tgt.Data(), eg[bin]);
		sim = new TFile( name);
		name = "MMissPi0_TGG_ThetaCM_P";
		h3c = (TH3D*) sim->Get( name);

		bin_lo = h3c->GetXaxis()->FindBin( th_lo);
		bin_hi = h3c->GetXaxis()->FindBin( th_hi)-1;
		h3c->GetXaxis()->SetRange( bin_lo, bin_hi);

		Tgg_c12 = Tgg_Min( qT_max(eg[bin], kM_C12_MEV, kMPI0_MEV),
				kMPI0_MEV)/kD2R;
		h3c->GetYaxis()->SetRangeUser( Tgg_c12-10, 180);
		h3c->GetZaxis()->SetRangeUser( left, right);

		proj = "yz";
		TH2D *h2c = h3c->Project3D( proj);

		h2c->SetTitle();

		h2c->GetXaxis()->SetTitle("M_{miss} (MeV)");
		h2c->GetXaxis()->CenterTitle();
		h2c->GetYaxis()->SetTitle("#theta_{#gamma#gamma} (deg)");
		h2c->GetYaxis()->SetTitleOffset( 0.9);
		h2c->GetYaxis()->CenterTitle();
		h2c->DrawClone("col");

		pt = new TPaveText( 0.15, 0.7, 0.35, 0.9, "NDC");
		pt->SetBorderSize(0);
		pt->SetFillStyle(0);
		pt->SetTextAlign(12);
		pt->SetTextSize(0.05);
		name = Form( "E_{#gamma} = %d MeV\n", eg[bin]);
		pt->AddText( name);
		pt->Draw();

		if ( bin == 0) {
			pl = new TPaveLabel( 0.65, 0.9, 0.85, 0.98, "Carbon", "NDC");
			pl->SetBorderSize(0);
			pl->SetFillStyle(0);
			pl->SetTextSize(0.7);
			pl->Draw();
		}

		l1 = new TLine( msep[bin][tbin], Tgg_c12 - 10, msep[bin][tbin], 180);
		l2 = new TLine( left, Tgg, right, Tgg);
		l1->Draw();
		l2->Draw();
	}


	name = Form( "efficiency/eff/TGG_v_Mmiss_%d.pdf", th[tbin]);
	c1->Print( name);
}

void AronOut( Int_t bin = 0, Int_t tbin = 4)
{
	Int_t i;
	Int_t bin_lo, bin_hi;
	Int_t th_lo, th_hi;
	Double_t pc, cc;
	Double_t tgg_cut;
	TString name, proj;

	th_lo = th[tbin]-10;
	th_hi = th[tbin]+10;

	name = Form( "histograms/MCH_pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MCH_pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);

	name = "MMissPi0_TGG_ThetaCM_NC_P";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

	bin_lo = h3p->GetXaxis()->FindBin(th_lo);
	bin_hi = h3p->GetXaxis()->FindBin(th_hi)-1;
	h3p->GetXaxis()->SetRange( bin_lo, bin_hi);
	h3c->GetXaxis()->SetRange( bin_lo, bin_hi);

	h3p->GetZaxis()->SetRangeUser( 900, 960);
	h3c->GetZaxis()->SetRangeUser( 900, 960);

	tgg_cut = 137;
	h3p->GetYaxis()->SetRangeUser( tgg_cut, 180);
	h3c->GetYaxis()->SetRangeUser( tgg_cut, 180);

	// Plot Missing Mass
	proj = "z";

	TH1D *h1p = h3p->Project3D( proj);
	h1p->Scale( 1/fact[bin]);

	proj.Append("2");
	TH1D *h1c = h3c->Project3D( proj);

	pc = h1p->Integral();
	cc = h1c->Integral();

	ofstream outFile( "aron.out");
	outFile << "Eg = ";
	outFile << eg[bin];
	outFile << " MeV";
	outFile << endl;
	outFile << "cm pion theta = ";
	outFile << th_lo;
	outFile << " - ";
	outFile << th_hi;
	outFile << " deg";
	outFile << endl;
	outFile << "tgg min = ";
	outFile << tgg_cut;
	outFile << " deg";
	outFile << endl;
	outFile << "proton counts = ";
	outFile << pc;
	outFile << endl;
	outFile << "12C counts = ";
	outFile << cc;
	outFile << endl;

	Int_t nbins = h1c->GetXaxis()->GetNbins();
	for ( i = 0; i < nbins; i++)
	{
		if ( !( ( h1p->GetBinContent(i) == 0) && ( h1c->GetBinContent(i) == 0))) {
			outFile << h1p->GetBinCenter(i);
			outFile << " ";
			outFile.precision(2);
			outFile.width(5);
			outFile << h1p->GetBinContent(i);
			outFile << " ";
			outFile.precision(4);
			outFile.width(4);
			outFile << h1c->GetBinContent(i);
			outFile << endl;
		}
	}
	outFile.close();
}

void TGGCutSingle( Int_t bin = 0, Int_t tbin = 0)
{
	Int_t i;
	Int_t bin_lo, bin_hi;
	Int_t th_lo, th_hi;
	Double_t pc, cc, max;
	Double_t qm, pm, q_pi, T_pi, Tgg, Eg, qth, tgg_cut, tgg_inc;
	TString name, proj;

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.025,"xyz");
	gStyle->SetTitleOffset(1.35,"xyz");

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 800, 0, 700, 500);

	th_lo = th[tbin]-10;
	th_hi = th[tbin]+10;

	Eg = (double) eg[bin];
	qth = (double) th[tbin];
	qm = kMPI0_MEV;

	// For 12C
//	pm = kM_C12_MEV; 
	// For proton
	pm = kMP_MEV; 
	q_pi = qp_thcm( Eg, pm, qth, qm);
	T_pi = Energy( q_pi, qm) - qm;
	Tgg = Tgg_Min( T_pi, qm)/kD2R;

	tgg_inc = (180 - Tgg)/8;

	name = Form( "histograms/MCH_pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MCH_pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);

	name = "MMissPi0_TGG_ThetaCM_P";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

	bin_lo = h3p->GetXaxis()->FindBin(th_lo);
	bin_hi = h3p->GetXaxis()->FindBin(th_hi)-1;
	h3p->GetXaxis()->SetRange( bin_lo, bin_hi);
	h3c->GetXaxis()->SetRange( bin_lo, bin_hi);

	Int_t left, right;

	left = msep[bin][tbin] - 5;
	right = msep[bin][tbin] + 5;

	h3p->GetZaxis()->SetRangeUser( left, right);
	h3c->GetZaxis()->SetRangeUser( left, right);

	tgg_cut = Tgg + 2*tgg_inc;

	c1->cd(1)->SetLogy();

	h3p->GetYaxis()->SetRangeUser( tgg_cut, 180);
	h3c->GetYaxis()->SetRangeUser( tgg_cut, 180);

	// Plot Missing Mass
	proj = "z";

	TH1D *h1p = h3p->Project3D( proj);
	proj.Append("2");
	TH1D *h1c = h3c->Project3D( proj);

	name = Form( "E_{#gamma} = %3d MeV, #theta_{#pi}^{CM} = %d-%d deg,",
			eg[bin], th_lo, th_hi);
	tl = new TLegend( 0.75, 0.6, 0.90, 0.8);
	tl->SetFillStyle(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.04);

	tl->AddEntry( h1c, "^{12}C");
	tl->AddEntry( h1p, "p");

	tl->Draw();

	// Proton
	h1p->SetFillStyle(3001);
	h1p->SetLineWidth(2);
	h1p->SetFillColor(2);
	h1p->Scale( 1/fact[bin]);

	// Carbon
	h1c->SetFillStyle(3001);
	h1c->SetFillColor(4);
	h1c->SetLineWidth(2);

	if ( h1c->GetMaximum() > h1p->GetMaximum())
	{
		h1c->SetTitle( name);
		max = h1c->GetMaximum();
		max *= 1.10;
		h1c->SetMaximum( max);
		h1c->DrawClone();
		h1p->DrawClone("same");
	}
	else
	{
		h1p->SetTitle( name);
		max = h1p->GetMaximum();
		max *= 1.10;
		h1p->SetMaximum( max);
		h1p->DrawClone();
		h1c->DrawClone("same");
	}

	l1 = new TLine( msep[bin][tbin], 0, msep[bin][tbin], max);
	l1->SetLineColor(1);
	l1->SetLineStyle(2);
	l1->SetLineWidth(4);
	l1->Draw();

	pl = new TPaveText( 0.5, 0.7, 0.8, 0.9, "NDC");
	pl->SetBorderSize(0);
	pl->SetFillStyle(0);
	pl->SetTextAlign(12);
	pl->SetTextSize(0.04);
	name = Form( "#theta_{#gamma#gamma} = %d - 180 deg\n", (int) tgg_cut);
	pl->AddText( name);
	pl->Draw();

// Output pdf
//	name = Form( "plots/pdf/MmissCut_%d_%d.pdf", eg[bin], th[tbin]);
//	c1->Print( name);
}

void TGGCutProp( Int_t bin = 0, Int_t tbin = 4)
{
	Int_t i;
	Int_t bin_lo, bin_hi;
	Int_t th_lo, th_hi;
	Double_t pc, cc;
	Double_t qm, pm, q_pi, T_pi, Tgg, Eg, qth, tgg_cut, tgg_inc;
	Double_t max;
	Double_t left, right;
	TString name, proj;

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.035,"xyz");
	gStyle->SetTitleOffset(1.30,"xyz");

	TCanvas *c1 = new TCanvas ( "c1", "CH2 Comparison", 400, 0, 700, 500);

//	c1->cd(1)->SetLogy();

	th_lo = th[tbin]-10;
	th_hi = th[tbin]+10;

	Eg = (double) eg[bin];
	qth = (double) th[tbin];
	qm = kMPI0_MEV;

	// For 12C
//	pm = kM_C12_MEV; 
	// For proton
	pm = kMP_MEV; 
	q_pi = qp_thcm( Eg, pm, qth, qm);
	T_pi = Energy( q_pi, qm) - qm;
	Tgg = Tgg_Min( T_pi, qm)/kD2R - 7;

	tgg_inc = (180 - Tgg)/8;

	name = Form( "histograms/MCH_pi0_p_%d.root", eg[bin]);
	prot = new TFile( name);
	name = Form( "histograms/MCH_pi0_c_%d.root", eg[bin]);
	carb = new TFile( name);

	name = "MMissPi0_TGG_ThetaCM_P";
	h3p = (TH3D*) prot->Get( name);
	h3c = (TH3D*) carb->Get( name);

	bin_lo = h3p->GetXaxis()->FindBin(th_lo);
	bin_hi = h3p->GetXaxis()->FindBin(th_hi)-1;
	h3p->GetXaxis()->SetRange( bin_lo, bin_hi);
	h3c->GetXaxis()->SetRange( bin_lo, bin_hi);

	left = msep[bin][tbin] - 2.5;
	right = msep[bin][tbin] + 15;
//	left = 900;
//	right = 960;

	h3p->GetZaxis()->SetRangeUser( left, right);
	h3c->GetZaxis()->SetRangeUser( left, right);

	tgg_cut = Tgg;
//	c1->cd(1)->SetLogy();
//	c1->cd(1);

	h3p->GetYaxis()->SetRangeUser( tgg_cut, 180);
	h3c->GetYaxis()->SetRangeUser( tgg_cut, 180);

	// Plot Missing Mass
	proj = "z";

	TH1D *h1p = h3p->Project3D( proj);
	proj.Append("2");
	TH1D *h1c = h3c->Project3D( proj);

	// Proton
	h1p->SetFillStyle(3001);
	h1p->SetLineWidth(2);
//	h1p->SetLineColor(2);
	h1p->SetFillColor(17);
	h1p->Scale( 1/fact[bin]);

	// Carbon
	h1c->SetFillStyle(3008);
	h1c->SetFillColor(13);
//	h1c->SetLineColor(4);
	h1c->SetLineWidth(2);

	if ( h1c->GetMaximum() > h1p->GetMaximum())
	{
		h1c->SetTitle();
		h1c->GetXaxis()->SetTitle("M_{miss} (MeV)");
		h1c->GetXaxis()->CenterTitle();
		max = h1c->GetMaximum();
		max *= 1.10;
		h1c->SetMaximum( max);
		h1c->DrawClone();
		h1p->DrawClone("same");
	}
	else
	{
		h1p->SetTitle();
		h1p->GetXaxis()->SetTitle("M_{miss} (MeV)");
		h1p->GetXaxis()->CenterTitle();
		max = h1p->GetMaximum();
		max *= 1.10;
		h1p->SetMaximum( max);
		h1p->DrawClone();
		h1c->DrawClone("same");
	}
//	cout << "max = " << max << endl;

	l1 = new TLine( msep[bin][tbin], 0, msep[bin][tbin], max);
	l1->SetLineColor(1);
	l1->SetLineStyle(2);
	l1->SetLineWidth(2);
	l1->Draw();

	pl = new TPaveText( 0.65, 0.62, 0.9, 0.87, "NDC");
	pl->SetBorderSize(0);
	pl->SetFillStyle(0);
	pl->SetTextAlign(12);
	pl->SetTextSize(0.04);
	pl->SetTextFont(132);
	name = Form( "E_{#gamma} = %3d #pm 2.5 MeV", eg[bin]);
	pl->AddText( name);
	name = Form( "#theta_{#pi}^{CM} = %d #pm 10 deg", th[tbin]);
	pl->AddText( name);
	name = Form( "#theta_{#gamma#gamma} > %d deg\n", (int) tgg_cut);
	pl->AddText( name);
	name = Form( "M_{miss} > %5.1f MeV\n", msep[bin][tbin]);
	pl->AddText( name);
	pl->Draw();

	tl = new TLegend( 0.65, 0.38, 0.9, 0.58);
	tl->SetFillStyle(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.04);
	tl->AddEntry( h1c, "^{12}C");
	tl->AddEntry( h1p, "p");
	tl->Draw();

// Output pdf
	name = Form( "efficiency/eff/MmissCutProp_%d_%d.pdf", eg[bin], th[tbin]);
	c1->Print( name);
}
