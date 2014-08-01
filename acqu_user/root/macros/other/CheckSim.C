gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

void CheckSim( Int_t eg, Char_t hist[32] = "TGG", Char_t cond[32] = "masscorr")
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, x3, a, b, fact;
	Char_t name[64];

	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, eg);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	carb = new TFile( name);

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	if ( strcmp( hist, "TGG") == 0) {

		h3p = (TH3D*) prot->Get( "EmPi0_TGG_ThetaCM");
		TH1D *h1 = h3p->Project3D("y");
		h3c = (TH3D*) carb->Get( "EmPi0_TGG_ThetaCM");
		TH1D *h2 = h3c->Project3D("y2");

		sprintf( name, "#Theta^{min}_{#gamma#gamma} for E_{#gamma} = %5.1f MeV",
				eg);
		h1->SetTitle( name);
		x1 = Tgg_Min( qT_max(eg, kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
		x2 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
		col = 4;
		a = x2 - 10;
		b = 180;
	}
	else if ( strcmp( hist, "EmPi0") == 0) {

		h3p = (TH3D*) prot->Get( "EmPi0_TGG_ThetaCM");
		TH1D *h1 = h3p->Project3D("z");
		h3c = (TH3D*) carb->Get( "EmPi0_TGG_ThetaCM");
		TH1D *h2 = h3c->Project3D("z2");

		sprintf( name, "M_{miss} for E_{#gamma} = %5.1f MeV", eg);
		h1->SetTitle( name);
		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);
		col = 2;
		a = Linear( 145, 900, 400, 600, eg);
		b = Linear( 145, 950, 400, 1000, eg);
	}

	h1->SetLineWidth(2);
	h1->SetFillColor(2);
	h1->GetXaxis()->SetRangeUser(a,b);
	h1->Draw();

	ch1 = h1->Integral();
	ch2 = h2->Integral();
	fact = (double) ch1/ch2;

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

	if ( strcmp( hist, "EmPi0") == 0) {
		l3 = new TLine( x3, 0, x3, max);
		l3->SetLineColor(col);
		l3->SetLineStyle(2);
		l3->SetLineWidth(4);
		l3->Draw();
	}

	sprintf( name, "plots/eps/Sim_%s_%s_%d.eps", hist, cond, eg);
	c1->Print( name);
}

void CheckSim2D( Int_t eg, Char_t proj[3], Char_t cond[32] = "masscorr")
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, a, b, fact;
	Char_t name[64];

	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, eg);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	carb = new TFile( name);

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 500, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(1,2);

	h3p = (TH3D*) prot->Get( "EmPi0_TGG_ThetaCM_P");
	TH2D *h2p = h3p->Project3D( proj);
	h3c = (TH3D*) carb->Get( "EmPi0_TGG_ThetaCM_P");
	strcat( proj, proj, "2");
	TH2D *h2c = h3c->Project3D( proj);

	c1->cd(1);
	h2p->Draw("col");

	c1->cd(2);
	h2c->Draw("col");

//	sprintf( name, "plots/eps/Sim_%s_%s_%d.eps", hist, cond, eg);
//	c1->Print( name);
}

void CheckSimPhi( Int_t eg)
{
	Char_t name[64], hist[32];
	Char_t cond[] = "masscorr";

	sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, eg);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	carb = new TFile( name);

	TCanvas *c1 = new TCanvas ( "c1", "Phi Simulation Comparison", 200, 0, 500,
			700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(1,2);

	strcpy( hist, "THits_PhiCM_P");
	h2p = (TH2D*) prot->Get( hist);
	TH1D *h1 = h2p->ProjectionX( "prot");
//	h2c = (TH2D*) carb->Get( hist);
//	TH1D *h2 = h2c->ProjectionX( "carb");

	c1->cd(1);
	h1->SetLineWidth(2);
	h1->SetMinimum(0);
	h1->Draw("E1P");

//	c1->cd(2);
//	h2->SetLineWidth(2);
//	h2->SetMinimum(0);
//	h2->Draw("E1P");
}

void CheckSimTgt( Int_t eg, Char_t tgt[2] = "p", Char_t hist[32] = "TGG",
		Char_t cond[32] = "masscorr")
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, x3, a, b, fact;
	Char_t name[64];

	sprintf( name, "histograms/MonteCarlo/MCH_%s_%s_%d.root", tgt, cond, eg);
	dat = new TFile( name);

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	if ( strcmp( hist, "TGG") == 0) {

		h3 = (TH3D*) dat->Get( "EmPi0_TGG_ThetaCM");
		TH1D *h = h3->Project3D("y");

		sprintf( name, "#Theta^{min}_{#gamma#gamma} for E_{#gamma} = %5.1f MeV",
				eg);
		h->SetTitle( name);
		x1 = Tgg_Min( qT_max(eg, kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
		x2 = Tgg_Min( qT_max(eg, kM_C12_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;
		col = 4;
		a = x2 - 50;
		b = 180;
	}
	else if ( strcmp( hist, "EmPi0") == 0) {

		h3 = (TH3D*) dat->Get( "EmPi0_TGG_ThetaCM");
		TH1D *h = h3->Project3D("z");

		sprintf( name, "M_{miss} for E_{#gamma} = %5.1f MeV", eg);
		h->SetTitle( name);
		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);
		col = 2;
		a = Linear( 135, 900, 400, 600, eg);
		b = Linear( 135, 950, 400, 1000, eg);
	}

	h->SetFillStyle(3001);
	h->SetFillColor(4);
	h->SetLineWidth(2);
	h->SetLineStyle(1);
	h->GetXaxis()->SetRangeUser(a,b);
	h->Draw();

	max = h->GetMaximum();

	max *= 1.10;
	h->SetMaximum(max);

	l = new TLine( x1, 0, x1, max);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	l = new TLine( x2, 0, x2, max);
	l->SetLineColor(col);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	if ( strcmp( hist, "EmPi0") == 0) {
		l = new TLine( x3, 0, x3, max);
		l->SetLineColor(col);
		l->SetLineStyle(2);
		l->SetLineWidth(4);
		l->Draw();
	}

//	sprintf( name, "plots/eps/Sim_%s_%s_%d.eps", hist, cond, eg);
//	c1->Print( name);
}

void ContEff( Int_t pbin)
{
	Int_t tbin, plot;
	Double_t eff, tot;
	Char_t name[128];

	plot = 0;
	sprintf( name, "xs/eff/carbcont_%d.out", pbin);
	ofstream outFile( name);
	tot = 0;
	for ( tbin = 1; tbin <= 9; tbin++ ) {
		eff = CarbCont( tbin, pbin, plot);
		cout << tbin*20-10 << "  " << eff << endl;
		outFile << tbin*20-10 << "  " << eff << endl;
		tot += eff;
	}
	cout << tot << endl;
	outFile.close();
}

Double_t CarbCont( Int_t tbin, Int_t pbin, Int_t plot)
{
	Int_t theta_lo, theta_hi;
	Int_t zl, zh, yl, yh;
	Int_t zcp, xcp, ycp;
	Int_t zcc, xcc, ycc;
	Double_t x1, x2, a, b, eg, effc, effp, ac, ap, max;
	Char_t name[64], cond[32] = "masscorr";

	if ( pbin == 2) {
		eg = 147;
		ap = 0.9089;
	}
	else if ( pbin == 3) {
		eg = 153;
		ap = 0.8816;
	}
	else if ( pbin == 4) {
		eg = 157;
		ap = 0.8567;
	}
	else if ( pbin == 5) {
		eg = 163;
		ap = 0.8482;
	}
	else if ( pbin == 6) {
		eg = 168;
		ap = 0.8164;
	}
	else {
		printf( "<E> CarbCont: Invalid Channel Number (2-6)\n");
		exit(-1);
	}
	ac = 1-ap;
	
	x1 = Tgg_Min( qT_max(eg, kMP_MEV, kMPI0_MEV), kMPI0_MEV)/kD2R;

	a = Linear( 145, 933, 400, 880, eg);
	b = Linear( 145, 945, 400, 1000, eg);

	theta_lo = (tbin-1)*20;
	theta_hi = tbin*20;

	Int_t pflag = 0;
	if ( pflag == 1) {
		sprintf( name, "histograms/MonteCarlo/MCH_p_%s_%d.root", cond, eg);
		pdat = new TFile( name);
		h3p = (TH3D*) pdat->Get( "EmPi0_TGG_ThetaCM");

		TH2D *h2p = h3p->Project3D("xz1");

		h3p->GetZaxis()->SetRangeUser(a,b);
		TH1D *hzp = h3p->Project3D( "z1");
		zl = hzp->GetXaxis()->FindBin( a);
		zh = hzp->GetXaxis()->FindBin( b);
		zcp = hzp->Integral(zl,zh);

		// zcp is the total number of proton counts in the energy window, i.e.
		// y^{tot}_p.  It is NOT a function of theta.  It is a function of the
		// tagger energy bin, and the missing mass window (which should not
		// chance).

		h3p->GetXaxis()->SetRangeUser(theta_lo,theta_hi);
//		TH1D *hxp = h3p->Project3D("x1");
//		xcp = hxp->Integral();

		// xcp is the total number of proton counts in the energy window AND
		// theta window.  It is a function of theta, and obviously tagger energy
		// bin and missing mass window.

		TH1D *hyp = h3p->Project3D("y1");
//		yl = hyp->GetXaxis()->FindBin( x1);
//		yh = hyp->GetXaxis()->FindBin( 180);
		ycp = hyp->Integral(x1,180);

		// ycp is the number of proton counts in the energy window AND theta
		// window that pass the tgg cut, i.e. y_p(theta).  It is also a function
		// of theta.

		// The efficiency factor f_p(theta) is then ycp/zcp

		effp = ap*ycp/zcp;

		cout << zcp;
//		cout << "  " << xcp;
		cout << "  " << ycp;
		cout << "  " << effp;
		cout << endl;
	}

	sprintf( name, "histograms/MonteCarlo/MCH_c_%s_%d.root", cond, eg);
	cdat = new TFile( name);
	h3c = (TH3D*)cdat->Get( "EmPi0_TGG_ThetaCM");

	TH2D *h2c = h3c->Project3D("xz2");
	TH3D *h3cc = (TH3D*)h3c->Clone( "3d");
	TH3D *h3ccc = (TH3D*)h3c->Clone( "3d");
	h3cc->GetZaxis()->SetRangeUser(875,950);

	h3c->GetZaxis()->SetRangeUser(a,b);
	TH1D *hzc = h3c->Project3D( "z2");
	zl = hzc->GetXaxis()->FindBin( a);
	zh = hzc->GetXaxis()->FindBin( b);
	zcc = hzc->Integral(zl,zh);

	TH1D *hzcc = h3cc->Project3D( "z3");
	sprintf( name, "M_{miss} for E_{#gamma} = %3d MeV\n", eg);
	hzcc->SetTitle( name);
	TH3D *h3ccc = (TH3D*)h3c->Clone( "3d");
	h3ccc->GetZaxis()->SetRangeUser(a,b);
	TH1D *hxcc = h3ccc->Project3D("x3");
	sprintf( name, " #theta*_{#pi} = %d - %d deg\n", theta_lo, theta_hi);
	hxcc->SetTitle( name);

	// zcc is the total number of carbon counts in the energy window, i.e.
	// y^{tot}_c.  It is NOT a function of theta.  It is a function of the
	// tagger energy bin, and the missing mass window (which should not chance).

	h3c->GetXaxis()->SetRangeUser(theta_lo,theta_hi);
	TH1D *hxc = h3c->Project3D("x2");
	xcc = hxc->Integral();

	// xcc is the total number of carbon counts in the energy window AND theta
	// window.  It is a function of theta, and obviously tagger energy bin and
	// missing mass window.

	TH1D *hyc = h3c->Project3D("y2");
	sprintf( name, " #theta_{#gamma#gamma} #geq %3d deg\n", (int) x1);
	hyc->SetTitle( name);
	yl = hyc->GetXaxis()->FindBin( x1);
	yh = hyc->GetXaxis()->FindBin( 180);
	ycc = hyc->Integral(x1,180);

	// ycc is the number of carbon counts in the energy window AND theta window
	// that pass the tgg cut, i.e. y_c(theta).

	// The efficiency factor f_c(theta) is then ycc/zcc

	effc = ac*ycc/zcc;

//	cout << zcc;
////	cout << "  " << xcc;
//	cout << "  " << ycc;
//	cout << "  " << effc;
//	cout << endl;

	if ( plot == 1) {
		TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 1000,
				500);
	   c1->SetGrid();
	   c1->GetFrame()->SetFillColor( 21);
	   c1->GetFrame()->SetBorderSize( 12);
//		c1->Divide(3,2);
		c1->Divide(3,1);

		c1->cd(1);
		hzcc->SetLineWidth(2);
		hzcc->SetFillColor(29);
		hzcc->Draw();

		max = hzcc->GetMaximum();
		max *= 1.10;
		hzcc->SetMaximum(max);

		l = new TLine( a, 0, a, max);
		l->SetLineColor(2);
		l->SetLineStyle(2);
		l->SetLineWidth(4);
		l->Draw();

		l = new TLine( b, 0, b, max);
		l->SetLineColor(2);
		l->SetLineStyle(2);
		l->SetLineWidth(4);
		l->Draw();

		Char_t contstr[32];
		sprintf( contstr, "%d\n", zcc);
		pl = new TPaveLabel(0.2,0.4,0.4,0.5, contstr, "NDC");
		pl->SetBorderSize(0);
		pl->SetFillColor(0);
		pl->SetTextSize(0.4);
		pl->Draw();

		c1->cd(2);
		hxcc->Draw();
		hxcc->SetLineWidth(2);
		hxcc->SetFillColor(29);

		max = hxcc->GetMaximum();
		max *= 1.10;
		hxcc->SetMaximum(max);

		l = new TLine( theta_lo, 0, theta_lo, max);
		l->SetLineColor(2);
		l->SetLineStyle(2);
		l->SetLineWidth(4);
		l->Draw();

		l = new TLine( theta_hi, 0, theta_hi, max);
		l->SetLineColor(2);
		l->SetLineStyle(2);
		l->SetLineWidth(4);
		l->Draw();

		sprintf( contstr, "%d\n", xcc);
		pl = new TPaveLabel(0.7,0.4,0.9,0.5, contstr, "NDC");
		pl->SetBorderSize(0);
		pl->SetFillColor(0);
		pl->SetTextSize(0.4);
		pl->Draw();

		c1->cd(3);
		hyc->Draw();
		hyc->SetLineWidth(2);
		hyc->SetFillColor(29);

		max = hyc->GetMaximum();
		max *= 1.10;
		hyc->SetMaximum(max);

		l = new TLine( x1, 0, x1, max);
		l->SetLineColor(2);
		l->SetLineStyle(2);
		l->SetLineWidth(4);
		l->Draw();

		pt = new TPaveText(0.3,0.4,0.5,0.55, "brNDC");
		pt->SetBorderSize(0);
		pt->SetFillColor(0);
		pt->SetTextSize(0.06);
		sprintf( contstr, "%d\n", ycc);
		pt->AddText( contstr);
		sprintf( contstr, "f_{c} = %4.1f %%\n", effc*100);
		pt->AddText( contstr);
		pt->Draw();

//		c1->cd(4);
//		h2p->Draw();

//		c1->cd(5);
//		h2c->Draw();

		sprintf( name, "plots/eps/sim_eff_%d_%d.eps", eg,
				(theta_lo + theta_hi)/2);
//		c1->Print( name);
	}

	return( effc);
}

void EmPi02D( Int_t pbin)
{
	Int_t eg;
	Double_t a, b, c;
	Char_t name[64];

	if ( pbin == 2) eg = 147;
	else if ( pbin == 3) eg = 153;
	else if ( pbin == 4) eg = 157;
	else if ( pbin == 5) eg = 163;
	else if ( pbin == 6) eg = 168;
	else {
		printf( "<E> CarbCont: Invalid Channel Number (2-6)\n");
		exit(-1);
	}

	a = Linear( 145, 918, 400, 933, eg);
	b = Linear( 145, 933, 400, 880, eg);
	c = Linear( 145, 945, 400, 1000, eg);

	sprintf( name, "histograms/MonteCarlo/MCH_p_masscorr_%d.root", eg);
	pdat = new TFile( name);
	h3p = (TH3D*) pdat->Get( "EmPi0_TGG_ThetaCM");
	TH2D *h2p = h3p->Project3D("xz1");
	h2p->SetMarkerColor( 4);

	sprintf( name, "histograms/MonteCarlo/MCH_c_masscorr_%d.root", eg);
	cdat = new TFile( name);
	h3c = (TH3D*)cdat->Get( "EmPi0_TGG_ThetaCM");
	TH2D *h2c = h3c->Project3D("xz2");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	h2p->GetXaxis()->SetRangeUser( 880, 950);
	h2p->Draw();
	h2c->Draw("same");

	l = new TLine( a, 0, a, 180);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	l = new TLine( b, 0, b, 180);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

	l = new TLine( c, 0, c, 180);
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->SetLineWidth(4);
	l->Draw();

//	sprintf( name, "plots/eps/sim_eff_%d_%d.eps", eg, (theta_lo + theta_hi)/2);
//	c1->Print( name);
}

void EmPi02DDat()
{
	TFile full( "histograms/ARH_full.root");
//	TFile empty( "histograms/ARH_empty.root");

	h3 = (TH3D*) full.Get( "EmPi0_TGG_ThetaCM_P");
	TH2D *h2 = h3->Project3D("xz1");
//	h2->SetMarkerColor( 4);

//	h3c = (TH3D*) empty.Get( "EmPi0_TGG_ThetaCM_R");
//	TH2D *h2c = h3c->Project3D("xz2");

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	h2->GetXaxis()->SetRangeUser( 880, 950);
	h2->Draw();
//	h2c->Draw("same");

//	sprintf( name, "plots/eps/sim_eff_%d_%d.eps", eg, (theta_lo + theta_hi)/2);
//	c1->Print( name);
}

void CheckSim2( Int_t eg = 150)
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, x3, a, b, fact;
	Char_t hist[32], cond[32], name[64];

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 700);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(1,2);

	strcpy( hist, "EmPi0");

	c1->cd(1);
	strcpy( cond, "5cm");
	sprintf( name, "histograms/MonteCarlo/MCH_p_%d_%s.root", eg, cond);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%d_%s.root", eg, cond);
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
	cout << " 5 cm target" << endl;
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

	c1->cd(2);
	strcpy( cond, "10cm");
	sprintf( name, "histograms/MonteCarlo/MCH_p_%d_%s.root", eg, cond);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%d_%s.root", eg, cond);
	carb = new TFile( name);

	g3p = (TH3D*) prot->Get( "EmPi0_TGG_ThetaCM");
	TH1D *g1 = g3p->Project3D("z3");
	g3c = (TH3D*) carb->Get( "EmPi0_TGG_ThetaCM");
	TH1D *g2 = g3c->Project3D("z4");

	sprintf( name, "E_{#gamma} = %5.1f MeV, %s target", eg, cond);
	g1->SetTitle( name);
	g1->SetLineWidth(2);
	g1->SetFillColor(2);
	g1->GetXaxis()->SetRangeUser(a,b);
	g1->Draw();

	ch1 = g1->Integral();
	ch2 = g2->Integral();
	fact = (double) ch1/ch2;

	chan1 = g1->GetXaxis()->FindBin( x1);
	chan2 = g1->GetXaxis()->FindBin( x2);
	chan3 = g1->GetXaxis()->FindBin( x3);
	intp = g1->Integral( chan1, chan2);
	intc = g2->Integral( chan3, chan1);
	Double_t ratio = intp/intc;
	cout << " 10 cm target" << endl;
	cout << intp << " " << intc << endl;
	cout << "ratio = " << ratio << endl;

	g2->SetFillStyle(3001);
	g2->SetFillColor(4);
	g2->SetLineWidth(2);
	g2->SetLineStyle(1);
	g2->Scale( fact);
	g2->Draw("same");

	TH1D *g3 = (TH1D*)g1->Clone();
	g3->Sumw2();
	g3->Add(g2, 1.0);
	g3->Draw("same");

	max = g3->GetMaximum();

	max *= 1.10;
	g1->SetMaximum(max);

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

	sprintf( name, "plots/pdf/EmPi0_tgtchk_%d.pdf", eg);
//	c1->Print( name);
}

void CheckSim3( Int_t eg = 150)
{
	Int_t col, ch1, ch2;
	Double_t max, x1, x2, x3, a, b, fact, scale;
	Char_t hist[32], cond[32], name[64];

	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 900, 500);
//	TCanvas *c1 = new TCanvas ( "c1", "Simulation Comparison", 200, 0, 700, 500);
   c1->SetFillColor( 42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide(2,1);

/*
	if ( eg == 150) scale = 0.2314;
	else if ( eg == 160) scale = 0.5162;
	else if ( eg == 170) scale = 0.8194;
*/

	if ( eg == 150) scale = 0.3658;
	else if ( eg == 160) scale = 0.8886;
	else if ( eg == 170) scale = 1.2218;

	strcpy( hist, "EmPi0");

	c1->cd(1);
	strcpy( cond, "5cm");
	sprintf( name, "histograms/MonteCarlo/MCH_p_%d_%s.root", eg, cond);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%d_%s.root", eg, cond);
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

	// This equalizes the number of counts in the two simulated histograms.
	ch1 = h1->Integral();
	ch2 = h2->Integral();
	fact = (double) ch1/ch2;
	cout << "fact = " << fact << endl;

	// This corrects the number of proton counts to agree with the data
	// as obtained by SimComp.C
	TH1D *h4 = (TH1D*)h1->Clone();
	h4->Scale( scale);
	h4->Draw();

	h2->SetFillStyle(3001);
	h2->SetFillColor(4);
	h2->SetLineWidth(2);
	h2->SetLineStyle(1);
	h2->Scale( fact);
	h2->Draw("same");

	TH1D *h3 = (TH1D*)h4->Clone();
	h3->Sumw2();
	h3->Add(h2, 1.0);
	h3->Draw("same");

	Int_t chan1 = h1->GetXaxis()->FindBin( x1);
	Int_t chan2 = h1->GetXaxis()->FindBin( x2);
	Int_t chan3 = h1->GetXaxis()->FindBin( x3);
	Double_t intp = h4->Integral( chan1, chan2);
	Double_t intc = h2->Integral( chan1, chan2);
	Double_t intt = h3->Integral( chan1, chan2);
	Double_t ratio = intc/intt;
	cout << " 5 cm target" << endl;
	cout << intp << " " << intc << " " << intt << endl;
	cout << "contamination = " << ratio << endl;

	max = h3->GetMaximum();

	max *= 1.10;
	h4->SetMaximum(max);

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

	c1->cd(2);
	strcpy( cond, "10cm");
	sprintf( name, "histograms/MonteCarlo/MCH_p_%d_%s.root", eg, cond);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%d_%s.root", eg, cond);
	carb = new TFile( name);

	g3p = (TH3D*) prot->Get( "EmPi0_TGG_ThetaCM");
	TH1D *g1 = g3p->Project3D("z3");
	g3c = (TH3D*) carb->Get( "EmPi0_TGG_ThetaCM");
	TH1D *g2 = g3c->Project3D("z4");

	sprintf( name, "E_{#gamma} = %5.1f MeV, %s target", eg, cond);
	g1->SetTitle( name);
	g1->SetLineWidth(2);
	g1->SetFillColor(2);
	g1->GetXaxis()->SetRangeUser(a,b);

	ch1 = g1->Integral();
	ch2 = g2->Integral();
	fact = (double) ch1/ch2;

	TH1D *g4 = (TH1D*)g1->Clone();
	g4->Scale( 2*scale);
	g4->Draw();

	g2->SetFillStyle(3001);
	g2->SetFillColor(4);
	g2->SetLineWidth(2);
	g2->SetLineStyle(1);
	g2->Scale( fact);
	g2->Draw("same");

	TH1D *g3 = (TH1D*)g4->Clone();
	g3->Sumw2();
	g3->Add(g2, 1.0);
	g3->Draw("same");

	intp = g4->Integral( chan1, chan2);
	intc = g2->Integral( chan1, chan2);
	intt = g3->Integral( chan1, chan2);
	ratio = intc/intt;
	cout << " 10 cm target" << endl;
	cout << intp << " " << intc << " " << intt << endl;
	cout << "contamination = " << ratio << endl;

	max = g3->GetMaximum();

	max *= 1.10;
	g4->SetMaximum(max);

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

//	sprintf( name, "plots/pdf/EmPi0_tgtchk_%d.pdf", eg);
	sprintf( name, "plots/eps/EmPi0_tgtchk_%d.eps", eg);
	c1->Print( name);
}

void CheckSim4( Int_t eg = 150)
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
	sprintf( name, "histograms/MonteCarlo/MCH_p_%d_2cm.root", eg);
	prot = new TFile( name);
	sprintf( name, "histograms/MonteCarlo/MCH_c_%d_2cm.root", eg);
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
