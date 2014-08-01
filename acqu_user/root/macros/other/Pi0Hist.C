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

typedef struct {
	Int_t mid;
	Double_t lo;
	Double_t hi;
} ThetaBins;
ThetaBins theta[10];

Double_t pa = 0.0833;

TFile full( "histograms/Pi0/Full/All_20MeV.root");
TH1D *hfsc = (TH1D*)full.Get( "SumScalers152to503");
TH1D *hdf = (TH1D*)full.Get( "SumScalers504to535");

TFile empty( "histograms/Pi0/Empty/All_20MeV.root");
TH1D *hesc = (TH1D*)empty.Get( "SumScalers152to503");
TH1D *hde = (TH1D*)empty.Get( "SumScalers504to535");

TH1D* fhist;
TH1D* ehist;
TH1D* shist;

TH1D* phist;
TH1D* rhist;
TH1D* subhist;

void Init()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

	theta[0].mid = 0;
	theta[0].lo = 0;
	theta[0].hi = 180;

	for ( UInt_t i = 1; i <= 9; i++) {
		theta[i].mid = 20*i - 10;
		theta[i].lo = 20*(i-1);
		theta[i].hi = 20*i;
	}
}

// This scale factor is for correcting the tagger scalers.  You must DIVIDE
// the tagger scalers by this amount!
// See Section x.x on page xx of the Tagger Manual.
//
// NOTE: The deadtime scalers must be in the pointer 'hdtsc'.
//
Double_t DeadTimeSF( TFile* file)
{
	Double_t f_tagg, f_tot, f_dead;
	TString dtscalers;

	TH1D* hdtsc;

	dtscalers = "SumScalers504to535";
	hdtsc = (TH1D*) file->Get( dtscalers);

	f_tagg = hdtsc->GetBinContent( 32)/hdtsc->GetBinContent( 31);
	f_tot = hdtsc->GetBinContent( 25)/hdtsc->GetBinContent( 26);
	f_dead = f_tagg/f_tot;

	return( f_dead);
}

//
//	HistSF
//
// This factor is for scaling the histogram belonging to the second set of
// scalers to that belonging to the first.
//
// No deadtime correction!
//
Double_t HistSF( TH1D* hsc1, TH1D* hsc2, UInt_t bin_lo, UInt_t bin_hi)
{
	Double_t c1, c2, factor;

	c1 = hsc1->Integral( bin_lo, bin_hi);
	c2 = hsc2->Integral( bin_lo, bin_hi);
	factor = c1/c2;

	return( factor);
}

void Pi0HistChan( TString hist = "MMiss", UInt_t i = 308, UInt_t j = 0,
		Bool_t Cut = kFALSE, Bool_t save = kFALSE)
{
	UInt_t rebin;
	Int_t th, dth;
	Double_t eg, deg;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;

	TggCutOffset = 3;

	TString name;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 800, 500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	th = theta[j].mid;
	dth = 20;

	x[0] = 0.55;
	x[1] = 0.80;
	y[0] = 0.75;
	y[1] = 0.85;
	if ( hist == "MMiss")
	{
		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);

		left = x1 - 30;
		right = x2 + 20;

		rebin = 1;
	}
	else if ( hist == "EMiss")
	{
//		x1 = -5;
		x1 = Linear( 144, -1, 190, -15, eg);
//		x2 = 25;
		x2 = Linear( 144, 5, 190, 20, eg);

		x3 = -1;


//		left = -25;
//		right = 100;
//		left = x1 - 10;
//		right = x2 + 25;
		left = -20;
		right = 35;

		rebin = 2;
	}
	else if ( hist == "TGG")
	{
		// Proton Opening Angle
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

		// 12-C Opening Angle
		q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

		x3 = -1;

		left = x2 - 45;
		right = 180;

		rebin = 4;

		x[0] = 0.15;
		x[1] = 0.40;
	}
	else if ( hist == "KECM")
	{
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		
		x1 = T_pi;

		left = 0;
		right = x1 + 20;

		rebin = 1;
	}
	else if ( hist == "M2gAron")
	{
		x1 = 115;
		x2 = 155;
		x3 = -1;


		left = x1 - 10;
		right = x2 + 10;

		rebin = 4;
	}

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);
	l->SetLineWidth( 2);

	Pi0Hist3D( i, j, hist, Cut);

	c1->cd( 1);
	name = Form( "full%d", i);
	TH1D *h1 = (TH1D*) fhist->Clone( name);
	h1->SetMarkerStyle( 21);
	h1->SetMarkerColor( 4);
	h1->SetLineWidth( 3);
	h1->SetLineColor( 4);
	h1->SetTitle();
	h1->GetXaxis()->SetTitleOffset( 1.1);
	h1->GetXaxis()->SetLabelSize( 0.03);
	h1->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "EMiss") name = "E_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	else if ( hist == "M2gAron") name = "m_{#gamma#gamma} (MeV)";
	h1->GetXaxis()->SetTitle( name);
	h1->Rebin( rebin);
	h1->GetXaxis()->SetRangeUser( left, right);
	h1->Draw();
	c1->Update();
	Double_t max = gPad->GetUymax();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 3);
	l2->SetLineWidth( 2);
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 3);
	l3->SetLineWidth( 2);
//	if ( hist == "MMiss") l3->Draw();

	l->Draw();

	name = Form( "empty%d", i);
	TH1D *h2 = (TH1D*) ehist->Clone( name);
	h2->SetMarkerStyle( 22);
	h2->SetMarkerColor( 2);
	h2->SetLineWidth( 3);
	h2->SetLineColor( 2);
	h2->Rebin( rebin);
	h2->Draw( "same");

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h1, "Full");
	tl->AddEntry( h2, "Empty");
	tl->Draw();

	name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.5);
	pl->SetTextFont( 132);
	pl->Draw();

	c1->cd( 2);
	name = Form( "subt%d", i);
	TH1D *h3 = (TH1D*) shist->Clone( name);
	h3->SetMarkerStyle( 20);
	h3->SetLineWidth( 3);
	h3->SetTitle();
	h3->GetXaxis()->SetTitleOffset( 1.1);
	h3->GetXaxis()->SetLabelSize( 0.03);
	h3->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "EMiss") name = "E_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	else if ( hist == "M2gAron") name = "m_{#gamma#gamma} (MeV)";
	h3->GetXaxis()->SetTitle( name);
	h3->Rebin( rebin);
	h3->GetXaxis()->SetRangeUser( left, right);
	h3->Draw();
	c1->Update();
	Double_t max = gPad->GetUymax();

	if ( j == 0) name = "#theta_{CM} = 0 - 180 deg";
	else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
	pl2 = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.5);
	pl2->SetTextFont( 132);
	pl2->Draw();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 3);
	l2->SetLineWidth( 2);
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 3);
	l3->SetLineWidth( 2);
//	if ( hist == "MMiss") l3->Draw();

	l->Draw();
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h3, "Subtracted");
	tl->Draw();

	if ( Cut == kFALSE)
		name = Form( "plots/Pi0/%s_chan%d_t%d", hist.Data(), i, th);
	else
		name = Form( "plots/Pi0/%s_chan%d_t%d_cut", hist.Data(), i, th);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void Pi0Hist4Eg( TString hist = "MMiss", UInt_t i = 308, Bool_t save = kFALSE)
{
	UInt_t i, j, rebin;
	Int_t th, dth;
	Double_t eg, deg;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;
	TString name;

	TggCutOffset = 3;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 1200, 600);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 5, 2);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	for ( j = 1; j <= 9; j++)
	{
		th = theta[j].mid;
		dth = 20;

		if ( hist == "MMiss")
		{
			x1 = Linear( 145, 933, 400, 880, eg);
			x2 = Linear( 145, 945, 400, 1000, eg);
			x3 = Linear( 145, 918, 400, 880, eg);

			left = x1 - 30;
			right = x2 + 20;

			rebin = 1;
		}
		else if ( hist == "TGG")
		{
			// Proton Opening Angle
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

			// 12-C Opening Angle
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

			x3 = -1;

			left = 0;
			right = 180;

			rebin = 4;
		}
		else if ( hist == "KECM")
		{
			// Proton
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x1 = T_pi;

			// 12-C
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x2 = T_pi;

			x3 = -1;

			left = 0;
			right = x2 + 20;

			rebin = 1;
		}

		l = new TLine( left, 0, right, 0);
		l->SetLineStyle( 2);

		Pi0Hist3D( i, j, hist);

		c1->cd( j);
		name = Form( "full%d", j);
		TH1D *h1 = (TH1D*) fhist->Clone( name);
		h1->SetMarkerStyle( 21);
		h1->SetMarkerColor( 2);
		h1->SetLineWidth( 3);
		h1->SetLineColor( 2);
		h1->SetTitle();
		h1->GetXaxis()->SetTitleOffset( 1.1);
		h1->GetXaxis()->SetLabelSize( 0.03);
		h1->GetXaxis()->CenterTitle();
		if ( hist == "MMiss") name = "M_{miss} (MeV)";
		else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
		else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
		h1->GetXaxis()->SetTitle( name);
		h1->Rebin( rebin);
		h1->GetXaxis()->SetRangeUser( left, right);
		Double_t max = h1->GetMaximum();
		max *= 1.4;
		h1->SetMaximum( max);
		h1->Draw();

		l1 = new TLine( x1, 0, x1, max);
		l1->SetLineStyle( 2);
		l1->Draw();
		l2 = new TLine( x2, 0, x2, max);
		l2->SetLineStyle( 2);
		l2->Draw();
		l3 = new TLine( x3, 0, x3, max);
		l3->SetLineStyle( 2);
//		l3->Draw();

		l->Draw();
		name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg);
		pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
		pl->SetFillStyle( 0);
		pl->SetBorderSize( 0);
		pl->SetTextSize( 0.5);
		pl->SetTextFont( 132);
		pl->Draw();
		if ( tbin == 0) name = "#theta_{CM} = 0 - 180 deg";
		else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
		pl2 = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
		pl2->SetFillStyle( 0);
		pl2->SetBorderSize( 0);
		pl2->SetTextSize( 0.5);
		pl2->SetTextFont( 132);
		if ( j == 6) pl2->Draw();

		name = Form( "empty%d", j);
		TH1D *h2 = (TH1D*) ehist->Clone( name);
		h2->SetMarkerStyle( 22);
		h2->SetMarkerColor( 4);
		h2->SetLineWidth( 3);
		h2->SetLineColor( 4);
		h2->Rebin( rebin);
		h2->Draw( "same");

		x[0] = 0.15;
		x[1] = 0.40;
		y[0] = 0.6;
		y[1] = 0.75;
		tl = new TLegend( x[0], y[0], x[1], y[1]);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.06);
		tl->SetTextFont( 132);
		tl->AddEntry( h1, "Full");
		tl->AddEntry( h2, "Empty");
		tl->Draw();
	}

	name = Form( "plots/Pi0_%s_chan%d-%d_t%d", hist.Data(), lo,
			lo-4, theta[tbin].mid);
//	name.Append( ".pdf");
	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void Pi0Hist4Th( TString hist = "MMiss", UInt_t i = 308, Bool_t save = kFALSE)
{
	UInt_t j, rebin;
	Int_t th, dth;
	Double_t eg, deg, max;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;
	TString name;

	TggCutOffset = 3;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 800, 800);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 3, 3);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	for ( j = 1; j <= 9; j++)
	{

		th = theta[j].mid;
		dth = 20;

		if ( hist == "MMiss")
		{
			x1 = Linear( 145, 933, 400, 880, eg);
			x2 = Linear( 145, 945, 400, 1000, eg);
			x3 = Linear( 145, 918, 400, 880, eg);

			left = x1 - 30;
			right = x2 + 20;

			rebin = 1;
		}
		else if ( hist == "TGG")
		{
			// Proton Opening Angle
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

			// 12-C Opening Angle
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

			x3 = -1;

			left = 0;
			right = 180;

			rebin = 4;
		}
		else if ( hist == "KECM")
		{
			// Proton
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x1 = T_pi;

			// 12-C
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
			T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
			x2 = T_pi;

			x3 = -1;

			left = 0;
			right = x2 + 20;

			rebin = 1;
		}

		l = new TLine( left, 0, right, 0);
		l->SetLineStyle( 2);

		Pi0Hist3D( i, j, hist);

		c1->cd( j);
		name = Form( "full%d", j);
		TH1D *h1 = (TH1D*) fhist->Clone( name);
		h1->SetMarkerStyle( 21);
		h1->SetMarkerColor( 4);
		h1->SetLineWidth( 3);
		h1->SetLineColor( 4);
		h1->SetTitle();
		h1->GetXaxis()->SetTitleOffset( 1.1);
		h1->GetXaxis()->SetLabelSize( 0.03);
		h1->GetXaxis()->CenterTitle();
		if ( hist == "MMiss") name = "M_{miss} (MeV)";
		else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
		else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
		h1->GetXaxis()->SetTitle( name);
		h1->Rebin( rebin);
		h1->GetXaxis()->SetRangeUser( left, right);
		h1->Draw();
		c1->Update();
		max = gPad->GetUymax();

		l1 = new TLine( x1, 0, x1, max);
		l1->SetLineStyle( 2);
		l1->Draw();
		l2 = new TLine( x2, 0, x2, max);
		l2->SetLineStyle( 2);
		l2->Draw();
		l3 = new TLine( x3, 0, x3, max);
		l3->SetLineStyle( 2);
//		l3->Draw();

		l->Draw();
		name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg);
		pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
		pl->SetFillStyle( 0);
		pl->SetBorderSize( 0);
		pl->SetTextSize( 0.5);
		pl->SetTextFont( 132);
		if ( j == 1) pl->Draw();
		name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
		pl2 = new TPaveLabel( 0.15, 0.75, 0.45, 0.85, name, "NDC");
		pl2->SetFillStyle( 0);
		pl2->SetBorderSize( 0);
		pl2->SetTextSize( 0.5);
		pl2->SetTextFont( 132);
		pl2->Draw();

		name = Form( "empty%d", j);
		TH1D *h2 = (TH1D*) ehist->Clone( name);
		h2->SetMarkerStyle( 22);
		h2->SetMarkerColor( 2);
		h2->SetLineWidth( 3);
		h2->SetLineColor( 2);
		h2->Rebin( rebin);
		h2->Draw( "same");

		x[0] = 0.15;
		x[1] = 0.40;
		y[0] = 0.6;
		y[1] = 0.75;
		tl = new TLegend( x[0], y[0], x[1], y[1]);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.06);
		tl->SetTextFont( 132);
		tl->AddEntry( h1, "Full");
		tl->AddEntry( h2, "Empty");
		if ( j == 1) tl->Draw();

	}

	name = Form( "plots/Pi0/Pi0_%s_chan%d", hist.Data(), i);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void Pi0HistSimComp( TString tgt = "subt", TString hist = "MMiss",
		UInt_t i = 308, UInt_t j = 0, Bool_t save = kFALSE)
{
	UInt_t tbin_lo, tbin_hi, rebin;
	Int_t th, dth;
	Double_t eg, deg, max;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;
	Char_t tchar;
	TString name;

	TggCutOffset = 3;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 700, 500);
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	th = theta[j].mid;
	dth = 20;

	if ( hist == "MMiss")
	{
		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);

		left = x1 - 30;
		right = x2 + 20;

		rebin = 1;
	}
	else if ( hist == "EMiss")
	{
		x1 = Linear( 144, -1, 190, -15, eg);
		x2 = Linear( 144, 5, 190, 20, eg);

		x3 = -1;

		left = -20;
		right = 35;

		rebin = 2;
	}
	else if ( hist == "TGG")
	{
		// Proton Opening Angle
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;
		cout << "tgg_min = " << x1;
//		x1 -= TggCutOffset;
		cout << " tgg_cut = " << x1;

		// 12-C Opening Angle
		q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

		cout << " tgg_c12 = " << x2;
		cout << endl;

		x3 = 180;

		left = x2 - 45;
		right = 180;

		rebin = 4;
	}
	else if ( hist == "KECM")
	{
		if ( tgt == "subt")
			q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		else if ( tgt == "empty")
			q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);

		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;

		x1 = 0;
		x2 = T_pi + 20;
		x3 = T_pi;

		left = 0;
		right = x2;

		rebin = 1;
	}

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);

	Pi0Hist3D( i, j, hist);

	l->Draw();

	name = Form( "%s%d", tgt.Data(), i);
	if ( tgt == "subt") TH1D *h3 = (TH1D*) shist->Clone( name);
	else if ( tgt == "empty") TH1D *h3 = (TH1D*) ehist->Clone( name);
	h3->SetMarkerStyle( 20);
	h3->SetMarkerSize( 1.5);
	h3->SetLineWidth( 4);
	h3->SetTitle();
	h3->GetXaxis()->SetTitleOffset( 1.2);
	h3->GetXaxis()->SetLabelSize( 0.03);
	h3->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "EMiss") name = "E_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	h3->GetXaxis()->SetTitle( name);
	h3->Rebin( rebin);
	h3->GetXaxis()->SetRangeUser( left, right);
	if ( ( hist == "TGG") && ( i > 306)) h3->SetMinimum( -50);
	h3->Draw();
	c1->Update();
	max = gPad->GetUymax();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 2);
//	if ( hist == "MMiss") l2->Draw();
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 2);
//	if ( hist == "MMiss") l3->Draw();

	l->Draw();

	if ( ( hist == "KECM") || ( hist == "EMiss"))
	{
		x[0] = 0.60;
		x[1] = 0.85;
	}
	else
	{
		x[0] = 0.15;
		x[1] = 0.40;
	}
	y[0] = 0.7;
	y[1] = 0.85;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	if ( tgt == "subt") name = "Subtracted";
	else if ( tgt == "empty") name = "Windows";
	else if ( tgt == "full") name = "Unsubtracted";
	tl->AddEntry( h3, name);

	name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.5);
	pl->SetTextFont( 132);
	pl->Draw();

	if ( j == 0) name = "#theta_{CM} = 0 - 180 deg";
	else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
	pl2 = new TPaveLabel( 0.55, 0.90, 0.85, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.5);
	pl2->SetTextFont( 132);
	pl2->Draw();

	UInt_t bin_lo;
	UInt_t bin_hi;

	if ( ( hist == "MMiss") && ( tgt == "empty")) {
		bin_lo = h3->GetXaxis()->FindBin( left);
		bin_hi = h3->GetXaxis()->FindBin( right);
	}
	else {
		bin_lo = h3->GetXaxis()->FindBin( x1);
		bin_hi = h3->GetXaxis()->FindBin( x2);
	}
	Double_t cts1 = h3->Integral( bin_lo, bin_hi);

	if ( tgt == "subt") tchar = 'p';
	else if ( tgt == "empty") tchar = 'w';
	name = Form( "histograms/MonteCarlo/pi0/LH2_10cm/pi0_%c_chan%d.root", tchar,
			i);
	sim = new TFile( name);
	name = "THR_" + hist + "P_v_ThetaCMP";
	TH2D *h2sim = (TH2D*)sim->Get( name);
	tbin_lo = h2sim->GetXaxis()->FindBin( theta[j].lo);
	tbin_hi = h2sim->GetXaxis()->FindBin( theta[j].hi)-1;
	h2sim->GetXaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *h1sim = h2sim->ProjectionY( "simhist");
	h1sim->Rebin( rebin);

	Double_t cts2 = h1sim->Integral( bin_lo, bin_hi);

	Double_t fact = cts1/cts2;

	h1sim->Scale( fact);
	h1sim->SetFillStyle( 3004);
	h1sim->SetFillColor( 16);
	h1sim->SetLineWidth( 2);
	h1sim->Draw( "same");

	tl->AddEntry( h1sim, "Monte Carlo");
	tl->Draw();

	name = Form( "plots/Pi0/Pi0_%c_%s_SimComp_chan%d_t%d", tchar,
			hist.Data(), i, th);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}

void Pi0Hist3D( UInt_t chan, UInt_t tbin, TString var_y,
		Bool_t Cut = kFALSE)
{
	UInt_t cbin, tbin_lo, tbin_hi;
	TString name;
	TString var_x, var_z;
	TString hprompt, hrand;

	var_x = "TChan";
	var_z = "ThetaCM";

	if ( Cut == kTRUE) {
		var_x += "Cut1";
		var_y += "Cut1";
		var_z += "Cut1";
	}

	hprompt = var_z + "P_v_" + var_y + "P_v_" + var_x + "P";
	hrand = var_z + "R_v_" + var_y + "R_v_" + var_x + "R";

	// Target FULL
	// Prompt
	TH3D *hP = (TH3D*)full.Get( hprompt);
	cbin = hP->GetXaxis()->FindBin( chan);
	hP->GetXaxis()->SetRange( cbin, cbin);
	tbin_lo = hP->GetZaxis()->FindBin( theta[tbin].lo);
	tbin_hi = hP->GetZaxis()->FindBin( theta[tbin].hi)-1;
	hP->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hPf_y = hP->Project3D( "y");

	// Random
	TH3D *hR = (TH3D*)full.Get( hrand);
	hR->GetXaxis()->SetRange( cbin, cbin);
	hR->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hRf_y = hR->Project3D( "y2");
	hRf_y->Scale( pa);

	// Subtracted
	hSf = (TH1D*) hPf_y->Clone( "FullSubt");
	hSf->Sumw2();
	hSf->Add( hRf_y, -1.0);

	// Target EMPTY
	// Prompt
	TH3D *hP = (TH3D*)empty.Get( hprompt);
	hP->GetXaxis()->SetRange( cbin, cbin);
	hP->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hPe_y = hP->Project3D( "y");

	// Random
	TH3D *hR = (TH3D*)empty.Get( hrand);
	hR->GetXaxis()->SetRange( cbin, cbin);
	hR->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hRe_y = hR->Project3D( "y2");
	hRe_y->Scale( pa);

	// Subtracted
	hSe = (TH1D*) hPe_y->Clone( "EmptySubt");
	hSe->Sumw2();
	hSe->Add( hRe_y, -1.0);

	// COMPLETELY SUBTRACTED
	Double_t f_tagg_f, f_tot_f, f_dead_f, fscal;
	Double_t f_tagg_e, f_tot_e, f_dead_e, escal;
	Double_t factor;

	f_tagg_f = hdf->GetBinContent( 32)/hdf->GetBinContent( 31);
	f_tot_f = hdf->GetBinContent( 25)/hdf->GetBinContent( 26);
	f_dead_f = f_tagg_f/f_tot_f;
	fscal = hfsc->GetBinContent( chan)*f_dead_f;
//	fscal = hfsc->GetBinContent( chan);

	f_tagg_e = hde->GetBinContent( 32)/hde->GetBinContent( 31);
	f_tot_e = hde->GetBinContent( 25)/hde->GetBinContent( 26);
	f_dead_e = f_tagg_e/f_tot_e;
	escal = hesc->GetBinContent( chan)*f_dead_e;
//	escal = hesc->GetBinContent( chan);

	factor = escal/fscal;
	hSf->Scale( factor);
	hS = (TH1D*) hSf->Clone( "Subtracted");
	hS->Add( hSe, -1.0);

	fhist = (TH1D*) hSf->Clone( "Full");
	ehist = (TH1D*) hSe->Clone( "Empty");
	shist = (TH1D*) hS->Clone( "Subt");
}

void M2gHist( UInt_t chan, UInt_t rebin = 1)
{
	UInt_t cbin;
	Double_t factor;
	TString name;

	TFile* file;

	cbin = chan+1;

	// Target FULL
	// Prompt
	TH2D *hP = (TH2D*)full.Get( "THR_M2gTaggP_v_TChanM2gP");
	hP->GetXaxis()->SetRange( cbin, cbin);
	TH1D *hPf_y = hP->ProjectionY( "M2gP");

	// Random
	TH2D *hR = (TH2D*)full.Get( "THR_M2gTaggR_v_TChanM2gR");
	hR->GetXaxis()->SetRange( cbin, cbin);
	TH1D *hRf_y = hR->ProjectionY( "M2gR");
	hRf_y->Scale( 0.5);

	// Subtracted
	hSf = (TH1D*) hPf_y->Clone( "FullSubt");
	hSf->Sumw2();
	hSf->Add( hRf_y, -1.0);

	// Target EMPTY
	// Prompt
	TH2D *hP = (TH2D*)empty.Get( "THR_M2gTaggP_v_TChanM2gP");
	hP->GetXaxis()->SetRange( cbin, cbin);
	TH1D *hPe_y = hP->ProjectionY( "M2gP");

	// Random
	TH2D *hR = (TH2D*)empty.Get( "THR_M2gTaggR_v_TChanM2gR");
	hR->GetXaxis()->SetRange( cbin, cbin);
	TH1D *hRe_y = hR->ProjectionY( "M2gR");
	hRe_y->Scale( 0.5);

	// Subtracted
	hSe = (TH1D*) hPe_y->Clone( "EmptySubt");
	hSe->Sumw2();
	hSe->Add( hRe_y, -1.0);

	// COMPLETELY SUBTRACTED

	// Dead-time correction
	file = full;
	Double_t f_dead_f = DeadTimeSF( file);
	file = empty;
	Double_t f_dead_e = DeadTimeSF( file);

	// Scale full to empty
	factor = HistSF( hesc, hfsc, cbin, cbin)*f_dead_f/f_dead_e;

	hSf->Scale( factor);
	hS = (TH1D*) hSf->Clone( "Subtracted");
	hS->Add( hSe, -1.0);

	hSf->Rebin( rebin);
	hSe->Rebin( rebin);
	hS->Rebin( rebin);

	fhist = (TH1D*) hSf->Clone( "Full");
	ehist = (TH1D*) hSe->Clone( "Empty");
	shist = (TH1D*) hS->Clone( "Subt");
}

void M2gPlot( UInt_t chan, Bool_t save = kFALSE)
{

	Double_t eg, deg, left, right;
	Double_t x1, x2, max;
	TString name;

// Set up Canvas
	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 700, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);
	gStyle->SetOptStat( 0);

	eg = tcd[chan].energy;
	deg = tcd[chan].denergy;

	name = Form( "Egamma = %5.1f +/- %3.1f MeV", eg, deg);
	cout << name << endl;

	left = 60;
	right = 160;
	x1 = kMPI0_MEV - 20;
	x2 = kMPI0_MEV + 20;

// Get Histograms
	M2gHist( chan, 4);

	// Full
	name = Form( "full_%d", chan);
	if ( hf = (TH1D*)(gROOT->FindObject( name)) ) delete hf;
	TH1D *hf = (TH1D*) fhist->Clone( name);

	// Empty
	name = Form( "empty_%d", chan);
	if ( he = (TH1D*)(gROOT->FindObject( name)) ) delete he;
	TH1D *he = (TH1D*) ehist->Clone( name);

	// Subtracted
	name = Form( "subt_%d", chan);
	if ( hs = (TH1D*)(gROOT->FindObject( name)) ) delete hs;
	TH1D *hs = (TH1D*) shist->Clone( name);

	c1->cd( 1);
	// Full
	hf->SetLineWidth( 3);
	hf->SetLineColor( 4);
	name = Form( "2#gamma Invariant Mass for E_{#gamma} = %5.1f #pm %3.1f MeV",
			eg, deg);
	hf->SetTitle( name);
	hf->GetXaxis()->SetTitleOffset( 1.1);
	hf->GetXaxis()->SetLabelSize( 0.03);
	hf->GetXaxis()->CenterTitle();
	hf->GetXaxis()->SetTitle( "M_{#gamma#gamma} (MeV)");
	hf->GetXaxis()->SetRangeUser( left, right);
	hf->Draw();
	c1->Update();

	max = gPad->GetUymax();
	cout << max << endl;
	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 3);
	l2->SetLineWidth( 2);
	l2->Draw();

	// Empty
	he->SetTitle();
	he->SetLineWidth( 3);
	he->SetLineColor( 2);
	he->Draw( "same");

	Double_t x[0] = 0.15;
	Double_t x[1] = 0.40;
	Double_t y[0] = 0.75;
	Double_t y[1] = 0.85;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( hf, "Full");
	tl->AddEntry( he, "Empty");
	tl->Draw();

	c1->cd(2);
	// Subtracted
	hs->SetTitle();
	hs->SetMarkerStyle( 21);
	hs->SetMarkerColor( 1);
	hs->SetLineWidth( 3);
	hs->SetLineColor( 1);
	hs->SetTitle();
	hs->Draw( "E");
	hs->GetXaxis()->SetRangeUser( left, right);
	c1->Update();

	max = gPad->GetUymax();
	cout << max << endl;
	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);
	l->SetLineWidth( 2);
	l->Draw();
	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 3);
	l2->SetLineWidth( 2);
	l2->Draw();

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( hs, "Subtracted");
	tl->Draw();

	name = Form( "plots/Pi0/m2g_%d.pdf", chan);
	if ( save == kTRUE) c1->Print( name);
}

void M2gChan( UInt_t chan = 250, UInt_t rebin = 1, TString tgt = "full")
{
	UInt_t left, right, bin_lo, bin_hi, rebin;
	Double_t eg, deg;
	TString name;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 200, 20, 700, 400);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	eg = tcd[chan].energy;
	deg = tcd[chan].denergy;

	name = Form( "Egamma = %5.1f +/- %3.1f MeV", eg, deg);
	cout << name << endl;

	left = 0;
	right = 200;

	// Prompt
	name = "THR_M2gTaggP_v_TChanM2gP";
	if ( tgt == "full") TH2D *hp2d = (TH2D*)full.Get( name);
	else if ( tgt == "empty") TH2D *hp2d = (TH2D*)empty.Get( name);
	bin_lo = hp2d->GetXaxis()->FindBin( chan);
	bin_hi = hp2d->GetXaxis()->FindBin( chan);
	hp2d->GetXaxis()->SetRange( bin_lo, bin_hi);
	name = Form( "M2gP_r%d", rebin);
	TH1D *hp = hp2d->ProjectionY( name);
	hp->Rebin( rebin);

	// Random
	name = "THR_M2gTaggR_v_TChanM2gR";
	if ( tgt == "full") TH2D *hr2d = (TH2D*)full.Get( name);
	else if ( tgt == "empty") TH2D *hr2d = (TH2D*)empty.Get( name);
	hr2d->GetXaxis()->SetRange( bin_lo, bin_hi);
	name = Form( "M2gR_r%d", rebin);
	TH1D *hr = hr2d->ProjectionY( name);
	hr->Rebin( rebin);

	c1->cd( 1);
	// Prompt
	hp->SetLineWidth( 3);
	hp->SetLineColor( 4);
	name = Form( "2#gamma Invariant Mass for E_{#gamma} = %5.1f #pm %3.1f MeV",
			eg, deg);
	hp->SetTitle( name);
	hp->GetXaxis()->SetTitleOffset( 1.1);
	hp->GetXaxis()->SetLabelSize( 0.03);
	hp->GetXaxis()->CenterTitle();
	hp->GetXaxis()->SetTitle( "M_{#gamma#gamma} (MeV)");
	hp->GetXaxis()->SetRangeUser( left, right);
	hp->DrawCopy();

	// Random
	hr->Scale( 0.5);
	hr->SetLineWidth( 3);
	hr->SetLineColor( 2);
	hr->Draw( "same");

	// Subtracted
	hp->Sumw2();
	hp->Add( hr, -1.0);
	hp->SetMarkerStyle( 21);
	hp->SetMarkerColor( 1);
	hp->SetLineWidth( 3);
	hp->SetLineColor( 1);
	c1->cd(2);
	hp->SetTitle();
	hp->Draw( "E");
	hp->GetXaxis()->SetRangeUser( left, right);

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);
	l->Draw();
}

void Pi0Hist3D_PR( TString subt, UInt_t chan, UInt_t tbin, TString var_y,
		Bool_t Cut = kFALSE)
{
	UInt_t cbin, tbin_lo, tbin_hi;
	TString name;
	TString var_x, var_z;
	TString hprompt, hrand;

	var_x = "TChan";
	var_z = "ThetaCM";

	if ( Cut == kTRUE) {
		var_x += "Cut1";
		var_y += "Cut1";
		var_z += "Cut1";
	}

	hprompt = var_z + "P_v_" + var_y + "P_v_" + var_x + "P";
	hrand = var_z + "R_v_" + var_y + "R_v_" + var_x + "R";

	// Prompt
	if ( subt == "full") TH3D *hP = (TH3D*)full.Get( hprompt);
	else if ( subt == "empty") TH3D *hP = (TH3D*)empty.Get( hprompt);
	cbin = hP->GetXaxis()->FindBin( chan);
	hP->GetXaxis()->SetRange( cbin, cbin);
	tbin_lo = hP->GetZaxis()->FindBin( theta[tbin].lo);
	tbin_hi = hP->GetZaxis()->FindBin( theta[tbin].hi)-1;
	hP->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hP_y = hP->Project3D( "y");

	// Random
	if ( subt == "full") TH3D *hR = (TH3D*)full.Get( hrand);
	else if ( subt == "empty") TH3D *hR = (TH3D*)empty.Get( hrand);
	hR->GetXaxis()->SetRange( cbin, cbin);
	hR->GetZaxis()->SetRange( tbin_lo, tbin_hi);
	TH1D *hR_y = hR->Project3D( "y2");
	hR_y->Scale( pa);

	// Subtracted
	hS = (TH1D*) hP_y->Clone( "Subt");
	hS->Sumw2();
	hS->Add( hR_y, -1.0);

	phist = (TH1D*) hP_y->Clone( "Full");
	rhist = (TH1D*) hR_y->Clone( "Empty");
	subhist = (TH1D*) hS->Clone( "Subtracted");
}

void EMissCheck( TString subt = "full", TString hist = "EMiss",
		UInt_t i = 300, UInt_t j = 0, Bool_t Cut = kFALSE,
		Bool_t save = kFALSE)
{
	UInt_t rebin;
	Int_t th, dth;
	Double_t eg, deg;
	Double_t x1, x2, x3, left, right;
	Double_t x[2], y[2];
	Double_t q_pi, T_pi;
	Double_t TggCutOffset;

	TString name;

	gStyle->SetOptStat( 0);

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	TCanvas *c1 = new TCanvas ( "c1", "Subtraction", 300, 20, 800, 500);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->Divide( 2, 1);

	eg = tcd[i].energy;
	deg = tcd[i].denergy;

	th = theta[j].mid;
	dth = 20;

	x[0] = 0.55;
	x[1] = 0.80;
	y[0] = 0.75;
	y[1] = 0.85;
	if ( hist == "MMiss")
	{
		x[0] = 0.25;
		x[1] = 0.50;

		x1 = Linear( 145, 933, 400, 880, eg);
		x2 = Linear( 145, 945, 400, 1000, eg);
		x3 = Linear( 145, 918, 400, 880, eg);

//		left = x1 - 100;
		left = 800;
		right = x2 + 20;

		rebin = 1;
	}
	else if ( hist == "EMiss")
	{
		x1 = -5;
		x2 = 25;
		x3 = -1;

		left = -50;
		right = 200;
//		left = -25;
//		right = 100;
//		left = -10;
//		right = 25;

		rebin = 1;
	}
	else if ( hist == "TGG")
	{
		// Proton Opening Angle
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x1 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R - TggCutOffset;

		// 12-C Opening Angle
		q_pi = qp_thcm( eg, kM_C12_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		x2 = Tgg_Min( T_pi, kMPI0_MEV)/kD2R;

		x3 = -1;

		left = x2 - 45;
		right = 180;

		rebin = 4;

		x[0] = 0.15;
		x[1] = 0.40;
	}
	else if ( hist == "KECM")
	{
		q_pi = qp_thcm( eg, kMP_MEV, theta[j].lo, kMPI0_MEV);
		T_pi = Energy( q_pi, kMPI0_MEV) - kMPI0_MEV;
		
		x1 = T_pi;

		left = 0;
		right = x1 + 20;

		rebin = 1;
	}

	l = new TLine( left, 0, right, 0);
	l->SetLineStyle( 2);
	l->SetLineWidth( 2);

	Pi0Hist3D_PR( subt, i, j, hist, Cut);

	c1->cd( 1);
	name = Form( "prompt%d", i);
	TH1D *h1 = (TH1D*) phist->Clone( name);
	h1->SetMarkerStyle( 21);
	h1->SetMarkerColor( 4);
	h1->SetLineWidth( 3);
	h1->SetLineColor( 4);
	h1->SetTitle();
	h1->GetXaxis()->SetTitleOffset( 1.1);
	h1->GetXaxis()->SetLabelSize( 0.03);
	h1->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "EMiss") name = "E_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	h1->GetXaxis()->SetTitle( name);
	h1->Rebin( rebin);
	h1->GetXaxis()->SetRangeUser( left, right);
	h1->Draw();
	c1->Update();
	Double_t max = gPad->GetUymax();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 3);
	l2->SetLineWidth( 2);
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 3);
	l3->SetLineWidth( 2);
//	if ( hist == "MMiss") l3->Draw();

//	l->Draw();

	name = Form( "random%d", i);
	TH1D *h2 = (TH1D*) rhist->Clone( name);
	h2->SetMarkerStyle( 22);
	h2->SetMarkerColor( 2);
	h2->SetLineWidth( 3);
	h2->SetLineColor( 2);
	h2->Rebin( rebin);
	h2->Draw( "same");

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h1, "Prompt");
	tl->AddEntry( h2, "Random");
	tl->Draw();

	name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", eg, deg);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.5);
	pl->SetTextFont( 132);
	pl->Draw();

	c1->cd( 2);
	name = Form( "subtracted%d", i);
	TH1D *h3 = (TH1D*) subhist->Clone( name);
	h3->SetMarkerStyle( 20);
	h3->SetLineWidth( 3);
	h3->SetTitle();
	h3->GetXaxis()->SetTitleOffset( 1.1);
	h3->GetXaxis()->SetLabelSize( 0.03);
	h3->GetXaxis()->CenterTitle();
	if ( hist == "MMiss") name = "M_{miss} (MeV)";
	else if ( hist == "EMiss") name = "E_{miss} (MeV)";
	else if ( hist == "TGG") name = "#theta_{#gamma#gamma} (deg)";
	else if ( hist == "KECM") name = "CM T_{#pi} (MeV)";
	h3->GetXaxis()->SetTitle( name);
	h3->Rebin( rebin);
	h3->GetXaxis()->SetRangeUser( left, right);
	h3->Draw();
	c1->Update();
	Double_t max = gPad->GetUymax();

	if ( j == 0) name = "#theta_{CM} = 0 - 180 deg";
	else name = Form( "#theta_{CM} = %d #pm %d deg", th, dth/2);
	pl2 = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.5);
	pl2->SetTextFont( 132);
	pl2->Draw();

	l1 = new TLine( x1, 0, x1, max);
	l1->SetLineStyle( 3);
	l1->SetLineWidth( 2);
	l1->Draw();
	l2 = new TLine( x2, 0, x2, max);
	l2->SetLineStyle( 3);
	l2->SetLineWidth( 2);
	l2->Draw();
	l3 = new TLine( x3, 0, x3, max);
	l3->SetLineStyle( 3);
	l3->SetLineWidth( 2);
//	if ( hist == "MMiss") l3->Draw();

	l->Draw();

	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.06);
	tl->SetTextFont( 132);
	tl->AddEntry( h3, "Subtracted");
	tl->Draw();

	if ( Cut == kTRUE)
		name = Form( "plots/Pi0/CheckCut_%s_%s_chan%d_t%d", subt.Data(),
				hist.Data(), i, th);
	else
		name = Form( "plots/Pi0/Check_%s_%s_chan%d_t%d", subt.Data(),
				hist.Data(), i, th);
	name.Append( ".pdf");
//	name.Append( ".eps");
	if ( save == kTRUE) c1->Print( name);
}
