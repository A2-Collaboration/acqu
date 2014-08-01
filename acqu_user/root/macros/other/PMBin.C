gROOT->Reset();

typedef struct {

	Int_t egamma;
	Double_t energy;
	Double_t denergy;

} TData;
TData tcd[352];

typedef struct {
	Int_t mid;
	Double_t lo;
	Double_t hi;
} ThetaBins;
ThetaBins theta[10];

typedef struct {
	Int_t mid;
	Int_t lo;
	Int_t hi;
	Int_t eg;
} EBins;
EBins tch[11];

Int_t tchan[] = { 286, 282, 277, 273, 269, 265, 261, 257, 253, 249, 245};

Double_t fact[11][7];

TFile datafile("histograms/Compton/Full/Full.root");

TString hist;

hist = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP";
TH3D *hP = (TH3D*)datafile.Get( hist);
hist = "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR";
TH3D *hR = (TH3D*)datafile.Get( hist);

TH1D* hP_z;

TH1D* hcomp;
TH1D* hback;
TH1D* hsum;

Int_t rebin = 5;

Int_t left, right;

void Init()
{
	Int_t i;
	Double_t junk, pa;
	TString name;

	pa = -0.0833;

	i = 0;
	name = "includes/txs_ratio.dat";
	ifstream rFile( name);
	if ( !rFile.is_open()) 
	{
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		exit( -1);
	}

	for ( i = 0; i < 11; i++)
		rFile >> junk >> fact[i][0] >> fact[i][1] >> fact[i][2]
			>> fact[i][3] >> fact[i][4] >> fact[i][5] >> fact[i][6];
	rFile.close();

	left = 800;
	right = 1200;

	gROOT->ProcessLine( ".L root/macros/ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

	hP->Sumw2();
	hP->Add( hR, -pa);

	theta[0].mid = 0;
	theta[0].lo = 0;
	theta[0].hi = 180;

	for ( i = 1; i <= 9; i++) {
		theta[i].mid = 20*i - 10;
		theta[i].lo = 20*(i-1);
		theta[i].hi = 20*i;
	}

	for ( i = 0; i <= 10; i++) {
		tch[i].mid = tchan[i];
		tch[i].lo = tchan[i] - 1;
		tch[i].hi = tchan[i] + 1;
		tch[i].eg = 200 + i*10;
	}
}

void CompEbin( UInt_t ebin, Bool_t save = kFALSE)
{
	UInt_t tbin;
	Int_t dc, sc;
	Int_t lo, hi;
	Double_t factor;
	Double_t x[2], y[2];
	Double_t elo, ehi;
	TString name;

	gStyle->SetOptStat( 0);
	Init();
	
	TCanvas *canv;

	if ( !(canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		TCanvas* canv = new TCanvas( "c1", "Test Canvas", 400, 0, 1000, 700); 
	canv->Divide( 2, 2);

	UInt_t i = 1;
	for ( tbin = 5; tbin <= 8; tbin++ )
	{
		canv->cd( i++);

		GetData( ebin, tbin);
		name = Form( "data%d-%d", ebin, tbin);
		TH1D *dhist = (TH1D*) hP_z->Clone( name);
		dhist->Rebin( rebin);
		lo = dhist->GetXaxis()->FindBin( 920);
		hi = dhist->GetXaxis()->FindBin( 1050);
		dc = dhist->Integral( lo, hi);
		dhist->SetMarkerStyle( 21);
		dhist->SetMarkerColor( 4);
		dhist->SetMarkerSize( 1.0);
		dhist->SetLineWidth( 4);
		dhist->SetLineColor( 4);
		dhist->GetXaxis()->SetTitle( "M_{miss} (MeV)");
		dhist->GetXaxis()->CenterTitle();
		elo = tcd[tch[ebin].hi].energy - tcd[tch[ebin].hi].denergy;
		ehi = tcd[tch[ebin].lo].energy + tcd[tch[ebin].lo].denergy;
		name = Form( "E_{#gamma} = %5.1f - %5.1f MeV", elo, ehi);
		if ( i == 2) dhist->SetTitle( name);
		else dhist->SetTitle();
		Double_t max = dhist->GetMaximum();
		dhist->SetMaximum( 1.5*max);
		dhist->Draw();

		GetSim( ebin, tbin, 5);
		name = Form( "hc%d-%d", ebin, tbin);
		TH1D *hc = (TH1D*) hcomp->Clone( name);
		hc->Rebin( rebin);
		name = Form( "hb%d-%d", ebin, tbin);
		TH1D *hb = (TH1D*) hback->Clone( name);
		hb->Rebin( rebin);
		name = Form( "hs%d-%d", ebin, tbin);
		TH1D *hs = (TH1D*) hsum->Clone( name);
		hs->Rebin( rebin);
		sc = hs->Integral( lo, hi);
		factor = (double) dc/sc;

		hc->Scale( factor);
		hb->Scale( factor);
		hs->Scale( factor);

		hs->SetLineWidth( 4);
		hs->Draw( "same");

		hc->SetFillStyle( 3008);
		hc->SetFillColor( 13);
		hc->SetLineWidth( 0);
		hc->Draw( "same");

		hb->SetFillStyle( 3001);
		hb->SetLineWidth( 0);
		hb->SetFillColor( 17);
		hb->Draw( "same");

		x[0] = 0.2;
		x[1] = 0.5;
		y[0] = 0.6;
		y[1] = 0.9;
		tl = new TLegend( x[0], y[0], x[1], y[1]);
		tl->SetFillStyle( 0);
		tl->SetBorderSize( 0);
		tl->SetTextSize( 0.05);
		tl->SetTextFont( 132);
		tl->AddEntry( dhist, "Data", "p");
		tl->AddEntry( hc, "Compton");
		tl->AddEntry( hb, "#pi^{0}");
		tl->AddEntry( hs, "Sum");
		if ( i == 2) tl->Draw();

		x[0] = 0.6;
		x[1] = 0.8;
		y[0] = 0.8;
		y[1] = 0.9;
		name = Form( "#theta^{CM}_{#gamma'} = %d #pm 10 deg", theta[tbin].mid);
		pl = new TPaveLabel( x[0], y[0], x[1], y[1], name, "NDC");
		pl->SetFillStyle( 0);
		pl->SetBorderSize( 0);
		pl->SetTextSize( 0.45);
		pl->SetTextFont( 132);
		pl->Draw();
	}
	canv->Update();

	name = Form( "plots/PMBin_%d.pdf", tch[ebin].eg);
	if ( save == kTRUE) c1->Print( name);
}

void CompOneBin( UInt_t ebin, UInt_t tbin)
{
	UInt_t lo, hi;
	Double_t dc, sc;
	Double_t factor;
	Double_t x[2], y[2];
	Double_t elo, ehi;
	TString name;

	gStyle->SetOptStat( 0);
	Init();

	TCanvas *canv;

	if ( !(canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		TCanvas* canv = new TCanvas( "c1", "Test Canvas", 400, 0, 700, 500); 

	GetData( ebin, tbin);
	name = Form( "data%d-%d", ebin, tbin);
	TH1D *dhist = (TH1D*) hP_z->Clone( name);
	dhist->Rebin( rebin);
	lo = dhist->GetXaxis()->FindBin( 920);
	hi = dhist->GetXaxis()->FindBin( 980);
	dc = dhist->Integral( lo, hi);
	dhist->SetMarkerStyle( 21);
	dhist->SetMarkerColor( 4);
	dhist->SetMarkerSize( 1.0);
	dhist->SetLineWidth( 4);
	dhist->SetLineColor( 4);
	dhist->GetXaxis()->SetTitle( "M_{miss} (MeV)");
	dhist->GetXaxis()->CenterTitle();
	elo = tcd[tch[ebin].hi].energy - tcd[tch[ebin].hi].denergy;
	ehi = tcd[tch[ebin].lo].energy + tcd[tch[ebin].lo].denergy;
	name = Form( "E_{#gamma} = %5.1f - %5.1f MeV", elo, ehi);
	dhist->SetTitle( name);
	Double_t max = dhist->GetMaximum();
	dhist->SetMaximum( 1.5*max);
	dhist->Draw();

	GetSim( ebin, tbin, 5);
	name = Form( "hc%d-%d", ebin, tbin);
	TH1D *hc = (TH1D*) hcomp->Clone( name);
	hc->Rebin( rebin);
	name = Form( "hb%d-%d", ebin, tbin);
	TH1D *hb = (TH1D*) hback->Clone( name);
	hb->Rebin( rebin);
	name = Form( "hs%d-%d", ebin, tbin);
	TH1D *hs = (TH1D*) hsum->Clone( name);
	hs->Rebin( rebin);
	sc = hs->Integral( lo, hi);
	factor = (double) dc/sc;

	hc->Scale( factor);
	hb->Scale( factor);
	hs->Scale( factor);

	hs->SetLineWidth( 4);
	hs->Draw( "same");

	hc->SetFillStyle( 3008);
	hc->SetFillColor( 13);
	hc->SetLineWidth( 0);
	hc->Draw( "same");

	hb->SetFillStyle( 3001);
	hb->SetLineWidth( 0);
	hb->SetFillColor( 17);
	hb->Draw( "same");

	x[0] = 0.2;
	x[1] = 0.5;
	y[0] = 0.6;
	y[1] = 0.9;
	tl = new TLegend( x[0], y[0], x[1], y[1]);
	tl->SetFillStyle( 0);
	tl->SetBorderSize( 0);
	tl->SetTextSize( 0.05);
	tl->SetTextFont( 132);
	tl->AddEntry( dhist, "Data", "p");
	tl->AddEntry( hc, "Compton");
	tl->AddEntry( hb, "#pi^{0}");
	tl->AddEntry( hs, "Sum");
	tl->Draw();

	x[0] = 0.6;
	x[1] = 0.8;
	y[0] = 0.8;
	y[1] = 0.9;
	if ( tbin != 0)
		name = Form( "#theta^{CM}_{#gamma'} = %d #pm 10 deg", theta[tbin].mid);
	else
		name = Form( "#theta^{CM}_{#gamma'} = 0-180 deg", theta[tbin].mid);
	pl = new TPaveLabel( x[0], y[0], x[1], y[1], name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.45);
	pl->SetTextFont( 132);
	pl->Draw();

	lo = dhist->GetXaxis()->FindBin( 900);
	hi = dhist->GetXaxis()->FindBin( 950);
	dc = dhist->Integral( lo, hi);
	cout << " counts = " << dc << endl;
}

void GetData( UInt_t ebin, UInt_t tbin)
{
	Int_t bin_lo, bin_hi;
	TString name;

	if ( tbin == 0) {
		bin_lo = hP->GetYaxis()->FindBin( theta[tbin].lo);
		bin_hi = hP->GetYaxis()->FindBin( theta[tbin].hi);
	}
	else {
		bin_lo = hP->GetYaxis()->FindBin( theta[tbin].lo);
		bin_hi = hP->GetYaxis()->FindBin( theta[tbin].hi)-1;
	}

	hP->GetZaxis()->SetRangeUser( left, right);
	hP->GetYaxis()->SetRange( bin_lo, bin_hi);
	hP->GetXaxis()->SetRangeUser( tch[ebin].lo, tch[ebin].hi);

	name = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP_z";
	if ( (hP_z = (TH1D*)(gROOT->FindObject( name)))) delete hP_z;
	hP_z = (TH1D*) hP->Project3D( "z");
}

void GetSim( Int_t ebin, Int_t tbin, Double_t open_cut)
{
	Int_t bin_lo, bin_hi;
	TString name;

	name = Form( "histograms/MonteCarlo/compton/LH2_10cm/compton_p_%d.root",
			tch[ebin].eg);
	comp = new TFile( name);

	name = Form( "histograms/MonteCarlo/compton/LH2_10cm/pi0_p_%d.root",
			tch[ebin].eg);
	back = new TFile( name);

	// Compton
	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";
	h3comp = (TH3D*) comp->Get( name);
	h3comp->GetXaxis()->SetRangeUser( 0, open_cut);
	if ( tbin == 0) {
		bin_lo = h3comp->GetYaxis()->FindBin( theta[tbin].lo);
		bin_hi = h3comp->GetYaxis()->FindBin( theta[tbin].hi);
	}
	else {
		bin_lo = h3comp->GetYaxis()->FindBin( theta[tbin].lo);
		bin_hi = h3comp->GetYaxis()->FindBin( theta[tbin].hi)-1;
	}
	h3comp->GetYaxis()->SetRange( bin_lo, bin_hi);
	h3comp->GetZaxis()->SetRangeUser( left, right);
	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP_z";
	if ( (hcomp = (TH1D*)(gROOT->FindObject( name)))) delete hcomp;
	hcomp = (TH1D*) h3comp->Project3D( "z");

	// pi0 p background
	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP";
	h3back = (TH3D*) back->Get( name);
	h3back->GetXaxis()->SetRangeUser( 0, open_cut);
	h3back->GetYaxis()->SetRange( bin_lo, bin_hi);
	h3back->GetZaxis()->SetRangeUser( left, right);
	name = "PhotonMmissProtP_v_PhotonThetaCMProtP_v_PhotProtOAP_z2";
	if ( (hback = (TH1D*)(gROOT->FindObject( name)))) delete hback;
	hback = (TH1D*) h3back->Project3D( "z2");
	hback->Scale( fact[ebin][0]);

	hsum = (TH1D*) hcomp->Clone( "simsum");
	hsum->Add( hcomp, hback, 1, 1); 
}
