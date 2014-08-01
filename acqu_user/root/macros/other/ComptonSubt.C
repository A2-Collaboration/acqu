gROOT->Reset();

TFile datafile("histograms/Compton/Full/Full.root");

typedef struct {
	Int_t egamma;
	Double_t energy;
	Double_t denergy;
} TChanData;
TChanData tcd[352];

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
} EBins;
EBins tch[11];

Int_t tchan[] = { 286, 282, 277, 273, 269, 265, 261, 257, 253, 249, 245};

Double_t pa;

TH1D* hist;

void Init()
{
	UInt_t i;

	pa = -0.08333;

	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");

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
	}
}

void PhotonMmissHist( UInt_t ebin, UInt_t tbin)
{
	UInt_t cbin_lo, cbin_hi, tbin_lo, tbin_hi;
	TString name;

	// Random Subtraction
	name = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP";
	TH3D *hP = (TH3D*)datafile.Get( name);
	name = "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR";
	TH3D *hR = (TH3D*)datafile.Get( name);
	hP->Sumw2();
	hP->Add( hR, -pa);

	// Tagger Channel Bin
	cbin_lo = hP->GetXaxis()->FindBin( tch[ebin].lo);
	cbin_hi = hP->GetXaxis()->FindBin( tch[ebin].hi);
	hP->GetXaxis()->SetRange( cbin_lo, cbin_hi);

	// Theta CM Bin
	tbin_lo = hP->GetYaxis()->FindBin( theta[tbin].lo);
	tbin_hi = hP->GetYaxis()->FindBin( theta[tbin].hi)-1;
	hP->GetYaxis()->SetRange( tbin_lo, tbin_hi);

	TH1D *hP_z = hP->Project3D( "z");

	hist = (TH1D*) hP_z->Clone( "Subt");
}

void PhotonMmissBin( Int_t ebin = 8, Int_t tbin = 0, Int_t rebin = 2,
		Bool_t save = kFALSE)
{
	Int_t cbin_lo, cbin_hi, tbin_lo, tbin_hi;
	Double_t left, right;
	Double_t x[2], y[2];
	TString name;

	Double_t open_cut = 5;

	left = 800;
	right = 1100;
//	left = 900;
//	right = 1049;

	name = Form( "Tagger energies %5.1f - %5.1f MeV", tcd[tch[ebin].hi].energy,
			tcd[tch[ebin].lo].energy);
	cout << name << endl;
	if ( tbin == 0) name = "All Theta";
	else name = Form( "Photon CM Theta %d +/- 10 deg", theta[tbin].mid);
	cout << name << endl;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Random Subtraction", 200, 0, 700, 500);
	c1->SetFillColor( 18);

	gStyle->SetOptStat( 0);

	PhotonMmissHist( ebin, tbin);

	name = Form( "Subt%d", tch[ebin].mid);
	TH1D *hP_z = (TH1D*) hist->Clone( name);
	name = "";
	hP_z->SetTitle( name);

	hP_z->GetXaxis()->SetTitleOffset( 1.2);
	hP_z->GetXaxis()->SetTitle( "M_{miss} (MeV)");
	hP_z->GetXaxis()->SetLabelSize( 0.03);
	hP_z->GetXaxis()->CenterTitle();

	hP_z->Rebin( rebin);
	hP_z->Draw();

	Double_t max = hP_z->GetMaximum();
	max *= 1.4;
	hP_z->SetMaximum( max);

	l1 = new TLine( left, 0, right, 0);
	l1->SetLineStyle(1);
	l1->SetLineWidth(1);
	l1->Draw();

	name = Form( "E_{#gamma} = %5.1f - %5.1f MeV", tcd[tch[ebin].hi].energy,
			tcd[tch[ebin].lo].energy);
	pl = new TPaveLabel( 0.15, 0.90, 0.45, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.4);
	pl->SetTextFont( 132);
	pl->Draw();
	if ( tbin == 0) name = "#theta_{CM} = 0 - 180 deg";
	else name = Form( "#theta_{CM} = %3d #pm 10 deg", theta[tbin].mid);
	pl2 = new TPaveLabel( 0.55, 0.90, 0.85, 1.00, name, "NDC");
	pl2->SetFillStyle( 0);
	pl2->SetBorderSize( 0);
	pl2->SetTextSize( 0.4);
	pl2->SetTextFont( 132);
	pl2->Draw();

	y[0] = 0.75;
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
		name = Form( "#theta^{CM}_{#gamma} = %d #pm 10 deg", theta[tbin].mid);
		pt->AddText( name);
	}
	name = Form( "#theta_{open} < %d deg", (int) open_cut);
	pt->AddText( name);
//	name = Form( "M_{miss} < %d MeV", (int) msep);
//	pt->AddText( name);
	pt->Draw();

	name = Form( "plots/PhotonMmissData_%d_%d_%d.pdf", (int) open_cut,
			tch[ebin].mid, theta[tbin].mid);
	if ( save == kTRUE) c1->Print( name);
}

void PhotonMmiss( Int_t cut = 0, Bool_t save = kFALSE)
{
	TString h1, name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Data (With Random Subtraction)", 200, 10, 700, 500);
   c1->GetFrame()->SetFillColor( 21);
   c1->GetFrame()->SetBorderSize( 12);
	c1->SetFillColor(18);

	h1 = "COMP_PhotonMmiss";
	if ( cut == 1) h1 = "COMP_PhotonMmissProt";
	else if ( cut == 2) h1 = "COMP_PhotonMmissProtOA";
	name = h1 + "P";
	TH1F *MMiss_P = (TH1F*)datafile.Get( name);
	name = h1 + "R";
	TH1F *MMiss_R = (TH1F*)datafile.Get( name);

	name = "Proton Missing Mass";
	if ( cut == 1) name.Append( " with Proton Cut");
	else if ( cut == 2) name.Append( " with Proton OA Cut");

	MMiss_P->Sumw2();
	MMiss_P->Add( MMiss_R, -0.0833);
	MMiss_P->GetXaxis()->SetTitleOffset( 1.1);
	MMiss_P->GetXaxis()->SetTitle("M_{miss} (MeV)");
	MMiss_P->GetXaxis()->SetLabelSize( 0.03);
	MMiss_P->GetXaxis()->CenterTitle();
	MMiss_P->SetMarkerStyle( 21);
	MMiss_P->SetMarkerColor( 4);
	MMiss_P->SetMarkerSize( 0.5);
	MMiss_P->SetLineWidth( 2);
	MMiss_P->SetLineColor( 4);
	MMiss_P->Draw();

	l1 = new TLine( 800, 0, 1100, 0);
	l1->SetLineStyle(2);
	l1->SetLineWidth(2);
	l1->Draw();

	name = "plots/PhotonMmiss";
	if ( cut == 1) name.Append( "Prot");
	else if ( cut == 2) name.Append( "ProtOA");
	name.Append( ".pdf");
	if ( save == kTRUE) c1->Print( name);
}

void PhotonMissTChanCut( UInt_t flag = 0, Bool_t save = kFALSE)
{
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas( "c1", "Random Subtraction", 200, 10, 700, 500);
	c1->SetFillColor( 18);

	name = "COMP_PhotonMmissProtOAP_v_TChanPhotProtOAP";
	if ( flag == 1) name = "COMP_PhotonEmissProtOAP_v_TChanPhotProtOAP";
	TH2D *hP = (TH2D*)datafile.Get( name);

	name = "COMP_PhotonMmissProtOAR_v_TChanPhotProtOAR";
	if ( flag == 1) name = "COMP_PhotonEmissProtOAR_v_TChanPhotProtOAR";
	TH2D *hR = (TH2D*)datafile.Get( name);

	hP->Sumw2();
	hP->Add( hR, -0.0833);

	name = "Missing Mass vs. Tagger Channel";
	if ( flag == 1) name = "Missing Energy vs. Tagger Channel";
	hP->SetTitle( name);

	hP->GetXaxis()->SetTitleOffset( 1.2);
	hP->GetXaxis()->SetTitle( "Tagger Channel");
	hP->GetXaxis()->SetLabelSize( 0.03);
	hP->GetXaxis()->CenterTitle();

	hP->GetYaxis()->SetTitleOffset( 1.2);
	hP->GetYaxis()->SetTitle( "M_{miss} (MeV)");
	if ( flag == 1) hP->GetYaxis()->SetTitle( "E_{miss} (MeV)");
	hP->GetYaxis()->SetLabelSize( 0.03);
	hP->GetYaxis()->CenterTitle();

	hP->Draw("colz");

	name = "plots/Mmiss_TChan_Cut.pdf";
	if ( flag == 1) name = "plots/Emiss_TChan_Cut.pdf";
	if ( save == kTRUE) c1->Print( name);
}

void PhotonMissOA( Int_t flag = 0, Bool_t save = kFALSE)
{
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas("c1", "Random Subtraction", 800, 0, 700, 500);
	c1->SetFillColor(18);

	name = "COMP_PhotonMmissProtP_v_PhotProtOAP";
	if ( flag == 1) name = "COMP_PhotonEmissProtP_v_PhotProtOAP";
	TH2D *hP = (TH2D*)datafile.Get( name);

	name = "COMP_PhotonMmissProtR_v_PhotProtOAR";
	if ( flag == 1) name = "COMP_PhotonEmissProtR_v_PhotProtOAR";
	TH2D *hR = (TH2D*)datafile.Get( name);

	hP->Sumw2();
	hP->Add( hR, -0.0833);

	name = "Missing Mass vs. Opening Angle";
	if ( flag == 1) name = "Missing Energy vs. Opening Angle";
	hP->SetTitle( name);

	hP->GetXaxis()->SetTitleOffset( 1.2);
	hP->GetXaxis()->SetTitle("Opening Angle (deg)");
	hP->GetXaxis()->SetLabelSize( 0.03);
	hP->GetXaxis()->CenterTitle();

	hP->GetYaxis()->SetTitleOffset( 1.2);
	hP->GetYaxis()->SetTitle("M_{miss} (MeV)");
	if ( flag == 1) hP->GetYaxis()->SetTitle("E_{miss} (MeV)");
	hP->GetYaxis()->SetLabelSize( 0.03);
	hP->GetYaxis()->CenterTitle();

	hP->Draw( "colz");

	name = "plots/PhotonMmiss_OpenA.pdf";
	if ( flag == 1) name = "plots/PhotonEmiss_OpenA.pdf";
	if ( save == kTRUE) c1->Print( name);
}

void PhotonMmissThetaCM( Bool_t save = kFALSE)
{
	TString name;

	if ( !gROOT->GetListOfCanvases()->IsEmpty()) delete c1;
	c1 = new TCanvas("c1", "Random Subtraction", 200, 10, 700, 500);
	c1->SetFillColor( 18);

	name = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP";
	TH3D *hP = (TH3D*)datafile.Get( name);
	name = "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR";
	TH3D *hR = (TH3D*)datafile.Get( name);

	hP->Sumw2();
	hP->Add( hR, -0.0833);

	// This is currently 262 - 279 MeV
	hP->GetXaxis()->SetRangeUser( 253, 260);
	TH2D *hP_zy = hP->Project3D( "zy");

	name = "Missing Mass vs. Photon Scattering Angle";
	hP_zy->SetTitle( name);

	hP_zy->GetXaxis()->SetTitleOffset( 1.2);
	hP_zy->GetXaxis()->SetTitle( "#theta_{#gamma} (deg)");
	hP_zy->GetXaxis()->SetLabelSize( 0.03);
	hP_zy->GetXaxis()->CenterTitle();

	hP_zy->GetYaxis()->SetTitleOffset( 1.2);
	hP_zy->GetYaxis()->SetTitle( "M_{miss} (MeV)");
	hP_zy->GetYaxis()->SetLabelSize( 0.03);
	hP_zy->GetYaxis()->CenterTitle();

	hP_zy->Draw( "colz");

	name = "plots/Mmiss_Theta_Cut.pdf";
	if ( save == kTRUE) c1->Print( name);
}
