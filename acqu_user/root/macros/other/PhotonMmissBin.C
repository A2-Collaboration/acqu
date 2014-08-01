gROOT->Reset();

TFile inFile("histograms/Full/Full_100118.root");

typedef struct {

	Int_t egamma;
	Double_t energy;
	Double_t denergy;

} TData;

TData beam[352];

Int_t tch[] = { 286, 282, 277, 273, 269, 265, 261, 257, 253, 249, 245};
Int_t th[] = { 0, 10, 30, 50, 70, 90, 110, 130, 150, 170};

TCanvas *c1 = new TCanvas( "c1", "Random Subtraction", 200, 0, 1000, 700);
c1->SetFillColor( 18);

TH3D* hP;
TH3D* hR;

void Init()
{
	TString name;

	name = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP";
	hP = (TH3D*)inFile.Get( name);

	name = "PhotonMmissProtOAR_v_PhotonThetaCMProtOAR_v_TChanPhotProtOAR";
	hR = (TH3D*)inFile.Get( name);

	hP->Sumw2();
	hP->Add( hR, -0.0833);

	ReadTagEng( "xs/tageng855.dat");
}

void PMPlot( UInt_t tbin)
{
	UInt_t ebin;

	TString name;

	Double_t open_cut = 5;

	c1->Divide( 4, 3);

	for ( ebin = 0; ebin <= 10; ebin++)
	{
		c1->cd( ebin+1);
		PhotonMmissBin( ebin, tbin, kFALSE);
	}
	name = Form( "plots/PhotonMmissData_oa%d_t%d.pdf", (int) open_cut,
			(int) th[tbin]);
	c1->Print( name);
}

void PhotonMmissBin( UInt_t ebin = 8, UInt_t tbin = 0, Bool_t save = kFALSE)
{
	Int_t chan_lo, chan_hi;
	Int_t theta_lo, theta_hi;
	Int_t rebin;
	Double_t left, right;
	Double_t x[2], y[2];
	TString hist, name;

	rebin = 5;

	gStyle->SetOptStat( 0);

	Double_t open_cut = 5;

	chan_lo = tch[ebin] - 1;
	chan_hi = tch[ebin] + 1;

	if ( tbin != 0) {
		theta_lo = th[tbin] - 10;
		theta_hi = th[tbin] + 10;
	}
	else {
		theta_lo = 0;
		theta_hi = 180;
	}

	left = 800;
	right = 1100;

	name = Form( "Tagger energies %5.1f - %5.1f MeV", beam[chan_hi].energy,
			beam[chan_lo].energy);
	cout << name << endl;
	name = Form( "Photon CM Theta %d - %d deg", theta_lo, theta_hi);
	cout << name << endl;

	hP->GetZaxis()->SetRangeUser( left, right);
	hP->GetYaxis()->SetRangeUser( theta_lo, theta_hi);
	hP->GetXaxis()->SetRangeUser( chan_lo, chan_hi);

	TH1D* hP_z;
	name = "PhotonMmissProtOAP_v_PhotonThetaCMProtOAP_v_TChanPhotProtOAP_z";
	if ( (hP_z = (TH1D*)(gROOT->FindObject( name)))) delete hP_z;

	hP_z = (TH1D*) hP->Project3D( "z");
	TH1D *hdraw = (TH1D*) hP_z->Clone( "clone");

	name = "";
	hdraw->SetTitle( name);

	hdraw->GetXaxis()->SetTitleOffset( 1.2);
	hdraw->GetXaxis()->SetTitle( "M_{miss} (MeV)");
	hdraw->GetXaxis()->SetLabelSize( 0.03);
	hdraw->GetXaxis()->CenterTitle();

	hdraw->Rebin( rebin);

	hdraw->Draw();

	Double_t max = hdraw->GetMaximum();
	max *= 1.4;
	hdraw->SetMaximum( max);

	l1 = new TLine( left, 0, right, 0);
	l1->SetLineStyle(1);
	l1->SetLineWidth(1);
	l1->Draw();

//	name = Form( "E_{#gamma} = %5.1f - %5.1f MeV", beam[chan_hi].energy,
//			beam[chan_lo].energy);
	name = Form( "#theta^{CM}_{#gamma'} = %d #pm 10 deg", th[tbin]);
	pl = new TPaveLabel( 0.15, 0.90, 0.40, 1.00, name, "NDC");
	pl->SetFillStyle( 0);
	pl->SetBorderSize( 0);
	pl->SetTextSize( 0.45);
	pl->SetTextFont( 132);
	if ( ebin == 0) pl->Draw();

	y[0] = 0.80;
	y[1] = 0.92;
	x[0] = 0.2;
	x[1] = 0.5;
	pt = new TPaveText( x[0], y[0], x[1], y[1], "NDC");
	pt->SetBorderSize( 0);
	pt->SetFillStyle( 0);
	pt->SetTextSize( 0.045);
	pt->SetTextFont( 132);
	pt->SetTextAlign( 12);
	name = Form( "E_{#gamma} = %5.1f - %5.1f MeV", beam[chan_hi].energy,
		beam[chan_lo].energy);
	pt->AddText( name);
//	name = Form( "#theta_{open} < %d deg", (int) open_cut);
//	pt->AddText( name);
//	name = Form( "M_{miss} < %d MeV", (int) msep);
//	pt->AddText( name);
	pt->Draw();

	Int_t egamma = (int) ( (beam[chan_hi].energy + beam[chan_lo].energy)/2
			+ 0.5);
	name = Form( "plots/PhotonMmissData_%d_%d_%d.pdf", (int) open_cut, egamma,
			(int) th[tbin]);
	if ( save == kTRUE) c1->Print( name);
}

void ReadTagEng( TString te_file)
{
	UInt_t i;
	Double_t eff, deff;

	ifstream TEfile( te_file);
	if ( !TEfile.is_open()) {
		cout << "Error opening file ";
		cout << te_file;
		cout << endl;
		break;
	}
	while( !TEfile.eof()) {
		TEfile >> i >> eff >> deff;
		beam[i].energy = eff;
		beam[i].denergy = deff/2;
		beam[i].egamma = (int)(eff + 0.5);
	}
	TEfile.close();
}
