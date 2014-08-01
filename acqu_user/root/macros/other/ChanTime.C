gROOT->Reset();

typedef struct {
	Int_t egamma;
	Double_t energy;
	Double_t denergy;
} TData;
TData tcd[352];

TFile full( "histograms/Pi0/Full/All_20MeV.root");
TFile empty( "histograms/Pi0/Empty/All_20MeV.root");

Int_t broken_ladd_ch[] = { -1};
Int_t prompt[2], acc[4];

Double_t par[5];
Double_t red_chisq;

prompt[0] = 95;
prompt[1] = 105;
acc[0] = 20;
acc[1] = 80;
acc[2] = 120;
acc[3] = 180;

// Full
TH2D *h2d_f = (TH2D*)full.Get( "THR_Pi0TaggTime_v_TChanHit");
TH1D *hits_f = (TH1D*)full.Get( "THR_TChanHit");
TH1D *scal_f = (TH1D*)full.Get( "SumScalers152to503");

// Empty
TH2D *h2d_e = (TH2D*)empty.Get( "THR_Pi0TaggTime_v_TChanHit");
TH1D *hits_e = (TH1D*)empty.Get( "THR_TChanHit");
TH1D *scal_e = (TH1D*)empty.Get( "SumScalers152to503");

void Init()
{
	gROOT->ProcessLine( ".L ReadParams.C");
	ReadTagEng( "xs/tageng855.dat");
}

void ChanAll( TString subt = "full", Int_t min = 0, Int_t max = 351)
{
	Int_t chan, flag;
	Int_t* lch;
	Double_t ratio[352];
	TString name;

	TCanvas *c1 = new TCanvas ( "c1", "Prompt-Random", 20, 20, 700, 500);

	for ( chan = min; chan <= max; chan++) {
		lch = broken_ladd_ch;
		flag = 0;
		while ( *lch != -1) {
			if ( chan == *lch) flag = 1;
			lch++;
		}
		if ( flag == 0) {
			ratio[chan] = TimeFit( subt, chan);
			cout << "chan = " << chan << "  subtraction ratio = " << ratio[chan]
				<< endl;
			c1->WaitPrimitive("");
		}
		else {
			ratio[chan] = 0;
		}
	}

	name = Form( "xs/chan_subt/chan_subt_%s_new.out", subt.Data());
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( chan = min; chan <= max; chan++) {
			outFile << chan << "  " << ratio[chan] << endl ;
	}
	outFile.close();
}

Double_t TimeFit( TString subt, Int_t chan)
{
	Int_t bin, i;
	Double_t x[6], y[6], A[3], ratio;
	Double_t fpar[5];
	Double_t max, max1, max2;
	TString name;

	x[0] = acc[0];
	x[1] = acc[1];
	x[2] = prompt[0];
	x[3] = prompt[1];
	x[4] = acc[2];
	x[5] = acc[3];

	TH1D* h1;
	TH2D* h2;
	TF1* f1;
	TF1* f2;
	TF1* fall;

	if ( subt == "full") h2 = h2d_f;
	else if ( subt == "empty") h2 = h2d_e;

	if ( ( chan >= 0) && ( chan <= 351)) {
		name = Form( "Channel %d", chan);
		bin = h2->GetXaxis()->FindBin( chan);
		h1 = h2->ProjectionY( name, bin, bin);
	}
	else {
		cout << "Error: channels must be 0-351" << endl;
		break;
	}

	f1 = new TF1( "f1", "pol1", x[0], x[1]);
	f2 = new TF1( "f2", "gaus", x[2], x[3]);

	h1->Fit( "f1", "NQR");
	f1->GetParameters( &fpar[0]);

	h1->Fit( "f2", "NQR+");
	f2->GetParameters( &fpar[2]);

	fall = new TF1( "fall", "pol1(0) + gaus(2)", x[0], x[5]);
	fall->SetParameters( fpar);

	h1->Fit( "fall", "QR+");
	cout << fpar[0] << " " << fpar[1] << endl;

	Int_t bins = h1->GetXaxis()->FindBin( x[5]) - h1->GetXaxis()->FindBin( x[0]);
	red_chisq = fall->GetChisquare()/(bins-7);

	cout << red_chisq << endl;

	for ( i = 0; i <= 5; i++) y[i] = fpar[0] + fpar[1]*x[i];

	A[0] = 0.5*(y[0] + y[1])*(x[1] - x[0]);
	A[1] = 0.5*(y[2] + y[3])*(x[3] - x[2]);
	A[2] = 0.5*(y[4] + y[5])*(x[5] - x[4]);

	ratio = A[1]/(A[0] + A[2]);

/*
	name = Form( "%6.4f %6.4f   %6.4f %6.4f   %6.4f %6.4f\n", x[0], x[1], x[2],
			x[3], x[4], x[5]);
	cout << name << endl;

	name = Form( "%6.4f %6.4f   %6.4f %6.4f   %6.4f %6.4f\n", y[0], y[1], y[2],
			y[3], y[4], y[5]);
	cout << name << endl;

	name = Form( "%6.4f %6.4f %6.4f   %6.4f\n", A[0], A[1], A[2], ratio);
	cout << name << endl;
*/

	max = h1->GetMaximum();

	max *= 1.1;

	h1->SetMaximum( max);

	l1 = new TLine( x[0], 0., x[0], max);
	l1->SetLineWidth(2);
	l1->Draw();

	l2 = new TLine( x[1], 0., x[1], max);
	l2->SetLineWidth(2);
	l2->Draw();

	l3 = new TLine( x[0], y[0], x[1], y[1]);
	l3->SetLineWidth(2);
	l3->Draw();

	l4 = new TLine( x[2], 0., x[2], max);
	l4->SetLineWidth(2);
	l4->Draw();

	l5 = new TLine( x[3], 0., x[3], max);
	l5->SetLineWidth(2);
	l5->Draw();

	l6 = new TLine( x[2], y[2], x[3], y[3]);
	l6->SetLineWidth(2);
	l6->Draw();

	l7 = new TLine( x[4], 0., x[4], max);
	l7->SetLineWidth(2);
	l7->Draw();

	l8 = new TLine( x[5], 0., x[5], max);
	l8->SetLineWidth(2);
	l8->Draw();

	l9 = new TLine( x[4], y[4], x[5], y[5]);
	l9->SetLineWidth(2);
	l9->Draw();

	max1 = max*0.87;
	max2 = max*0.93;

	name = Form( "Subtraction Ratio = %6.4f", ratio);
	pt = new TPaveText( 20, max1, 90, max2);
	pt->AddText( name);
	pt->Draw();

	gPad->SetLogy();

	return ratio;
}

Double_t TimeFitBin( Char_t subt[8], Int_t tbin)
{
	Int_t bin_lo, bin_hi;
	Int_t tag_chan_lo[] = { 256, 252, 248, 244, 240, 236, 232};
	Int_t tch_lo, tch_hi;
	Double_t param[2], ratio, max, max1, max2;
	Char_t string[80];

	if ( ( tbin < 2) || ( tbin > 6)) {
		cout << "<E> bin number ";
		cout << tbin;
		cout << " is not valid.  Must be 2-6.";
		cout << endl;
		break;
	}

	tch_lo = tag_chan_lo[tbin];
	tch_hi = tch_lo + 3;

	bin_lo = h2d_f->GetYaxis()->FindBin( tch_lo);
	bin_hi = h2d_f->GetYaxis()->FindBin( tch_hi);

	sprintf( string, "Bin %d", tbin);
	if ( strcmp( subt, "full") == 0)
		TH1D *h1 = h2d_f->ProjectionX( string, bin_lo, bin_hi);
	else
		TH1D *h1 = h2d_e->ProjectionX( string, bin_lo, bin_hi);

	// The fit
	TF1 *f1 = new TF1( "f1", "pol1", acc[0], acc[1]+10);
	h1->Fit( "f1", "RLMQ");
	f1->GetParameters( &param[0]);

	ratio = (double) (prompt[1]-prompt[0])/(acc[1]-acc[0])
		* (2*param[0] + (prompt[1] + prompt[0])*param[1])/(2*param[0]
			+ (acc[1] + acc[0])*param[1]);

	max = h1->GetMaximum();

	max *= 1.1;

	h1->SetMaximum( max);

	l1 = new TLine( prompt[0], 0., prompt[0], max);
	l1->SetLineStyle(1);
	l1->SetLineWidth(2);
	l1->Draw();

	l2 = new TLine( prompt[1], 0., prompt[1], max);
	l2->SetLineStyle(1);
	l2->SetLineWidth(2);
	l2->Draw();

	l3 = new TLine( acc[0], 0., acc[0], max);
	l3->SetLineStyle(2);
	l3->SetLineWidth(2);
	l3->Draw();

	l4 = new TLine( acc[1], 0., acc[1], max);
	l4->SetLineStyle(2);
	l4->SetLineWidth(2);
	l4->Draw();

	max1 = max*0.87;
	max2 = max*0.93;

	sprintf( string, "Subtraction Ratio = %6.4f", ratio);
	pt = new TPaveText( 20, max1, 90, max2);
	pt->AddText( string);
	pt->Draw();

	return ratio;
}

void ChanZero( Char_t subt[8])
//void ChanZero( Char_t subt[8], Int_t chan)
{
	Int_t i, Ny;
	Int_t* lch = broken_ladd_ch;

	Ny = h2d_f->GetXaxis()->GetNbins();

	while ( *lch != -1) {
		for ( i = 0; i < Ny; i++) {
			if ( strcmp( subt, "full") == 0)
				h2d_f->SetBinContent( i, *lch, 0);
			else
				h2d_e->SetBinContent( i, *lch, 0);
//			cout << i << " " << h2d->GetBinContent( i, chan, 0);
//			h2d->SetBinContent( i, chan, 0);
//			cout << " " << h2d->GetBinContent( i, chan, 0) << endl;
		}
		lch++;
	}
}

void ChanBin( Char_t subt[8] = "full", Int_t min = 2, Int_t max = 6)
{
	Int_t bin, flag;
	Double_t ratio[7];
	TString filename;

	TCanvas *c1 = new TCanvas ( "c1", "Prompt-Random", 20, 20, 700, 500);

	ChanZero( subt);

	for ( bin = min; bin <= max; bin++) {
		ratio[bin] = TimeFitBin( subt, bin);
		cout << "bin = " <<  bin << "  subtraction ratio = " << ratio[bin]
			<< endl;
//		c1->WaitPrimitive("");
	}

	if ( strcmp( subt, "full") == 0)
		filename = "xs/chan_subt/bin_full_new.out";
	else
		filename = "xs/chan_subt/bin_empty_new.out";
	ofstream outFile( filename);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}
	for ( bin = min; bin <= max; bin++) {
			outFile << bin << "  " << ratio[bin] << endl ;
	}
	outFile.close();
}

void ChanScan( TString tgt = "full", UInt_t min = 0, UInt_t max = 351,
		UInt_t rebin = 1)
{
	UInt_t chan;

	TCanvas *c1 = new TCanvas ( "c1", "Pi0-Tagger Time", 20, 20, 700, 500);

	for ( chan = min; chan <= max; chan++) {
		Scan( tgt, chan, rebin);
		c1->WaitPrimitive( "");
	}
}

void ChanScanData( TString tgt = "full")
{
	UInt_t i;

	name = "xs/chan_subt/tdcpeaks.dat";
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	for ( i = 0; i <= 351; i++)
	{
		Scan( tgt, i, kFALSE, 1); 
		name = Form( "%3d %5.1f %5.1f %4.1f %5.2f", i, tcd[i].energy, par[3],
				par[4], red_chisq);
		outFile << name << endl;
	}
	outFile.close();
}

void Scan( TString tgt = "full", Int_t chan = 300, Bool_t plot = kTRUE,
		UInt_t rebin = 1)
{
	UInt_t bin;
	TString name;

	TH1D* h1;
	TH2D* h2;

	if ( tgt == "full") h2 = h2d_f;
	else if ( tgt == "empty") h2 = h2d_e;

	if ( ( chan >= 0) && ( chan <= 351)) {
		bin = h2->GetXaxis()->FindBin( chan);
		name = Form( "Channel %d", chan);
		h1 = h2->ProjectionY( name, bin, bin);
		name = Form( "E_{#gamma} = %5.1f #pm %3.1f MeV", tcd[chan].energy,
				tcd[chan].denergy);
	}
	else if ( chan == -1) {
		name = "All Channels";
		h1 = h2->ProjectionY( name, 0, 351);
	}
	else {
		cout << "Error: channels must be 0-351 or -1" << endl;
		break;
	}

	h1->SetTitle( name);

	TString opt;
	opt = "NQR+";
	if ( plot == kTRUE) {
		h1->Draw();
		opt = "QR+";
	}

	TF1 *f1 = new TF1( "f1", "pol1", 20, 80);
	TF1 *f2 = new TF1( "f2", "gaus", 95, 105);
	h1->Fit( "f1", "NQR");
	f1->GetParameters( &par[0]);
	h1->Fit( "f2", "NQR+");
	f2->GetParameters( &par[2]);
	TF1 *fall = new TF1( "fall", "pol1(0) + gaus(2)", 20, 180);
	fall->SetParameters( par);
	h1->Fit( "fall", opt);
	Int_t bins = h1->GetXaxis()->FindBin( 180) - h1->GetXaxis()->FindBin( 20);
	red_chisq = fall->GetChisquare()/(bins-5);

//	gStyle->SetOptStat( 0);
	gStyle->SetOptFit( 1);

	cout << chan;
	cout << " " << par[3];
	cout << " " << par[4];
	cout << " " << red_chisq;
	cout << endl;

}

void Eff()
{
	TCanvas *c1 = new TCanvas ( "c1", "Crude Efficiency", 20, 20, 700, 500);
	c1->Divide(1,3);

	c1->cd(1);
	hits->GetXaxis()->SetRangeUser( 228, 258);
	hits->DrawCopy();

	c1->cd(2);
	scal->Draw();
	scal->GetXaxis()->SetRangeUser( 228, 258);

	c1->cd(3);
	hits->Divide( scal);
	hits->Draw();
	hits->GetXaxis()->SetRangeUser( 228, 258);

	TCanvas *c2 = new TCanvas ( "c2", "Crude Efficiency", 800, 20, 700, 500);
	hits->Draw();
	hits->GetXaxis()->SetRangeUser( 228, 258);
}
