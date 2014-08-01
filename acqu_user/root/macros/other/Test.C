void Plot4()
{
	TCanvas *canvas = new TCanvas ( "c1", "Plot4", 100, 0, 700, 500);
	canvas->Divide( 2, 2);

	for ( UInt_t tbin = 0; tbin < 4; tbin++) {
		canvas->cd( tbin+1);
		PlotSingle( tbin);
	}
}

void PlotSingle( UInt_t tbin)
{
	Double_t tcut, max;

	TCanvas* canv;
	if ( !( canv = (TCanvas*)(gROOT->FindObject( "c1"))))
		TCanvas* canv = new TCanvas ( "c1", "PlotSingle", 100, 0, 700, 500);

	tcut = tbin*20;

	// Histogram 1.
	proton = new TFile( "proton.root");
	h3p = (TH3D*) proton->Get( "MMiss_v_TGG_v_Theta");
	h3p->GetXaxis()->SetRangeUser( tcut, 180);
	TH1D *h1 = h3p->Project3D( "z");

	// Histogram 2.
	carbon = new TFile( "carbon.root");
	h3c = (TH3D*) carbon->Get( "MMiss_v_TGG_v_Theta");
	h3c->GetXaxis()->SetRangeUser( tcut, 180);
	TH1D *h2 = h3c->Project3D( "z2");

	// Add h1 and h2 and plot it as well.
	TH1D *h3 = (TH1D*)h1->Clone();
	h3->Sumw2();
	h3->Add( h2, 1.0);

	// Set maximum for frame and draw it.
	if ( h1->GetMaximum() > h2->GetMaximum()) max = 1.10*h1->GetMaximum();
	else max = 1.10*h2->GetMaximum();
	TH1F *frame = gPad->DrawFrame( 900, 0, 1000, max);
	frame->Draw();

	// Plot h1 and h2 on the frame.
	h1->Draw( "same");
	h2->Draw( "same");
	h3->Draw( "same");
}
