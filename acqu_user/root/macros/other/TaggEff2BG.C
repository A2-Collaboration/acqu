//
// TaggEff2BG.C
//
// Extract tagging efficiency values and statistical uncertainties
// using background runs.
//
// 2008.10.02		J.R.M. Annand			First Version
// 2009.05.31		D.L. Hornidge			A few modifications
// 2010.01.07		D.L. Hornidge			Added TaggEffThresh macro
//

TGraphErrors* TaggEffPlot;
TH1F* histEelec;             // electron energy
TH1F* histHits;              // hits histogram
TH1D* histSc;                // scaler histogram
TH1D* histScB;               // background scaler histogram
TH1D* histScB1;              // background scaler histogram
TH1D* histS;                 // trigger scalers histogram
TH1D* histSB;                // background trigger scalers histogram
TH1D* histSB1;               // background trigger scalers histogram

Double_t* Eelec;             // electron energy
Double_t* Hits;              // hits array
Double_t* Scalers;           // accumulated scalers
Double_t* ScalersB;          // accumulated scalers background run
Double_t* TaggEff;           // tagging efficiency
Double_t* dTaggEff;          // uncertainty
Double_t* TaggEffB;          // background corrected tagging efficiency
Double_t* dTaggEffB;         // uncertainty

Double_t meanEff, meanEffB, dmeanEff, dmeanEffB;

void TaggEff2BG( Int_t rfore, Int_t rback1=0, Int_t rback2=0,
		Int_t chan_lo=0, Int_t chan_hi=351, Int_t jclk=31, Int_t jinh=32,
		Int_t liveCorr=0)
{

	Int_t nchan = 352;
	TString ffore, fback1, fback2;
	TString string;

	ffore = Form( "histograms/TaggEff/TaggEff_%d.root", rfore);
	if ( rback1 != 0)
		fback1 = Form( "histograms/TaggEff/TaggEff_%d.root", rback1);
	if ( rback2 != 0)
		fback2 = Form( "histograms/TaggEff/TaggEff_%d.root", rback2);

	// calculate tagging efficiency
	Eelec = new Double_t[nchan];
	Hits = new Double_t[nchan];
	Scalers = new Double_t[nchan];
	ScalersB = new Double_t[nchan];
	TaggEff = new Double_t[nchan];
	dTaggEff = new Double_t[nchan];
	TaggEffB = new Double_t[nchan];
	dTaggEffB = new Double_t[nchan];

	TFile* fore = new TFile( ffore);
	TFile* back1 = NULL;
	TFile* back2 = NULL;
	if ( !fore )
	{
		printf("File %s not found\n", ffore.Data());
		return;
	}
	if ( rback1 != 0 ) back1 = new TFile( fback1);
	if ( !back1 ) cout << "Assume no 1st background correction run" << endl;
	if ( rback2 != 0 ) back2 = new TFile( fback2);
	if ( !back2 ) cout << "Assume no 2nd background correction run" << endl;

	TString logfile = Form( "xs/tag_eff/tageff_%d.dat", rfore);
	cout << "Text output saved to file: " << logfile << endl;

	ofstream TaggEffLog( logfile);
	if ( !TaggEffLog.is_open()) 
	{
	  cout << "Unable to open " << logfile << endl;
	  return;
	}

	histEelec =  (TH1F*)fore->Get("FPD_Eelect");
	histHits = (TH1F*)fore->Get("FPD_Hits");
	//  histSc = (TH1D*)fore->Get("FPD_ScalerAcc");
	histSc = (TH1D*)fore->Get("SumScalers12to363");
	histS = (TH1D*)fore->Get("SumScalers364to405");
	if( back1 ){
	  //    histScB = (TH1D*)back1->Get("FPD_ScalerAcc");
	  histScB = (TH1D*)back1->Get("SumScalers12to363");
	  histSB = (TH1D*)back1->Get("SumScalers364to405");
	}
	if( back2 ){
	  //    histScB1 = (TH1D*)back2->Get("FPD_ScalerAcc");
	  histScB1 = (TH1D*)back2->Get("SumScalers12to363");
	  histSB1 = (TH1D*)back2->Get("SumScalers364to405");
	}

	if( !histEelec ) printf("Electron energy hist not found\n");
	if( !histHits ) printf("Hits hist not found\n");
	if( !histSc ) printf("Scalers hist not found\n");
	if( !histS ) printf("Trig-Scalers hist not found\n");
	if( back1 ){
	  if( !histScB ) printf("1st Background Scalers hist not found\n");
	  if( !histSB ) printf("1st Background Trig-Scalers hist not found\n");
	}
	if( back2 ){
	  if( !histScB1 ) printf("2nd Background Scalers hist not found\n");
	  if( !histSB1 ) printf("2nd Background Trig-Scalers hist not found\n");
	}

	Double_t clk = (Double_t)histS->GetBinContent(jclk);
	Double_t inh = (Double_t)histS->GetBinContent(jinh);
	Double_t live = inh/clk;
	if( back1 ){
	  Double_t clkB = (Double_t)histSB->GetBinContent(jclk);
	  Double_t inhB = (Double_t)histSB->GetBinContent(jinh);
	  Double_t liveB = inhB/clkB;
	}
	else {
	  Double_t clkB = 0;
	  Double_t inhB = 0;
	  Double_t liveB = 0;
	}
	if( back2 ){
	  Double_t clkB = clkB + (Double_t)histSB1->GetBinContent(jclk);
	  Double_t inhB = inhB + (Double_t)histSB1->GetBinContent(jinh);
	  Double_t liveB = inhB/clkB;
	}  
	//  printf("%d %d %5.3f %d %d %5.3f\n",clk,inh,live,clkB,inhB,liveB);

	Double_t bscale;
	if( inhB ) bscale = inh/inhB;
	else bscale = 1.0;
	Double_t liveC;
	if( liveCorr ) liveC = live;
	else liveC = 1.0;
	Double_t ScB, dScB;

	string = " Tagging efficiency calculation\n J.R.M.Annand 2008.10.02";
	TaggEffLog << string << endl;

	string = Form( "Clock = %10g   Inhib.Clock = %10g   Live = %8g\n"
		"Backgr.clock = %10g    Backgr.Inhib.clock = %10g  Live = %8g\n"
		"Background scale factor = %8g",
		clk, inh, live, clkB, inhB, liveB, bscale);

	TaggEffLog << string << endl << endl;

	string = "  Chan  Hits  Scaler Background     TaggEff        "
		"TaggEff(back-sub)";
	TaggEffLog << string << endl;

	for ( UInt_t i = chan_lo; i <= chan_hi; i++)
	{

	  // Eelec[i] = (Double_t)histEelec->GetBinContent(i);
	  Eelec[i] = (Double_t)i;
	  Hits[i]  = (Double_t)histHits->GetBinContent(i);
	  Scalers[i] = (Double_t)histSc->GetBinContent(i);

	  if ( back1 ) ScalersB[i] = (Double_t)histScB->GetBinContent(i);
	  else ScalersB[i] = 0;

	  if ( back2 )
		  ScalersB[i] = ScalersB[i] + (Double_t)histScB1->GetBinContent(i);
	  if ( (Scalers[i] == 0) || (Hits[i] == 0) ) {
	    TaggEff[i] = TaggEffB[i] = dTaggEff[i] = dTaggEffB[i] = 0.0;
	  }
	  else {
	    ScB = Scalers[i] - bscale*ScalersB[i];
	    dScB = TMath::Sqrt(ScalersB[i]) * bscale;
	    dScB = TMath::Sqrt( Scalers[i] + dScB*dScB )/ScB;
	    
	    TaggEff[i] = liveC * Hits[i]/Scalers[i];
	    TaggEffB[i] = liveC * Hits[i]/ScB;
	    dTaggEff[i] = TaggEff[i] * TMath::Sqrt( 1./Hits[i] + 1./Scalers[i] );
	    dTaggEffB[i] = TaggEffB[i] * TMath::Sqrt( 1./Hits[i] + dScB*dScB );
	  }

		string = Form( "%4.0f %5.0f %8.0f %8.0f   %7.4f +- %6.4f  %7.4f +- %6.4f",
			Eelec[i], Hits[i], Scalers[i], ScalersB[i], TaggEff[i], dTaggEff[i],
			TaggEffB[i], dTaggEffB[i]);
		TaggEffLog << string << endl;
	}
	TaggEffDraw( chan_lo, chan_hi, rfore);

	string = Form( "\nMean tagging efficiency %7.4f +- %7.4f\n\
	                      %7.4f +- %7.4f without background subtraction\n",
		meanEffB, dmeanEffB, meanEff, dmeanEff);
	TaggEffLog << string << endl;

	TaggEffLog.close();

	fore->Close();

	if ( rback1 != 0) back1->Close();
	if ( rback2 != 0) back2->Close();
}

void TaggEffDraw( Int_t chan_lo, Int_t chan_hi, Int_t run)
{
	TString title, name;

	// refills the TeffData into gr_data avoiding "bad" channels and
	// plots the graph.

	TCanvas* canv;

	if ( !( canv = (TCanvas*)(gROOT->FindObject("Teff")) ) )
	{
	  TCanvas* canv = new TCanvas( "Teff", "Tagging Efficiency", 240, 0, 700,
			  700);
	  canv->SetFillStyle( 4000);
	  canv->Divide( 1, 2, 0.01, 0.01);
	  canv->cd(1);
	}
	else canv->cd(1);

	UInt_t j = 0;
	for ( UInt_t i = chan_lo; i <= chan_hi; i++)
	{ 
	  if ( TaggEff[i] > 0 ) j++;
	}

	Double_t* ee = new Double_t[j];
	Double_t* tf = new Double_t[j];
	Double_t* dtf = new Double_t[j];
	j = 0;
	meanEff = dmeanEff = 0.0;
	for ( UInt_t i = chan_lo; i <= chan_hi; i++)
	{
	  if( TaggEff[i] > 0 ){
	    ee[j] = Eelec[i];
	    tf[j] = TaggEff[i];
	    meanEff = meanEff + tf[j];
	    dtf[j] = dTaggEff[i];
	    dmeanEff = dmeanEff + dtf[j]*dtf[j];
	    j++;
	  }
	}
	meanEff = meanEff/j;
	dmeanEff = TMath::Sqrt( dmeanEff)/j;

	TaggEffPlot = new TGraphErrors( j, ee, tf, NULL, dtf);
	title = Form( "Tagging-Efficiency of Run %d (Preliminary)", run);
	TaggEffPlot->SetTitle( title);
	TaggEffPlot->SetMaximum( 0.6);
	TaggEffPlot->SetMinimum( 0.0);
	TaggEffPlot->GetXaxis()->SetTitle( "Tagger Channel");
	TaggEffPlot->GetYaxis()->SetTitle( "#epsilon_{tagg}");
	TaggEffPlot->GetXaxis()->CenterTitle();
	TaggEffPlot->GetYaxis()->CenterTitle();
	TaggEffPlot->SetMarkerColor( kBlue);
	TaggEffPlot->SetMarkerStyle( 21);
	TaggEffPlot->SetMarkerSize( 0.7);
	TaggEffPlot->Draw( "AP");

	canv->cd(2);
	j = 0;
	meanEffB = 0.0;
	for ( UInt_t i = chan_lo; i <= chan_hi; i++){
	  if ( TaggEff[i] > 0 ){
	    ee[j] = Eelec[i];
	    tf[j] = TaggEffB[i];
	    meanEffB = meanEffB + tf[j];
	    dtf[j] = dTaggEffB[i];
	    dmeanEffB = dmeanEffB + dtf[j]*dtf[j];
	    j++;
	  }
	}
	meanEffB = meanEffB/j;
	dmeanEffB = TMath::Sqrt( dmeanEffB)/j;

	TaggEffPlotB = new TGraphErrors( j, ee, tf, NULL, dtf);
	title = Form( "Background-Corrected Tagging-Efficiency of Run %d", run);
	title.Append( " (Preliminary)");
	TaggEffPlotB->SetTitle( title);
	TaggEffPlotB->SetMaximum( 0.6);
	TaggEffPlotB->SetMinimum( 0.0);
	TaggEffPlotB->GetXaxis()->SetTitle( "Tagger Channel");
	TaggEffPlotB->GetYaxis()->SetTitle( "#epsilon_{tagg}");
	TaggEffPlotB->GetXaxis()->CenterTitle();
	TaggEffPlotB->GetYaxis()->CenterTitle();
	TaggEffPlotB->SetMarkerColor( kRed);
	TaggEffPlotB->SetMarkerStyle( 20);
	TaggEffPlotB->SetMarkerSize( 0.7);
	TaggEffPlotB->Draw( "AP");

	name = Form( "plots/TaggEff_%d", run);
	name.Append( ".pdf");
	canv->Print( name);
}

void TaggEffThresh( Int_t n = 0)
{
	Int_t i;
	Int_t m[] = { 20842, 20887, 20911, 20951, 20978, 21057};
	Int_t rfore[3], rback1, rback2;
	Int_t chan_lo, chan_hi, jclk, jinh, liveCorr;

	chan_lo = 0;
	chan_hi = 351;
	jclk = 31;
	jinh = 32;
	liveCorr = 0;

	rback1 = m[n]-2;
	rfore[0] = m[n]-1;
	rfore[1] = m[n];
	rfore[2] = m[n]+1;
	rback2 = m[n]+2;

	// Correction for fucked-up order for 20911 series.
	if ( n == 2) {
		rback2 = m[n] + 1;
		rfore[2] = m[n] + 2;
	}

	cout << "Analysing Tag Eff series " << m[n] << endl;
	for ( i = 0; i <= 2; i++)
	{
		TaggEff2BG( rfore[i], rback1, rback2, chan_lo, chan_hi, jclk, jinh,
				liveCorr);
	}
}
