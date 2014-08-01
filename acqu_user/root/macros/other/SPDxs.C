gROOT->Reset();

typedef struct {
	UInt_t chan;
	UInt_t chan_sp;
	Double_t e_lo;
	Double_t e_hi;
	Double_t eg;
	Double_t deg;
	Double_t cth[21];
	Double_t th[21];
	Double_t xs[21];
	Double_t dxs[21];
	Double_t tot_xs;
	Double_t tot_dxs;
	Double_t eff;
} TBins;

TBins tbin[353];

Double_t eg_as[] = { 145.1, 146.2, 147.4, 148.5, 149.6, 150.8, 151.9, 153.1,
	154.2, 155.3, 156.5, 157.6, 158.8, 159.9, 161.0, 162.2, 163.3, 164.4, 165.6,
	166.7, 167.8, 169.0, 202.0, 210.0 };

void SortChanData()
{
	UInt_t i, j, imax, chan;
	Double_t cth, th, xs, dxs, eff;
	TString name, string;
	Char_t chname[10];

//	name = "xs/pi0/xs_sergey/sergey.txt";
	name = "xs/pi0/xs_sergey/sergey_v3.txt";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}

//	imax = 352;
	imax = 90;

	for ( i = 0; i < imax; i++) {

		tbin[i].chan = 308 - i;

		// Tagger channel
		name.ReadLine( inFile);
		cout << name << endl;
		sscanf( name.Data(), "%*s%*s%s%*s%lf", chname, &eff);
		cout << "chname = " << chname << endl;
		UInt_t n = sizeof( chname);
		chname[n] = '\0';
		chan = atoi( chname);

		cout << "chan = " << chan << endl;
		cout << "eff = " << eff << endl;

		tbin[i].chan_sp = chan;
		tbin[i].eff = eff;

		name.ReadLine( inFile);
		// Bin lower limit
		string = name;
		string.Remove( 0, string.First('=')+1);
		string.Remove( string.First('-'), string.Length());
		tbin[i].e_lo = string.Atof();
		// Bin upper limit
		string = name;
		string.Remove( 0, string.First('-')+1);
		string.Remove( string.First(' ')+1, string.Length());
		tbin[i].e_hi = string.Atof();

		// Junk
		name.ReadLine( inFile);

		// Cos theta, theta, xs, dxs
		for ( j = 0; j < 20; j++) {
			name.ReadLine( inFile);
			sscanf( name.Data(), "%lf%lf%lf%lf", &cth, &th, &xs, &dxs);
			tbin[i].cth[j] = cth;
			tbin[i].th[j] = th;
			tbin[i].xs[j] = xs;
			tbin[i].dxs[j] = dxs;
		}

		name.ReadLine( inFile);
		// Total cross section
		string = name;
		string.Remove( 0, string.First('=')+1);
		string.Remove( string.First('+'), string.Length());
		tbin[i].tot_xs = string.Atof();
		// Error in total cross section
		string = name;
		string.Remove( 0, string.First('-')+1);
		string.Remove( string.Length()-4, string.Length());
		tbin[i].tot_dxs = string.Atof();

		tbin[i].eg = (tbin[i].e_lo + tbin[i].e_hi)/2;
		tbin[i].deg = (tbin[i].e_hi - tbin[i].e_lo)/2;

		// as a function of theta_cm
		name = Form( "xs/pi0/xs_sergey/dxs_chan%d.dat", tbin[i].chan);
		ofstream outFile( name);
		if ( !outFile.is_open()) {
			cout << "Error opening file ";
			cout << name;
			cout << endl;
			break;
		}

		name = Form( "#\n# Channel %d\n#\n# Egamma = %5.1f +/- %4.1f MeV\n",
				tbin[i].chan, tbin[i].eg, tbin[i].deg);
		outFile << name;
		name = Form( "# TotXS = %7.3f +/- %6.3f ub\n#\n", tbin[i].tot_xs,
				tbin[i].tot_dxs);
		outFile << name;
		name = "# theta   xs     dxs\n#\n";
		outFile << name;
		for ( j = 0; j < 20; j++) {
			name = Form( "%6.2f  %6.4f  %6.4f\n", tbin[i].th[j], tbin[i].xs[j],
					tbin[i].dxs[j]);
			outFile << name;
		}
		outFile.close();

		// as a function of cos(theta_cm)
		name = Form( "xs/pi0/xs_sergey/dxs_cth_chan%d.dat", tbin[i].chan);
		ofstream outFile( name);
		if ( !outFile.is_open()) {
			cout << "Error opening file ";
			cout << name;
			cout << endl;
			break;
		}

		name = Form( "#\n# Channel %d\n#\n# Egamma = %5.1f +/- %4.1f MeV\n",
				tbin[i].chan, tbin[i].eg, tbin[i].deg);
		outFile << name;
		name = Form( "# TotXS = %7.3f +/- %6.3f ub\n#\n", tbin[i].tot_xs,
				tbin[i].tot_dxs);
		outFile << name;
		name = "# theta   xs     dxs\n#\n";
		outFile << name;
		for ( j = 0; j < 20; j++) {
			name = Form( "%5.2f  %6.4f  %6.4f\n", tbin[i].cth[j], tbin[i].xs[j],
					tbin[i].dxs[j]);
			outFile << name;
		}
		outFile.close();
	}
	inFile.close();

	name = "xs/pi0/xs_sergey/tot_xs.dat";
	ofstream outFile( name);
	if ( !outFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}

	for ( i = 0; i < 28; i++) {
		name = Form( "%5.1f %7.3f %6.3f\n", tbin[i].eg, tbin[i].tot_xs,
				tbin[i].tot_dxs);
		outFile << name;
	}
	outFile.close();

}

//
// NOTE: you MUST run SortChanData() first!!!
//
void PlotCompDXS( UInt_t bin, Bool_t save = kFALSE)
{
	UInt_t i, count, ct;
	Double_t diff, ecomp;
	Double_t th[21], dth[21], xs[21], dxs[21];
	Double_t theta, xxs, dxxs;
	TString name;

	count = 0;
	for ( i = 0; i < 24; i++) {
		diff = fabs( tbin[bin].eg - eg_as[i]);
		if ( diff <= 0.5 ) {
			ecomp = eg_as[i];
			count++;
		}
	}
	
	if ( count != 1) {
		cout << "No TAPS match for bin " << bin << " corresponding to ";
		cout << tbin[bin].eg << " MeV" << endl;
	}

	c1 = new TCanvas( "c1", "Differential Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

	for ( i = 0; i < 20; i++) {
		th[i] = tbin[bin].th[i];
		dth[i] = 0;
		xs[i] = tbin[bin].xs[i];
		dxs[i] = tbin[bin].dxs[i];
	}

	// Plot Results
	gr = new TGraphErrors( 20, th, xs, dth, dxs);
	name = Form( "Preliminary #gammap#rightarrowp#pi^{0} "
			"Differential Cross Section for E_{#gamma} = %6.1f #pm %4.1f MeV",
			tbin[bin].eg, tbin[bin].deg);
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.1);
	gr->GetYaxis()->SetTitleOffset( 1.0);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "#theta^{cm}_{#pi^{0}} (deg)");
	gr->GetYaxis()->SetTitle( "d#sigma/d#Omega (#mub/sr)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetMinimum( 0);
	gr->GetXaxis()->SetRangeUser( 0, 180);
	gr->Draw( "AP");

	if ( count == 1) {
		// Read in Schmidt data
		name = Form( "xs/pi0/previous/dxs_%5.1fMeV.dat", ecomp);
		ifstream inFile( name);
		if ( !inFile.is_open()) {
			cout << "Error opening file ";
			cout << name;
			cout << endl;
			break;
		}
		i = 0;
		while( !inFile.eof()) {
			inFile >> theta >> xxs >> dxxs;
			th[i] = theta;
			dth[i] = 0;
			xs[i] = xxs;
			dxs[i++] = dxxs;
		}
		inFile.close();
		ct = i-1;

		// Plot Results
		gr1 = new TGraphErrors( ct, th, xs, dth, dxs);
		gr1->SetMarkerColor( 2);
		gr1->SetMarkerSize( 1.2);
		gr1->SetLineWidth( 2);
		gr1->SetLineColor( 2);
		gr1->SetMarkerStyle( 20);
		gr1->Draw( "Psame");

		pt = new TLegend( 0.6, 0.15, 0.8, 0.30);
		pt->SetTextSize( 0.04);
		pt->SetFillColor( 0);
		pt->SetBorderSize( 0);
		pt->AddEntry( gr, "This Work", "p");
		pt->AddEntry( gr1, "TAPS", "p");
		pt->Draw();
	}

	name = Form( "plots/Pi0/dxs_comp_%d", bin);
	name.Append( ".pdf");
//	c1->Print( name);
}

//
// NOTE: you MUST run SortChanData() first!!!
//
void PlotCompTotXS( Bool_t save = kFALSE)
{
	UInt_t i, count, ct, imax;
	Double_t eg[353], deg[353], xs[353], dxs[353];
	Double_t egg, xxs, dxxs;
	TString name;

	c1 = new TCanvas( "c1", "Total Cross Sections", 200, 10, 700, 500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor( 21);
	c1->GetFrame()->SetBorderSize( 12);

//	imax = 28;
	imax = 90;

	for ( i = 0; i < imax; i++) {
		eg[i] = tbin[i].eg;
		deg[i] = 0;
		xs[i] = tbin[i].tot_xs;
		dxs[i] = tbin[i].tot_dxs;
	}

	// Plot Results
	gr = new TGraphErrors( imax, eg, xs, deg, dxs);
	name = "Preliminary #gammap#rightarrowp#pi^{0} Total Cross Sections";
	gr->SetTitle( name);
	gr->SetMarkerColor( 4);
	gr->SetMarkerSize( 1.2);
	gr->SetMarkerStyle( 21);
	gr->SetLineWidth( 2);
	gr->SetLineColor( 4);
	gr->GetXaxis()->SetTitleOffset( 1.2);
	gr->GetYaxis()->SetTitleOffset( 0.8);
	gr->GetYaxis()->SetTitleSize( 0.05);
	gr->GetXaxis()->SetTitle( "E_{#gamma} (MeV)");
	gr->GetYaxis()->SetTitle( "#sigma (#mub)");
	gr->GetXaxis()->SetLabelSize( 0.03);
	gr->GetYaxis()->SetLabelSize( 0.03);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->SetMinimum( 0);
	gr->SetMaximum( 350);
	gr->GetXaxis()->SetRangeUser( 140, 400);
	gr->Draw( "AP");

	// Read in TAPS data
	name = "xs/pi0/previous/tot_xs.dat";
	ifstream inFile( name);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << name;
		cout << endl;
		break;
	}
	i = 0;
	while( !inFile.eof()) {
		inFile >> egg >> xxs >> dxxs;
		eg[i] = egg;
		deg[i] = 0;
		xs[i] = xxs;
		dxs[i++] = dxxs;
	}
	inFile.close();
	ct = i-1;

	// Plot Results
	gr1 = new TGraphErrors( ct, eg, xs, deg, dxs);
	gr1->SetMarkerColor( 2);
	gr1->SetMarkerSize( 1.2);
	gr1->SetLineWidth( 2);
	gr1->SetLineColor( 2);
	gr1->SetMarkerStyle( 20);
	gr1->Draw( "Psame");

	pt = new TLegend( 0.2, 0.6, 0.4, 0.75);
	pt->SetTextSize( 0.04);
	pt->SetFillColor( 0);
	pt->SetBorderSize( 0);
	pt->AddEntry( gr, "This Work", "p");
	pt->AddEntry( gr1, "TAPS", "p");
	pt->Draw();

	name = "plots/Pi0/tot_xs_comp";
	name.Append( ".pdf");
	if ( save == kTRUE) c1->Print( name);
}
