gROOT->Reset();

void TotPhot()
{
	UInt_t i, nphotons;
	Int_t binphot, photsbin, totphot;

	TH1D* nphot;

	nphot = (TH1D*) gROOT->FindObject( "COMP_Nphoton");

	totphot = 0;
	for ( i = 1; i <= 10; i++)
	{
		nphotons = i - 1;
		binphot = nphot->GetBinContent( i);
		photsbin = nphotons*binphot;
		totphot += photsbin;
		cout << nphotons;
		cout << " " << binphot;
		cout << " " << photsbin;
		cout << " " << totphot;
		cout << endl;
	}
	cout << "tot photons = " << totphot << endl;
}
