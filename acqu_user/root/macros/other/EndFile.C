// Stuff to do at the end of each file.

void EndData( TString exp, TString tgt, TString name)
{
	TString path;

	cout << "End-of-File macro executing";
	cout << endl;

	path = "histograms/" + exp + "/" + tgt + "/";

	path.Append( name);
	path.Append( ".root");
	TFile f( path, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	path.Prepend( "All histograms saved to \"");
	path.Append( "\"");
	cout << path;
	cout << endl;
}

void EndFile( TString name)
{
	cout << "End-of-File macro executing" << endl;

	name.Prepend( "histograms/");
	name.Append( ".root");
	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to ");
	cout << name << endl;
}

void EndFileCheck( TString subt = "full")
{
	TString name;

	cout << "End-of-File macro executing" << endl;

	gROOT->ProcessLine( ".L Check.C");

//	YieldCheck( subt);
//	PolCheck( subt);

/*
	name = "histograms/ARH.root";
	TFile f( name, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	name.Prepend( "All histograms saved to ");
	cout << name << endl;
*/

	TA2ThreshPi0* tpi0 = (TA2ThreshPi0*)gAN->GetPhysics();
	tpi0->ZeroAll();
	gAN->ZeroSumScalers();

	name = Form( "Zeroing Sum Scalers and Physics histograms...");
	cout << name << endl;

}

void EndFileMC(Char_t tgt[4], Char_t cond[32], Int_t eg)
{
	Char_t dir[128];

	sprintf( dir, "histograms/MonteCarlo/MCH_%s_%s_%d.root", tgt, cond, eg);

	printf( "End-of-File macro executing\n");

	TFile f( dir, "recreate");
	gROOT->GetList()->Write();
	f.Close();

	printf( "All histograms saved to %s\n", dir);

	gROOT->ProcessLine( ".L ThetaEff.C");

	ThetaEffBatch( tgt, cond, eg);

	printf( "Theta efficiencies calculated for %s %s %d\n", tgt, cond, eg);
}

void TDCDrift( )
{
	Int_t i, run, min, max, peak;
	Char_t file[128], name[128];

//	TA2PhotoPhysics* ph = (TA2PhotoPhysics*)gAN->GetPhysics();
//	run = ph->GetNrun();

//	sprintf( file, "peak_pos%d.out", run);
	sprintf( file, "test.out");

	min = 21;
	max = 283;

	ofstream outFile( file);

	for ( i = min; i <= max; i++) {

		sprintf( name, "THR_Pi0TaggTime%d", i);
		TH1D *h = (TH1D*)gROOT->FindObject( name);
		peak=h->GetMaximumBin();

		outFile << i << "  " << peak << endl;
	}

	outFile.close();
}
