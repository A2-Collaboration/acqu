#include "/home/dave/Code/lib/physics.h"
#include "/home/dave/Code/lib/functions.h"

TH2D* hS;

void SubtHist()
{
	TString name;

	file = new TFile("histograms/ARH.root");

	TH2D* hP;
	TH2D* hR;

	name = "THR_NPi0PhotP_v_TChanP";
	hP = (TH2D*) file->Get( name);
	name = "THR_NPi0PhotR_v_TChanR";
	hR = (TH2D*) file->Get( name);

	hS = (TH2D*) hP->Clone( "Subt2D");

	hS->Sumw2();
	hS->Add( hR, -0.0833);
}

void NPi0Phot( UInt_t chan)
{
	UInt_t bin[3];
	Double_t cts[3], dcts[3];
	Double_t ratio[2], dratio[2];
	TString name;

	TH1D* h1;

	hS->GetYaxis()->SetRange( 0, 8);
	hS->GetXaxis()->SetRange( chan+1, chan+1);

	name = Form( "npi0phot%d", chan);
	h1 = (TH1D*) hS->ProjectionY( name);

	bin[0] = h1->GetXaxis()->FindBin( 2);
	bin[1] = h1->GetXaxis()->FindBin( 3);
	bin[2] = h1->GetXaxis()->FindBin( 4);
	cts[0] = h1->GetBinContent( bin[0]);
	cts[1] = h1->GetBinContent( bin[1]);
	cts[2] = h1->GetBinError( bin[2]);
	dcts[0] = h1->GetBinError( bin[0]);
	dcts[1] = h1->GetBinError( bin[1]);
	dcts[2] = h1->GetBinError( bin[2]);

//	cout << bin[2];
//	cout << " " << cts[2];
//	cout << endl;

	ratio[0] = cts[1]/cts[0]*100;
	dratio[0] =
		100*sqrt( Sqr( dcts[1]/cts[0]) + Sqr( cts[1]*dcts[0]/Sqr( cts[0])));
	ratio[1] = cts[2]/cts[0]*100;
	dratio[1] =
		100*sqrt( Sqr( dcts[2]/cts[0]) + Sqr( cts[2]*dcts[0]/Sqr( cts[0])));
	name = Form( "Ratio of counts 3:2 = %5.1f +/- %5.1f\n", ratio[0], dratio[0]);
	cout << name;
	name = Form( "Ratio of counts 4:2 = %5.1f +/- %5.1f\n", ratio[1], dratio[1]);
	cout << name;

	TCanvas *c1 = new TCanvas ( "c1", "Projection", 20, 20, 500, 700);
	h1->Draw( "e");
}
