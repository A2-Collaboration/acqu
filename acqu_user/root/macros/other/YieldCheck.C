gROOT->Reset();

#include "includes/physics.h"
#include "includes/functions.h"

TFile full( "histograms/ARH_full.root");

void GetYield()
{
	Int_t xmin, xmax, ymin, ymax, zmin, zmax;
	Double_t count;
	Char_t name[64];


	TH3D *h = (TH3D*)full.Get( "THits_TGG_ThetaCM_P");
	xmin = h->GetXaxis()->FindBin( 0);
	xmax = h->GetXaxis()->FindBin( 180);
	ymin = h->GetYaxis()->FindBin( 0);
	ymax = h->GetYaxis()->FindBin( 180);
	zmin = h->GetZaxis()->FindBin( 0);
	zmax = h->GetZaxis()->FindBin( 352);
	count = h->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;

	count = 0;
	TH3D *h1 = (TH3D*)full.Get( "THits_TGG_Phi_1perp_P");
	xmin = h1->GetXaxis()->FindBin( -180);
	xmax = h1->GetXaxis()->FindBin( 180);
	ymin = h1->GetYaxis()->FindBin( 0);
	ymax = h1->GetYaxis()->FindBin( 180);
	zmin = h1->GetZaxis()->FindBin( 0);
	zmax = h1->GetZaxis()->FindBin( 352);
	count += h1->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h2 = (TH3D*)full.Get( "THits_TGG_Phi_2perp_P");
	count += h2->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h3 = (TH3D*)full.Get( "THits_TGG_Phi_3perp_P");
	count += h3->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h4 = (TH3D*)full.Get( "THits_TGG_Phi_4perp_P");
	count += h4->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h5 = (TH3D*)full.Get( "THits_TGG_Phi_5perp_P");
	count += h5->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h6 = (TH3D*)full.Get( "THits_TGG_Phi_6perp_P");
	count += h6->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h7 = (TH3D*)full.Get( "THits_TGG_Phi_7perp_P");
	count += h7->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h8 = (TH3D*)full.Get( "THits_TGG_Phi_8perp_P");
	count += h8->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h9 = (TH3D*)full.Get( "THits_TGG_Phi_9perp_P");
	count += h9->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;

	TH3D *h1 = (TH3D*)full.Get( "THits_TGG_Phi_1para_P");
	count += h1->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h2 = (TH3D*)full.Get( "THits_TGG_Phi_2para_P");
	count += h2->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h3 = (TH3D*)full.Get( "THits_TGG_Phi_3para_P");
	count += h3->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h4 = (TH3D*)full.Get( "THits_TGG_Phi_4para_P");
	count += h4->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h5 = (TH3D*)full.Get( "THits_TGG_Phi_5para_P");
	count += h5->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h6 = (TH3D*)full.Get( "THits_TGG_Phi_6para_P");
	count += h6->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h7 = (TH3D*)full.Get( "THits_TGG_Phi_7para_P");
	count += h7->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h8 = (TH3D*)full.Get( "THits_TGG_Phi_8para_P");
	count += h8->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;
	TH3D *h9 = (TH3D*)full.Get( "THits_TGG_Phi_9para_P");
	count += h9->Integral(xmin, xmax, ymin, ymax, zmin, zmax);
	cout << count << endl;

}
