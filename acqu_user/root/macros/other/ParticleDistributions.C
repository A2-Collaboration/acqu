#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"

void Compare(Char_t* filename1 = "TA2Pi0ComptonMC_evgen_Compton.root", Char_t* filename2 = "TA2Pi0ComptonMC_acquMC_Compton.root"){

	canvas1 = new TCanvas("canvas1","PhotonDistributions");
	canvas1->Divide(2,4);

	PhotonDistributions1(filename1);
	PhotonDistributions2(filename2);

	canvas2 = new TCanvas("canvas2","Pi0Distributions");
	canvas2->Divide(2,4);

	Pi0Distributions1(filename1);
	Pi0Distributions2(filename2);

	canvas3 = new TCanvas("canvas3","ProtonDistributions");
	canvas3->Divide(2,4);

	ProtonDistributions1(filename1);
	ProtonDistributions2(filename2);

}

void PhotonDistributions1(Char_t* name = "TA2Pi0Compton2012.root") {

	Char_t* filename = Form("/work0/cristina/%s",name);
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NPhoton;
	Double_t PhotonTheta[100], PhotonPhi[100], PhotonEnergy[100];

	tree1->SetBranchAddress("NPhoton",		&NPhoton);
	tree1->SetBranchAddress("PhotonPhi",		&PhotonPhi);
	tree1->SetBranchAddress("PhotonTheta",		&PhotonTheta);
	tree1->SetBranchAddress("PhotonEnergy",		&PhotonEnergy);

	TH1D *a0 = new TH1D("NPhoton",		"NPhoton",	5,	   0,	  5);
	TH1D *a1 = new TH1D("PhotonPhi",	"PhotonPhi",	360,	-180,	180);
	TH1D *a2 = new TH1D("PhotonTheta",	"PhotonTheta",	180,	   0,	180);
	TH1D *a3 = new TH1D("PhotonEnergy",	"PhotonEnergy",	400,	   0,	400);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		a0->Fill(NPhoton);
		
		for (Int_t a = 0; a < NPhoton; a++){

			a1->Fill(PhotonPhi[a]);
			a2->Fill(PhotonTheta[a]);
			a3->Fill(PhotonEnergy[a]);
		}	
	}

	canvas1->cd(1); a0->Draw();
	canvas1->cd(3); a1->Draw();
	canvas1->cd(5); a2->Draw();
	canvas1->cd(7); a3->Draw();

}

void PhotonDistributions2(Char_t* name = "TA2Pi0Compton2012.root") {

	Char_t* filename = Form("/work0/cristina/%s",name);
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NPhoton;
	Double_t PhotonTheta[100], PhotonPhi[100], PhotonEnergy[100];

	tree1->SetBranchAddress("NPhoton",		&NPhoton);
	tree1->SetBranchAddress("PhotonPhi",		&PhotonPhi);
	tree1->SetBranchAddress("PhotonTheta",		&PhotonTheta);
	tree1->SetBranchAddress("PhotonEnergy",		&PhotonEnergy);

	TH1D *a0 = new TH1D("NPhoton",		"NPhoton",	5,	   0,	  5);
	TH1D *a1 = new TH1D("PhotonPhi",	"PhotonPhi",	360,	-180,	180);
	TH1D *a2 = new TH1D("PhotonTheta",	"PhotonTheta",	180,	   0,	180);
	TH1D *a3 = new TH1D("PhotonEnergy",	"PhotonEnergy",	400,	   0,	400);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		a0->Fill(NPhoton);
		
		for (Int_t a = 0; a < NPhoton; a++){

			a1->Fill(PhotonPhi[a]);
			a2->Fill(PhotonTheta[a]);
			a3->Fill(PhotonEnergy[a]);
		}	
	}

	canvas1->cd(2); a0->Draw();
	canvas1->cd(4); a1->Draw();
	canvas1->cd(6); a2->Draw();
	canvas1->cd(8); a3->Draw();

}

void Pi0Distributions1(Char_t* name = "TA2Pi0Compton2012.root") {

	Char_t* filename = Form("/work0/cristina/%s",name);
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NPi0;
	Double_t Pi0Theta[100], Pi0Phi[100], Pi0Energy[100], Pi0Time[100];

	tree1->SetBranchAddress("NPi0",		&NPi0);
	tree1->SetBranchAddress("Pi0Phi",	&Pi0Phi);
	tree1->SetBranchAddress("Pi0Theta",	&Pi0Theta);
	tree1->SetBranchAddress("Pi0Energy",	&Pi0Energy);

	TH1D *b0 = new TH1D("NPi0",		"NPi0",		5,	   0,	  5);
	TH1D *b1 = new TH1D("Pi0Phi",		"Pi0Phi",	360,	-180,	180);
	TH1D *b2 = new TH1D("Pi0Theta",		"Pi0Theta",	180,	   0,	180);
	TH1D *b3 = new TH1D("Pi0Energy",	"Pi0Energy",	400,	   0,	400);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		b0->Fill(NPi0);

		for (Int_t a = 0; a < NPi0; a++){

			b1->Fill(Pi0Phi[a]);
			b2->Fill(Pi0Theta[a]);
			b3->Fill(Pi0Energy[a]);
		}	
	}

	canvas2->cd(1); b0->Draw();
	canvas2->cd(3); b1->Draw();
	canvas2->cd(5); b2->Draw();
	canvas2->cd(7); b3->Draw();

}

void Pi0Distributions2(Char_t* name = "TA2Pi0Compton2012.root") {

	Char_t* filename = Form("/work0/cristina/%s",name);
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NPi0;
	Double_t Pi0Theta[100], Pi0Phi[100], Pi0Energy[100], Pi0Time[100];

	tree1->SetBranchAddress("NPi0",		&NPi0);
	tree1->SetBranchAddress("Pi0Phi",	&Pi0Phi);
	tree1->SetBranchAddress("Pi0Theta",	&Pi0Theta);
	tree1->SetBranchAddress("Pi0Energy",	&Pi0Energy);

	TH1D *b0 = new TH1D("NPi0",		"NPi0",		5,	   0,	  5);
	TH1D *b1 = new TH1D("Pi0Phi",		"Pi0Phi",	360,	-180,	180);
	TH1D *b2 = new TH1D("Pi0Theta",		"Pi0Theta",	180,	   0,	180);
	TH1D *b3 = new TH1D("Pi0Energy",	"Pi0Energy",	400,	   0,	400);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		b0->Fill(NPi0);

		for (Int_t a = 0; a < NPi0; a++){

			b1->Fill(Pi0Phi[a]);
			b2->Fill(Pi0Theta[a]);
			b3->Fill(Pi0Energy[a]);
		}	
	}

	canvas2->cd(2); b0->Draw();
	canvas2->cd(4); b1->Draw();
	canvas2->cd(6); b2->Draw();
	canvas2->cd(8); b3->Draw();

}

void ProtonDistributions1(Char_t* name = "TA2Pi0Compton2012.root") {

	Char_t* filename = Form("/work0/cristina/%s",name);
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NProton;
	Double_t ProtonTheta[100], ProtonPhi[100], ProtonEnergy[100];

	tree1->SetBranchAddress("NProton",		&NProton);
	tree1->SetBranchAddress("ProtonPhi",		&ProtonPhi);
	tree1->SetBranchAddress("ProtonTheta",		&ProtonTheta);
	tree1->SetBranchAddress("ProtonEnergy",		&ProtonEnergy);

	TH1D *c0 = new TH1D("NProton",		"NProton",	  5,	   0,	  5);
	TH1D *c1 = new TH1D("ProtonPhi",	"ProtonPhi",	360,	-180,	180);
	TH1D *c2 = new TH1D("ProtonTheta",	"ProtonTheta",	180,	   0,	180);
	TH1D *c3 = new TH1D("ProtonEnergy",	"ProtonEnergy",	400,	   0,	400);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		c0->Fill(NProton);

		for (Int_t a = 0; a < NProton; a++){

			c1->Fill(ProtonPhi[a]);
			c2->Fill(ProtonTheta[a]);
			c3->Fill(ProtonEnergy[a]);
		}	
	}

	canvas3->cd(1); c0->Draw();
	canvas3->cd(3); c1->Draw();
	canvas3->cd(5); c2->Draw();
	canvas3->cd(7); c3->Draw();

}

void ProtonDistributions2(Char_t* name = "TA2Pi0Compton2012.root") {

	Char_t* filename = Form("/work0/cristina/%s",name);
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NProton;
	Double_t ProtonTheta[100], ProtonPhi[100], ProtonEnergy[100];

	tree1->SetBranchAddress("NProton",		&NProton);
	tree1->SetBranchAddress("ProtonPhi",		&ProtonPhi);
	tree1->SetBranchAddress("ProtonTheta",		&ProtonTheta);
	tree1->SetBranchAddress("ProtonEnergy",		&ProtonEnergy);

	TH1D *c0 = new TH1D("NProton",		"NProton",	  5,	   0,	  5);
	TH1D *c1 = new TH1D("ProtonPhi",	"ProtonPhi",	360,	-180,	180);
	TH1D *c2 = new TH1D("ProtonTheta",	"ProtonTheta",	180,	   0,	180);
	TH1D *c3 = new TH1D("ProtonEnergy",	"ProtonEnergy",	400,	   0,	400);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		c0->Fill(NProton);

		for (Int_t a = 0; a < NProton; a++){

			c1->Fill(ProtonPhi[a]);
			c2->Fill(ProtonTheta[a]);
			c3->Fill(ProtonEnergy[a]);
		}	
	}

	canvas3->cd(2); c0->Draw();
	canvas3->cd(4); c1->Draw();
	canvas3->cd(6); c2->Draw();
	canvas3->cd(8); c3->Draw();

}
