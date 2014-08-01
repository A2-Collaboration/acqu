#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"

void MissingMass() {

	Char_t* filename = Form("/work0/cristina/TA2Pi0Compton2008.root");
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	gStyle->SetOptStat(1111);

	Int_t NPromptPi0, NRandomPi0;

	Double_t MissingMassPromptPi0[1000], MissingMassRandomPi0[1000];

	tree1->SetBranchAddress("NPromptPi0",		&NPromptPi0);
	tree1->SetBranchAddress("NRandomPi0",		&NRandomPi0);

	tree1->SetBranchAddress("MissingMassPromptPi0",	&MissingMassPromptPi0);
	tree1->SetBranchAddress("MissingMassRandomPi0",	&MissingMassRandomPi0);

	TH1D *b0 = new TH1D("Missing-Mass-Pi0-Prompt","Missing-Mass-Pi0-Prompt", 400,800,1200);
	TH1D *b1 = new TH1D("Missing-Mass-Pi0-Random","Missing-Mass-Pi0-Random", 400,800,1200);
	TH1D *b2 = new TH1D("Missing-Mass-Pi0",	      "Missing-Mass-Pi0",	 400,800,1200);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		for (Int_t a = 0; a < NPromptPi0; a++){
				b0->Fill(MissingMassPromptPi0[a]);
		}	

		for (Int_t a = 0; a < NRandomPi0; a++){
				b1->Fill(MissingMassRandomPi0[a]);
		}
	}
	
	b2->Add(b0,1);
//	b2->Add(b1,0);	    // DataMC
	b2->Add(b1,-0.0714);// Data2008
//	b2->Add(b1,-0.025); // Data2012

	canvas1 = new TCanvas("canvas1","Missing-Mass-Pi0-Prompt");
	b0->Draw();

	canvas2 = new TCanvas("canvas2","Missing-Mass-Pi0-Random");	
	b1->Draw();

	canvas3 = new TCanvas("canvas3","Missing-Mass-Pi0-Subtracted");	
	b2->Draw();
}
