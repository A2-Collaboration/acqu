//###################################################################################################
//
// Cross Section macro, Author: Cristina Collicott
//
// Subroutines:
// Fill3DHist -- Fills random subtracted 3D histogram 
// 		 (Missing Mass vs. Tagger Channel vs. Pi0 Theta)
//
// CalculateYield -- Integrates 3D histogram over specified range of
// 		  Missing mass, Tagger channels and Pi0 Theta.	 
//
// CalculateDeff -- Calculates yield/total events for a simulated file.
//
// GetTarget -- Returns Target density (for LH2 target, 10cm) in protons/cm^2
//
// ##################################################################################################

#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"

	gStyle->SetOptStat(1111);

	TH3D *h2    = new TH3D("hist",  "hist",	 1200,0,1200,   352,0,352,  180,0,180);
	TH3D *data  = new TH3D("data",  "data",  1200,0,1200,   352,0,352,  180,0,180);
	TH3D *sim   = new TH3D("sim",   "sim",   1200,0,1200,   352,0,352,  180,0,180);

	TH1D *projection = new TH1D("Projection","Projection",1200,0,1200);

	Int_t yield, nentries;
	Double_t Deff, Target;

void CrossSection(){

	Char_t* filename = Form("/work0/cristina/TA2Pi0Compton2012.root");
	Fill3DHist(filename, -0.025);
	data = (TH3D*) h2->Clone();

	Int_t MM_low, MM_high, TC_low, TC_high, Theta_low, Theta_high;
	MM_low = 900; MM_high = 1000; TC_low = 250; TC_high = 300; Theta_low = 0; Theta_high = 180;

	for (Int_t TC = TC_low; TC <= TC_high; TC++){
		CalculateYield(MM_low, MM_high, TC, TC, Theta_low, Theta_high);
	}
}

void Fill3DHist(Char_t* filename, Double_t Ratio) {

	h2->Clear();

	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	Int_t NPromptPi0, NRandomPi0;
	Double_t MissingMassPromptPi0[1000],  	MissingMassRandomPi0[1000];
	Int_t TaggerChannelPromptPi0[1000], 	TaggerChannelRandomPi0[1000];
	Double_t Pi0ThetaPrompt[1000],      	Pi0ThetaRandom[1000];

	tree1->SetBranchAddress("NPromptPi0",		&NPromptPi0);
	tree1->SetBranchAddress("NRandomPi0",		&NRandomPi0);

	tree1->SetBranchAddress("MissingMassPromptPi0",	&MissingMassPromptPi0);
	tree1->SetBranchAddress("MissingMassRandomPi0",	&MissingMassRandomPi0);

	tree1->SetBranchAddress("TaggerChannelPromptPi0",&TaggerChannelPromptPi0);
	tree1->SetBranchAddress("TaggerChannelRandomPi0",&TaggerChannelRandomPi0);

	tree1->SetBranchAddress("Pi0ThetaPrompt",	&Pi0ThetaPrompt);
	tree1->SetBranchAddress("Pi0ThetaRandom",	&Pi0ThetaRandom);

// Missing Mass vs. Tagger Channel vs. Pi0 Theta

	TH3D *h0 = new TH3D("Hist-Prompt","Hist-Prompt", 1200,0,1200, 	352,0,352,  180,0,180);
	TH3D *h1 = new TH3D("Hist-Random","Hist-Random", 1200,0,1200, 	352,0,352,  180,0,180);

	nentries = (Int_t)tree1->GetEntries();
	nentries = 10000;
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		for (Int_t a = 0; a < NPromptPi0; a++){
			h0->Fill(MissingMassPromptPi0[a], TaggerChannelPromptPi0[a], Pi0ThetaPrompt[a]);
		}	

		for (Int_t a = 0; a < NRandomPi0; a++){
			h1->Fill(MissingMassRandomPi0[a], TaggerChannelRandomPi0[a], Pi0ThetaRandom[a]);
		}
	}
	
	h2->Add(h0,1);
	h2->Add(h1,Ratio);

}

void CalculateYield(Int_t MM_low, Int_t MM_high, Int_t TC_low, Int_t TC_high, Int_t Theta_low, Int_t Theta_high){

	projection->Clear();

	h = (TH3D*) data->Clone();

	h->GetYaxis()->SetRange( TC_low, TC_high);
	h->GetZaxis()->SetRange( Theta_low, Theta_high);
	projection = (TH1D*) h->Project3D("x");
	
	yield = projection->Integral(MM_low,MM_high);
	printf("TC: %d - %d , Yield: %d \n", TC_low, TC_high, yield);

}


void CalculateDeff(Int_t MM_low, Int_t MM_high, Char_t* filename, Double_t Ratio) {

	CalculateYield(MM_low, MM_high, filename, Ratio);
	Deff = Double_t(yield)/Double_t(nentries);
	printf("Detection Efficiency: %f \n", Deff);

}

void GetTarget() {

	Target = 4.249*10**(23);

}

