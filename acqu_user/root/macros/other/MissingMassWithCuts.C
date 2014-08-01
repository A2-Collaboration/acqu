#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"

	TH1D *a2 = new TH1D("Missing-Mass-Compton","Missing-Mass-Compton",200,850,1050);
	TH1D *b2 = new TH1D("Missing-Mass-Pi0",	   "Missing-Mass-Pi0",	  200,850,1050);


void Main(){

	MissingMass(2008); // Get real data
	TH1D *a2012 = new TH1D("Missing-Mass-Compton-2012","Missing-Mass-Compton-2012",200,850,1050);
	TH1D *b2012 = new TH1D("Missing-Mass-Pi0-2012",    "Missing-Mass-Pi0-2012",    200,850,1050);

	a2012->Add(a2,1);
	b2012->Add(b2,1);

	MissingMass(0);   // Get Compton simulated data
	TH1D *aMC_C = new TH1D("Missing-Mass-Compton-MC-C","Missing-Mass-Compton-MC-C",200,850,1050);
	TH1D *bMC_C = new TH1D("Missing-Mass-Pi0-MC-C",    "Missing-Mass-Pi0-MC-C",    200,850,1050);

	aMC_C->Add(a2,1);
	bMC_C->Add(b2,1);

	MissingMass(1);   // Get Pi0 simulated data
	TH1D *aMC_P = new TH1D("Missing-Mass-Compton-MC-P","Missing-Mass-Compton-MC-P",200,850,1050);
	TH1D *bMC_P = new TH1D("Missing-Mass-Pi0-MC-P",    "Missing-Mass-Pi0-MC-P",    200,850,1050);

	aMC_P->Add(a2,1);
	bMC_P->Add(b2,1);

	canvas1 = new TCanvas("canvas1","Missing-Mass-Compton");
	canvas1->Divide(3,1);
	canvas1->cd(1); a2012->Draw();
	canvas1->cd(2); aMC_C->Draw();
	canvas1->cd(3); aMC_P->Draw();

	canvas2 = new TCanvas("canvas2","Missing-Mass-Pi0");
	canvas2->Divide(3,1);
	canvas2->cd(1); b2012->Draw();
	canvas2->cd(2); bMC_C->Draw();
	canvas2->cd(3); bMC_P->Draw();

	THStack *ahs = new THStack("ahs", "Simulated data (Blue = Compton, Purple = Pi0)");
	TH1D *aMC_Cs = new TH1D("Missing-Mass-Compton-MC-C-S","Missing-Mass-Compton-MC-C-S",200,850,1050);
	TH1D *aMC_Ps = new TH1D("Missing-Mass-Compton-MC-P-S","Missing-Mass-Compton-MC-P-S",200,850,1050);

	aMC_Cs->Add(aMC_C,0.0025);
	aMC_Ps->Add(aMC_P,0.25);

	aMC_Cs->SetMarkerStyle(21);
	aMC_Cs->SetMarkerColor(kBlue);
	aMC_Cs->SetFillColor(kBlue);

	aMC_Ps->SetMarkerStyle(21);
	aMC_Ps->SetMarkerColor(kMagenta+1);
	aMC_Ps->SetFillColor(kMagenta+1);

	ahs->Add(aMC_Cs);
	ahs->Add(aMC_Ps);

	THStack *bhs = new THStack("bhs", "Simulated data (Blue = Compton, Purple = Pi0)");
	TH1D *bMC_Ps = new TH1D("Missing-Mass-Pi0-MC-P-S",    "Missing-Mass-Pi0-MC-P-S",    200,850,1050);
	TH1D *bMC_Cs = new TH1D("Missing-Mass-Pi0-MC-C-S",    "Missing-Mass-Pi0-MC-C-S",    200,850,1050);

	bMC_Cs->Add(bMC_C,0.00125);
	bMC_Ps->Add(bMC_P,.125);

	bMC_Cs->SetMarkerStyle(21);
	bMC_Cs->SetMarkerColor(kBlue);
	bMC_Cs->SetFillColor(kBlue);

	bMC_Ps->SetMarkerStyle(21);
	bMC_Ps->SetMarkerColor(kMagenta+1);
	bMC_Ps->SetFillColor(kMagenta+1);

	bhs->Add(bMC_Cs);
	bhs->Add(bMC_Ps);

	canvas3 = new TCanvas("canvas3","Stacked Simulated Data");
	canvas3->Divide(2,1);
	canvas3->cd(1); ahs->Draw(); a2012->Draw("sames");
	canvas3->cd(2); bhs->Draw(); b2012->Draw("sames");

}

void MissingMass(UInt_t DATA = 2012) {

	a2->Reset();
	b2->Reset();

	if (DATA == 2008) {
	Char_t* filename = Form("/work0/cristina/TA2Pi0Compton2008.root");
	}
	else if (DATA == 2012) {
	Char_t* filename = Form("/work0/cristina/TA2Pi0Compton2012_oldcommit.root");
	}
	else if (DATA == 0) { 
	Char_t* filename = Form("/work0/cristina/FromOldComp/TA2Pi0ComptonMC_Compton.root");
	}
	else if (DATA == 1) { 
	Char_t* filename = Form("/work0/cristina/FromOldComp/TA2Pi0ComptonMC_Pi0.root");
	}

	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

	gStyle->SetOptStat(1111);

	if (DATA == 2008) {
	Int_t ScalingFactorCompton = -0.0714; //Data2008
	Int_t ScalingFactorPi0 	   = -0.0714; //Data2008
	}

	else if (DATA == 2012) {
	Int_t ScalingFactorCompton = -0.025;   //Data2012
	Int_t ScalingFactorPi0 	   = -0.025;   //Data2012
	}

	else if (DATA == 0) {
	Int_t ScalingFactorCompton = 0;       //DataMC
	Int_t ScalingFactorPi0 	   = 0;       //DataMC
	}

	else if (DATA == 1) {
	Int_t ScalingFactorCompton = 0;       //DataMC
	Int_t ScalingFactorPi0 	   = 0;       //DataMC
	}

	Int_t NPhoton, NProton, NPi0;
//	Int_t CBNParticle;
//	Double_t CBESum;

	Int_t NPrompt, NRandom;
	Int_t NPromptPi0, NRandomPi0;

	Int_t TaggerChannelPrompt[1000], TaggerChannelRandom[1000];
	Int_t TaggerChannelPromptPi0[1000], TaggerChannelRandomPi0[1000];

	Double_t MissingMassPrompt[1000], MissingMassRandom[1000];
	Double_t MissingMassPromptPi0[1000], MissingMassRandomPi0[1000];

	tree1->SetBranchAddress("NPhoton",		&NPhoton);
	tree1->SetBranchAddress("NProton",		&NProton);
	tree1->SetBranchAddress("NPi0",			&NPi0);

//	tree1->SetBranchAddress("CBNParticle",		&CBNParticle);
//	tree1->SetBranchAddress("CBESum",		&CBESum);

	tree1->SetBranchAddress("NPrompt",		&NPrompt);
	tree1->SetBranchAddress("NRandom",		&NRandom);

	tree1->SetBranchAddress("NPromptPi0",		&NPromptPi0);
	tree1->SetBranchAddress("NRandomPi0",		&NRandomPi0);

	tree1->SetBranchAddress("TaggerChannelPrompt",	&TaggerChannelPrompt);
	tree1->SetBranchAddress("TaggerChannelRandom",	&TaggerChannelRandom);

	tree1->SetBranchAddress("TaggerChannelPromptPi0",&TaggerChannelPromptPi0);
	tree1->SetBranchAddress("TaggerChannelRandomPi0",&TaggerChannelRandomPi0);

	tree1->SetBranchAddress("MissingMassPrompt",	&MissingMassPrompt);
	tree1->SetBranchAddress("MissingMassRandom",	&MissingMassRandom);

	tree1->SetBranchAddress("MissingMassPromptPi0",	&MissingMassPromptPi0);
	tree1->SetBranchAddress("MissingMassRandomPi0",	&MissingMassRandomPi0);

	TH1D *a0 = new TH1D("Missing-Mass-Prompt", "Missing-Mass-Prompt", 200,850,1050);
	TH1D *a1 = new TH1D("Missing-Mass-Random", "Missing-Mass-Random", 200,850,1050);

	TH1D *b0 = new TH1D("Missing-Mass-Pi0-Prompt","Missing-Mass-Pi0-Prompt", 200,850,1050);
	TH1D *b1 = new TH1D("Missing-Mass-Pi0-Random","Missing-Mass-Pi0-Random", 200,850,1050);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		Bool_t KinCutsCompton == kFALSE;
		if (NPhoton == 1 && NProton == 1){
//		if (NPhoton == 1 && NProton == 1 && CBNParticle == 2){
			KinCutsCompton = kTRUE;
		}

		Bool_t KinCutsPi0 == kFALSE;
		if (NPi0 == 1 && NProton == 1){
//		if (NPi0 == 1 && NProton == 1 && CBNParticle == 3){
			KinCutsPi0 = kTRUE;
		}

		Bool_t Trigger == kFALSE;
//		if (CBESum > 80.0){
			Trigger = kTRUE;
//		}

		if (Trigger == kTRUE) {

			if (DATA == 2012 || DATA == 2008) {

			if (KinCutsCompton == kTRUE){
				for (Int_t a = 0; a < NPrompt; a++){
					if (TaggerChannelPrompt[a] >= 240 && TaggerChannelPrompt[a] <= 260){
						a0->Fill(MissingMassPrompt[a]);
					}
				}	

				for (Int_t a = 0; a < NRandom; a++){
					if (TaggerChannelRandom[a] >= 240 && TaggerChannelRandom[a] <= 260){
						a1->Fill(MissingMassRandom[a]);
					}
				}	
			}

			if (KinCutsPi0 == kTRUE){
				for (Int_t a = 0; a < NPromptPi0; a++){
					if (TaggerChannelPromptPi0[a] >= 240 && TaggerChannelPromptPi0[a] <= 260){
						b0->Fill(MissingMassPromptPi0[a]);
					}
				}	
	
			for (Int_t a = 0; a < NRandomPi0; a++){
					if (TaggerChannelRandomPi0[a] >= 240 && TaggerChannelRandomPi0[a] <= 260){
						b1->Fill(MissingMassRandomPi0[a]);
					}	
				}
			}
			}

			else if (DATA == 0 || DATA == 1) {

				if (KinCutsCompton == kTRUE){

				for (Int_t a = 0; a < NPrompt; a++){
						a0->Fill(MissingMassPrompt[a]);
				}	

				for (Int_t a = 0; a < NRandom; a++){
						a1->Fill(MissingMassRandom[a]);
				}	
			}

			if (KinCutsPi0 == kTRUE){

				for (Int_t a = 0; a < NPromptPi0; a++){
						b0->Fill(MissingMassPromptPi0[a]);
				}	

				for (Int_t a = 0; a < NRandomPi0; a++){
						b1->Fill(MissingMassRandomPi0[a]);
				}
			}
			}

		}

	}
	
	a2->Add(a0,1);
	a2->Add(a1,ScalingFactorCompton); 

	b2->Add(b0,1);
	b2->Add(b1,ScalingFactorPi0); 

}


