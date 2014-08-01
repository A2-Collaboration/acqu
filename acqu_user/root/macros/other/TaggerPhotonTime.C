#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"

void TaggerPhotonTime() {

	Char_t* filename = Form("TA2Cristina.root");
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("TA2CristinaTree");

	gStyle->SetOptStat(1111);

	Int_t NTagg, NPhoton;
	Int_t TaggerChannel[352];
	Double_t TaggerPhotonTime[2000];

	tree1->SetBranchAddress("NTagg",		&NTagg);
	tree1->SetBranchAddress("NPhoton",		&NPhoton);

	tree1->SetBranchAddress("TaggerChannel",	&TaggerChannel);
	tree1->SetBranchAddress("TaggerPhotonTime",	&TaggerPhotonTime);

	TH2F *a0 = new TH2F("TaggerPhotonTime","TaggerPhotonTime",352,0,352,600,-300,300);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		Int_t NTaggNPhot = 0;
		for (Int_t a = 0; a < NPhoton; a++){
			for (Int_t b = 0; b < NTagg; b++){

				a0->Fill(TaggerChannel[b],TaggerPhotonTime[NTaggNPhot]);
				NTaggNPhot++;		
			}
		}	
	}
	
	canvas1 = new TCanvas("canvas1","TaggerPhotonTime");
	a0->Draw("surf2");

}

void TaggerTime() {

	Char_t* filename = Form("TA2Cristina.root");
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("TA2CristinaTree");

	gStyle->SetOptStat(1111);

	Int_t NTagg, NPhoton;
	Int_t TaggerChannel[352];
	Double_t TaggerTime[352];
	Double_t PhotonTime[20];
	Double_t TaggerPhotonTime;

	tree1->SetBranchAddress("NTagg",		&NTagg);
	tree1->SetBranchAddress("NPhoton",		&NPhoton);

	tree1->SetBranchAddress("TaggerChannel",	&TaggerChannel);
	tree1->SetBranchAddress("TaggerTime",		&TaggerTime);
	tree1->SetBranchAddress("PhotonTime",		&PhotonTime);

	TH2F *a0 = new TH2F("TaggerPhotonTime","TaggerPhotonTime",352,0,352,600,-300,300);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		Int_t NTaggNPhot = 0;
		for (Int_t a = 0; a < NPhoton; a++){
			for (Int_t b = 0; b < NTagg; b++){
				TaggerPhotonTime = TaggerTime[b] - PhotonTime[a];
				a0->Fill(TaggerChannel[b],TaggerPhotonTime);
				NTaggNPhot++;		
			}
		}	
	}
	
	canvas1 = new TCanvas("canvas1","TaggerPhotonTime");
	a0->Draw("surf2");

}


