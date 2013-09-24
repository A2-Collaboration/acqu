//###################################################################################################
//
// Asymmetry macro, Author: Cristina Collicott
//
// Subroutines:
// Fill3DHist -- Fills random subtracted 3D histogram 
// 		 (Missing Mass vs. Tagger Channel vs. Pi0 Theta)
//
// CalculateYield -- Integrates 3D histogram over specified range of
// 		  Missing mass, Tagger channels and Pi0 Theta.	 
//
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


//###################################################################################################
// Files
//###################################################################################################

	Char_t* PERP = Form("TA2Pi0Compton2012_PERP.root");
	Char_t* PARA = Form("TA2Pi0Compton2012_PARA.root");
	Char_t* SAVE = Form("Asymmetry.root");

//###################################################################################################
// Shared Variables
//###################################################################################################
	Int_t yield;

void DisplayProjection() {

	// Look up asym histograms
	TFile f(SAVE,"READ"); 
	if (f.IsOpen()) {
		printf("Asymmetry histograms found -- using %s\n",SAVE);

		TH3D *asym = (TH3D*)f.Get("asym");
		asym->SetDirectory(0);
		f.Close();
	}	
	else {
		printf("Asymmetry histograms not found -- generating from scratch\n",SAVE);
		Asymmetry();
	} 

}

void Asymmetry(Int_t TC_binsize = 5, Int_t Theta_binsize = 20) {

	// Set default parameters
	Int_t MM_low = 900; 	Int_t MM_high = 1000; 
	Int_t TC_min = 240; 	Int_t TC_max  = 250;
	Int_t Theta_min = 0; 	Int_t Theta_max = 180;
	
	// Look up, or generate 3d histograms
	TFile f(SAVE,"READ"); 
	if (f.IsOpen()) {
		printf("Previously generated 3d histograms found -- using %s\n",SAVE);

		TH3D *perp3d = (TH3D*)f.Get("perp3d");
		TH3D *para3d = (TH3D*)f.Get("para3d");	

		perp3d->SetDirectory(0);
		para3d->SetDirectory(0);
		
		f.Close();
	}	
	else {
		printf("Previously generated 3d histograms not found -- generating now... \n");

		TH3D *perp3d  	= new TH3D("perp3d", 	"perp3d",	1200,0,1200,   	352,0,352,	180,0,180);
		TH3D *para3d  	= new TH3D("para3d", 	"para3d",  	1200,0,1200,   	352,0,352,  	180,0,180);
	
		printf("Filling 3D histograms for PERP file... \n");
		Fill3DHist(PERP, -0.025, perp3d);

		printf("Filling 3D histograms for PARA file... \n");
		Fill3DHist(PARA, -0.025, para3d);
	}

	// Calculate Asymmetry
	Int_t 	TC_low, TC_high, Theta_low, Theta_high;
	Int_t 	Y_perp, Y_para, Diff, Sum;
	Double_t Asymmetry, Theta_mid, TC_mid;

	TH3D *asym 	= new TH3D("asym",   	"asym", 	352,0,352,  	180,0,180,	1000, -0.2, 0.2);

	// Create histograms (and clean up any asymmetries previously produced)
	TFile f(SAVE,"UPDATE"); 

	for (Int_t i = 0; i < 20; i++) {
		Char_t* keyname  = Form("TC_histo_%d;1", i);
		TKey *key = f.FindKey(keyname);
		if (key !=0) f.Delete(keyname);

		Char_t* keyname  = Form("Theta_histo_%d;1", i);
		TKey *key = f.FindKey(keyname);
		if (key !=0) f.Delete(keyname);

		Char_t* keyname  = Form("asym;1");
		TKey *key = f.FindKey(keyname);
		if (key !=0) f.Delete(keyname);

	}
	f.Close();
	
	Int_t *TaggerChannel[20]; Int_t NTC = 0;
	for (Int_t i = TC_min; i < TC_max; i+=TC_binsize){ TaggerChannel[NTC] = i; NTC++; }

	Int_t *ThetaBin[20]; Int_t NTheta = 0;
	for (Int_t i = Theta_min; i < Theta_max; i+=Theta_binsize){ ThetaBin[NTheta] = i; NTheta++; }

	// Loop over tagger channels and calculate yield, store in histograms
	Double_t  FoundTaggerChannels[200],	 FoundThetaBins[200], 	 FoundAsym[200];
	Double_t dFoundTaggerChannels[200],	dFoundThetaBins[200], 	dFoundAsym[200];
	Double_t dAsymmetry;
	TGraphErrors *AsymByTaggerChannel[20];

	Int_t NAsymPoints = 0;

	for (Int_t i = 0; i < NTC; i++){

		TC_low  = TaggerChannel[i];
		TC_high = TC_low + TC_binsize - 1;
		TC_mid  = double(TC_high + TC_low)/2.0;

		for (Int_t j = 0; j < NTheta; j++) {

			Theta_low  = ThetaBin[j];
			Theta_high = Theta_low + Theta_binsize - 1; // Might miss 180
			Theta_mid  = double(Theta_high + Theta_low)/2.0;

//printf("j = %d, ThetaBin[j] = %d, Theta_binsize = %d, Theta_low = %d, Theta_high = %d, Theta_mid =%f\n",j, ThetaBin[j], Theta_binsize, Theta_low, Theta_high, Theta_mid);

			CalculateYield(perp3d, MM_low, MM_high, TC_low, TC_high, Theta_low, Theta_high); Y_perp = yield;
			CalculateYield(para3d, MM_low, MM_high, TC_low, TC_high, Theta_low, Theta_high); Y_para = yield;

			Diff = Y_perp - Y_para;
			Sum  = Y_perp + Y_para;
			Asymmetry = double(Diff)/double(Sum);
			dAsymmetry = (((((Y_perp)**(0.5) + (Y_para)**(0.5) )/(Y_perp - Y_para))**2) + ((((Y_perp)**(0.5) + (Y_para)**(0.5))/(Y_perp + Y_para))**2))**(0.5);
			printf("%f\n",dAsymmetry);

			if (Sum != 0) {
//				printf("TC: %f, Theta: %f, PERP: %d, PARA: %d, Diff: %d, Sum: %d, Asym: %f \n", TC_mid, Theta_mid, Y_perp, Y_para, Diff, Sum, Asymmetry);

				asym->Fill(TC_mid,Theta_mid,Asymmetry);

				FoundTaggerChannels[NAsymPoints] = TC_mid;
				FoundThetaBins[NAsymPoints] 	 = Theta_mid;
				FoundAsym[NAsymPoints] 	 	 = Asymmetry;

				dFoundTaggerChannels[NAsymPoints] = double(TC_binsize)/2.0;
				dFoundThetaBins[NAsymPoints] 	  = double(Theta_binsize)/2.0;
				dFoundAsym[NAsymPoints] 	  = dAsymmetry;

				NAsymPoints++;
			}
		}
	}

	Double_t StoreTaggerChannels[200],   StoreThetaBins[200],  StoreAsym[200]; Int_t Count;
	Double_t dStoreTaggerChannels[200], dStoreThetaBins[200], dStoreAsym[200]; 
	Double_t TC_compare;

	for (Int_t i = 0; i < NTC; i++) {
		Count = 0;
		TC_compare = double(TaggerChannel[i]) + (double(TC_binsize - 1)/2.0);
		for (Int_t j = 0; j < NAsymPoints; j++) {
			if (FoundTaggerChannels[j] == TC_compare){
				StoreThetaBins[Count] 	= FoundThetaBins[j];
				StoreAsym[Count] 	= FoundAsym[j];

				dStoreThetaBins[Count] 	= dFoundThetaBins[j];
				dStoreAsym[Count] 	= dFoundAsym[j];

				Count++;
			}
		}
		printf("%d\n",Count);
		AsymByTaggerChannel[i] = new TGraphErrors(Count,StoreThetaBins,StoreAsym,dStoreThetaBins, dStoreAsym);	
	}

	AsymByTaggerChannel[0]->Draw("ALP");

//	for (Int_t i = 0; i < (NTC-1); i++) { TChist[i]->Write(); }
//	for (Int_t i = 0; i < (NTheta-1); i++) { Thetahist[i]->Write(); }
	
//	asym->Write();
//	f.Close(); 

}

void Fill3DHist(Char_t* filename, Double_t Ratio, TH3D* h) {

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

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		for (Int_t a = 0; a < NPromptPi0; a++){
			h0->Fill(MissingMassPromptPi0[a], TaggerChannelPromptPi0[a], Pi0ThetaPrompt[a]);
		}	

		for (Int_t a = 0; a < NRandomPi0; a++){
			h1->Fill(MissingMassRandomPi0[a], TaggerChannelRandomPi0[a], Pi0ThetaRandom[a]);
		}
	}
	
	h->Add(h0,1);
	h->Add(h1,Ratio);

	TFile f(SAVE,"UPDATE"); 
		h->Write();
	f.Close(); 

	h0->Delete();
	h1->Delete();

}

void CalculateYield(TH3D* data, Int_t MM_low, Int_t MM_high, Int_t TC_low, Int_t TC_high, Int_t Theta_low, Int_t Theta_high){

	h = (TH3D*) data->Clone();
	h->GetYaxis()->SetRange( TC_low, TC_high);
	h->GetZaxis()->SetRange( Theta_low, Theta_high);
	proj = (TH1D*) h->Project3D("x");
	
	yield = proj->Integral(MM_low,MM_high);

	proj->Delete();
	h->Delete();
}
