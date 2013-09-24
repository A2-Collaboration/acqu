void FinishTree(Char_t* fileAR = NULL)
{
	TString name;
	printf("\nEnd-of-Run macro executing:\n");

        printf("Closing tree files...");
	TA2Pi0Compton* comp = (TA2Pi0Compton*)(gAN->GetPhysics());
	comp->CloseTrees();
	printf(" done. \n");

	printf("Saving defined histograms...");
	if( !fileAR) fileAR = "ARHistograms.root";
	TFile f1(fileAR,"RECREATE");
	gROOT->GetList()->Write();
	f1.Close();
	printf("done.\n",fileAR);
	printf("All histograms saved to %s\n",fileAR);

	// Get some histos from full histogram file and save them into tree file
	Int_t IsTreeFileOn  = ((TA2Pi0Compton*)gAR->GetAnalysis()->GetPhysics())->IsTreeFileOn();
	if (IsTreeFileOn == 1) {
 
		TFile f1(fileAR,  "READ"); 

		TH1D *FPD_ScalerAcc 	 = (TH1D*)f1.Get("FPD_ScalerAcc");	FPD_ScalerAcc->SetDirectory(0);
		TH1D *SumScalers268to619 = (TH1D*)f1.Get("SumScalers268to619"); SumScalers268to619->SetDirectory(0);
		TH1D *Trigger		 = (TH1D*)f1.Get("Trigger");		Trigger->SetDirectory(0);

		f1.Close();

		Char_t* fileTree = ((TA2Pi0Compton*)gAR->GetAnalysis()->GetPhysics())->GetTreeFileName();
		TFile f2(fileTree,"UPDATE");
		
		FPD_ScalerAcc->Write();
		SumScalers268to619->Write();
		Trigger->Write();

		f2.Close();		
	}
}


void FinishTreeWorker()
{
	// Establish file name from Input file name
	Char_t namehold[256]; Int_t RUN; 
	Char_t* fInputName = gAR->GetFileName();
	sscanf( fInputName, "scratch/%[^_]_%d.dat\n", namehold, &RUN);
	Char_t* fileAR = Form("output/ARHistograms_%s_%d.root",namehold,RUN);

    	TString name;
	printf("\nEnd-of-Run macro executing:\n");

        printf("Closing tree files..."); 
 	TA2Pi0Compton* comp = (TA2Pi0Compton*)(gAN->GetPhysics());
 	comp->CloseTrees();
	printf(" done. \n");

	printf("Saving defined histograms...");
	TFile f1(fileAR,"RECREATE");
	gROOT->GetList()->Write();
	f1.Close();
  	printf("done.\n",fileAR);
  	printf("All histograms saved to %s\n",fileAR);

	// Get some histos from full histogram file and save them into tree file
	Int_t IsTreeFileOn  = ((TA2Pi0Compton*)gAR->GetAnalysis()->GetPhysics())->IsTreeFileOn();
	if (IsTreeFileOn == 1) {
 
		TFile f1(fileAR,  "READ"); 

		TH1D *FPD_ScalerAcc 	 = (TH1D*)f1.Get("FPD_ScalerAcc");	FPD_ScalerAcc->SetDirectory(0);
		TH1D *SumScalers268to619 = (TH1D*)f1.Get("SumScalers268to619"); SumScalers268to619->SetDirectory(0);
		TH1D *Trigger		 = (TH1D*)f1.Get("Trigger");		Trigger->SetDirectory(0);

		f1.Close();

		Char_t* fileTree = ((TA2Pi0Compton*)gAR->GetAnalysis()->GetPhysics())->GetTreeFileName();
		TFile f2(fileTree,"UPDATE");
		
		FPD_ScalerAcc->Write();
		SumScalers268to619->Write();
		Trigger->Write();

		f2.Close();		
	}
}

