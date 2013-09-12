void FinishTree(Char_t* file = NULL)
{
	TString name;
	printf("\nEnd-of-Run macro executing:\n");

        printf("Closing tree files...");
	TA2Pi0Compton* comp = (TA2Pi0Compton*)(gAN->GetPhysics());
	comp->CloseTrees();
	printf(" done. \n");

	printf("Saving defined histograms...");
	if( !file) file = "ARHistograms.root";
	TFile f1(file,"RECREATE");
	gROOT->GetList()->Write();
	f1.Close();
	printf("done.\n",file);
	printf("All histograms saved to %s\n",file);
}


void FinishTreeWorker()
{
	// Establish file name from Input file name
	Char_t namehold[256]; Int_t RUN; 
	Char_t* fInputName = gAR->GetFileName();
	sscanf( fInputName, "scratch/%[^_]_%d.dat\n", namehold, &RUN);
	Char_t* file = Form("output/ARHistograms_%s_%d.root",namehold,RUN);

    	TString name;
	printf("\nEnd-of-Run macro executing:\n");

        printf("Closing tree files..."); 
 	TA2Pi0Compton* comp = (TA2Pi0Compton*)(gAN->GetPhysics());
 	comp->CloseTrees();
	printf(" done. \n");

	printf("Saving defined histograms...");
	TFile f1(file,"RECREATE");
	gROOT->GetList()->Write();
	f1.Close();
  	printf("done.\n",file);
  	printf("All histograms saved to %s\n",file);
}
