void SaveOpenSpectras() {
char *spectra_names[] =    
	{"CB-Spectra",
	"TAPS-BaF2-Spectra",
        "TAPS-PbWO4-Spectra",
        "TAPS-Veto-Spectra",
	"Tagger-Spectra",
	"EPT-Spectra",
	"PID-Spectra",
	"MWPC-Spectra",
	"Scaler-Spectra",
	"Livetimes-Spectra",
	"Physics-Spectra",
	0};

	printf("\nDelete all PNG images on a2kallisti:/home/a2cb/OnlineSpectra\n");

	stringstream CmdStr;
        CmdStr << "rm /home/a2cb/acqu/OnlineSpectra/*.png";
        system(CmdStr.str().c_str());


//Save Spectra
UInt_t i = 0;
char spectra_title[50];
char spectra_filename[100];
while (spectra_names[i]){
	sprintf(spectra_title, "%s", spectra_names[i]);
	sprintf(spectra_filename, "/home/a2cb/acqu/OnlineSpectra/%i %s.png", i, spectra_names[i]);
   gSystem->ProcessEvents();

   if (gROOT->FindObject(spectra_title)) {

   	TImage *img = TImage::Create();
   	img->FromPad((TVirtualPad*)gROOT->FindObject(spectra_title));
   	img->WriteImage(spectra_filename);
	delete img;
	printf("Spectra saved: %s\n", spectra_title);
   } else {
	printf("Spectra not found: %s\n", spectra_title);
   }
	i++;
}


	stringstream CmdStr2;
        CmdStr2 << "scp /home/a2cb/acqu/OnlineSpectra/*.png a2cb@macrobusy:/home/a2cb/";
        system(CmdStr2.str().c_str());

	printf("Finished saving all open spectras as PNG on a2cb@macrobusy:/home/a2cb/\n\n");

}
