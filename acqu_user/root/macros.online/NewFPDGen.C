void NewFPDGen(Int_t tdcStart=800, Int_t scalerStart=2000,
	       Int_t n=328, Char_t* calFile=NULL, Double_t eBeam=0.0)
{
  // Generate New FPD ladder file. Input parameters and default values:
  // tdcstart    start index of TDCs
  // scalerStart start index of Scalers
  // n           number of ladder elements
  // calFile     name of electron energy calibration file 
  // eBeam       e- beam energy
  //
  // JRM Annand 8th Sept 2017 
  //
  // Ordering of TDC sections for V1190A
  Int_t OrderTDC[] =
    {19,  2,  3, 18, 17,  0,  1, 16, 23,  6,  7, 22, 21,  4,  5, 20,
     27, 10, 11, 26, 25,  8,  9, 24, 31, 14, 15, 30, 29, 12, 13, 28};
  Int_t OrderScaler[] =
    { 7,  4,  6,  5,  3,  0,  2,  1, 15, 12, 14, 13, 11,  8, 10,  9,
     23, 20, 22, 21, 19, 16, 18, 17, 31, 28, 30, 29, 27, 24, 26, 25};
  Int_t BaseScaler[] = {864, 64, 576, 288, 0, 352, 640, 32, 320, 608, 896};
  //
  Char_t* file_name = "/home/a2cb/acqu/data/NewFPDGen.dat"; // output file name
  Double_t th = 700.0;                 // time upper limit
  Double_t tl = -700.0;               // time lower limit
  Double_t tc = 0.10 ;                 // TDC ns/channel
  Double_t tp = -4050;                 // TDC zero point channel;
  Double_t eCal = 0.0;
  Double_t goosey, eWidth;
  Int_t tdc,scaler;
  FILE* lfile;
  FILE* cfile;
  //
  // Open output file name
  printf("Tagger New FPD ladder configuration saved to file %s\n",file_name);
  if( (lfile = fopen(file_name,"w")) == NULL ){
    printf("File open failed\n");
    return;
  }
  // Open energy calibration file if name supplied
  if( calFile ){
    if( (cfile = fopen(calFile,"r")) == NULL ){
      printf("File %s open failed\n", calFile);
    }
    else printf("Using electron-energy calibration file %s\n", calFile);
  }
  else cfile = NULL;
  //
  // Print some standard lines    
  fprintf(lfile,"##  Auto-Generated Tagger-Ladder Config File\n");
  fprintf(lfile,"##  ROOT Macro NewFPDGen.C\n");
  fprintf(lfile,"##  J.R.M. Annand 8th September 2017\n");
  fprintf(lfile,"Size: %d 0 1 1\n", n);
  fprintf(lfile,"Misc: 1 1 0 0 0   0.0	0.0\n");
  fprintf(lfile,"Window:   -20. 20. -40. -20.\n");
  fprintf(lfile,"Window:   -20. 20.  20.  40.\n");
  fprintf(lfile,"##\n## N.B. No QDC information supplied\n");
  if(cfile)
    fprintf(lfile,"## Electron energy calibration file: %s\n",calFile);
  else
    fprintf(lfile,"## No electron energy calibration supplied\n");
  //
  Int_t ielem = 0;       // element counter
  Int_t itdc;            // TDC ID counter
  Int_t isec;            // section counter
  Int_t k = 1;
  Int_t j, jmax;
  for(itdc=0,isec=0;; itdc+=32,isec++){
    if( k < 4 )
      jmax = 32;
    else{
      jmax = 24;
      k = 0;
    }
    fprintf(lfile,"##\n## Section: %d, FPD Channels: %d - %d  TDC Channels: %d - %d\n##\n",
	    isec, ielem, ielem+jmax-1, itdc+tdcStart, itdc+tdcStart+31);
    for (j=0; j<jmax; j++){
      tdc = tdcStart + OrderTDC[j] + itdc;
      scaler = scalerStart + BaseScaler[isec] + OrderScaler[j];
      if( cfile ){
	fscanf( cfile, "%lf%lf%lf", &goosey, &eCal, &eWidth );
	if( (Int_t)goosey != (ielem+1) ){
	  printf(" Data out of synch in calibration file %s %d %d\n",
		 calFile, goosey, i );
	  return;
	}
	eCal = eBeam - eCal;
      }
      fprintf(lfile,
"Element: %3d 0.0 0.0 0.0 0.0   %4dM0 %.1f %.1f %.1f %.4f 0 0 0 %6.1f 0.0 %d\n",
	      ielem,tdc,tl,th,tp,tc,eCal,scaler);
      if( ielem >= n-1 ) break;
      ielem++;
    }
    if(jmax == 24)
      fprintf(lfile,"## Last 8 TDC channels reserved for reference TDC and misc. applications\n");
    k++;
    if( ielem >= n-1 ) break;
  }
  //
  // Print some standard display lines
  fprintf(lfile,"Initialise: \n");
  fprintf(lfile,"Display: 1D Nhits       64 0 64\n");
  fprintf(lfile,"Display: 1D NhitsPrompt 64 0 64\n");
  fprintf(lfile,"Display: 1D NhitsRand   64 0 64\n");
  fprintf(lfile,"Display: 1D Hits        %d 0 %d\n", n,n);
  fprintf(lfile,"Display: 1D HitsPrompt  %d 0 %d\n", n,n);
  fprintf(lfile,"Display: 1D HitsRand    %d 0 %d\n", n,n);
  fprintf(lfile,"Display: 1D TimeOR     1000 %.1f %.1f\n", tl,th);
  fprintf(lfile,"#Display: 1D Time  1000 %.1f %.1f   0 %d\n", tl,th,n);
  fprintf(lfile,"Display: 1D Eelect   1000 0 %.1f\n", eBeam);
  fprintf(lfile,"Display: 1D ScalerCurr %d 0 %d\n",n,n);
  fprintf(lfile,"Display: 1D ScalerAcc  %d 0 %d\n",n,n);
  fprintf(lfile,"Display: 1D Doubles  	%d 0 %d\n",n,n);
  fprintf(lfile,"Display: 1D NDoubles  	64 0 64\n");
  fprintf(lfile,"Display: 2D TimeOR   200 %.1f %.1f -1 Hits %d 0 %d -1\n",tl,th,n,n);
  fprintf(lfile,"Display: 2D TimeOR   200 %.1f %.1f -1 Doubles %d 0 %d -1\n",tl,th,n,n);
  fclose(lfile);
  if(cfile)fclose(cfile);
  return;
}
//                                           
