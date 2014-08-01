
/* -dk- 26.02.2004: this first version.
           Attempt to create shiftfile for tagger TDCs.
	   
	   Loop through tagger TDCs declared as elements in
	   in_file="FPladder.cal", calulate shift form offet and mean
	   and save new FPladder.cal-like file to out_file

	   - Scaler-Index is not supported.

	   To Do:
	   * cleanup
	   * make global offset (-800) a parameter or constant
	   * buffer output such, that in/out-file can be the same (.sav!)
	   * crate .root-file with "output" spectra.
	   * add some error-handling.
*/

int TaggCal(char *in_file, char* out_file){
  ifstream f_in;                     // input-file
  FILE *f_out;                       // outut-file
  char line[255];                    // current f_in-line
  UInt_t ncols;                      // #columns in current f_in-line
  char tag[255], A1[255], Ecal[255]; // i/o data. See FPladder.cal
  UInt_t    TDC, LTh, HTh;           //  for more info.
  Double_t  A0;
  
  char hist_name[20];                // name of histo to look at
  TH1F *hist;                        // the histo itself

  // the canvas; disable this and Draw()/Update() to speed up this macro
  tdc_canv = new TCanvas("tdc_canv", "Tagger: RAW TDCs", 
			 200, 60, 700, 900);  
  
  f_out = fopen(out_file, "w");
  f_in.open(in_file);
  
  if (!f_in) {printf("<E> TagCal - File not found\n"); return -1;}
  while(1) {
    if (!f_in.good()) break;
    f_in.getline(line, 255);
    // input format: 
    // " Element: TDC#  LTh HTh   A0    A1   Ecal     Scaler#"
    // think aobout takin out the spaces in format-string.
    ncols = sscanf(line, "%s %d %d %d %lf %s %s", 
		   tag, &TDC, &LTh, &HTh, &A0, A1, Ecal);
    if( ncols == 7 && strncmp(tag, "Element:", (size_t)8) == 0 ){
      sprintf(hist_name, "ADC%d", TDC);
      hist = (TH1F)gROOT->Get(hist_name);
      // --> This is the main thing to do
      A0 = -800 +1 * hist->GetMean(); // get weighted mean
      // <--
      fprintf(f_out, " %s  %d  %d  %d  %8.3lf  %s  %s\n", 
	      tag, TDC, LTh, HTh, A0, A1, Ecal);
      hist->Draw();                  // disable to speed up the macro
      tdc_canv->Update();            // disable to speed up the macro
    } else
      {
	fprintf(f_out, "%s\n", line);
      }
  }

  f_in.close();
  fclose(f_out);
  return 0;
}
