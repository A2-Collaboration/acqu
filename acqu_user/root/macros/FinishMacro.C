void FinishMacro(Char_t* file = "ARHistograms.root")
{
  // Stuff to do at the end of an analysis run
  // Here all spectra are saved to disk
  printf("End-of-Run macro executing\n");
  TFile f(file,"recreate");
  if( !f ){
    printf("Open file %s for histogram save FAILED!!\n",file);
    return;
  }
  gROOT->GetList()->Write();
  f.Close();
  printf("All histograms saved to %s\n\n",file);
}
