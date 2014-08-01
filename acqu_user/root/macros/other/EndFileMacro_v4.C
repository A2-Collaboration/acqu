void EndFileMacro()
{
  // Stuff to do at the end of a data file
  // Here all spectra are saved to disk
  printf("End-of-File macro executing\n");
  Char_t* n = strrchr( gAR->GetFileName(), '/' ) + 1;
  Char_t name[128];
  strcpy( name, n );
  Char_t* n1 = strrchr( name, '.' );
  *n1 = '\0';
  strcat( name, "_hist.root" );
  TFile f( name,"recreate");
  gROOT->GetList()->Write();
  f.Close();
  gAN->ZeroAll();
  printf("All histograms saved to %s and zero'ed\n\n", name);
}
