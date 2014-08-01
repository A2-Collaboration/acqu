void SaveCut(Char_t* name, Char_t* file = NULL)
{
  // Save a polygon cut made with the graphical editor
  //
  TCutG* cut = (TCutG*)gPad->GetPrimitive("CUTG");
  if( !cut ){
    printf(" No graphical cut defined for displayed pad\n");
    return;
  }
  cut->SetName( name );
  cut->SetTitle(" AcquRoot 2D polygon cut");
  Char_t fname[256];
  TFile* f;
  strcpy( fname, "/home/dave/Desktop/" );
  if( !file ){
    strcat( fname, name );
    strcat( fname, ".root" );
    f = new TFile( fname, "RECREATE" );
  }
  else{
    strcat( fname, file );
    f = new TFile( fname, "UPDATE" );
  }
  cut->Write();
  f->Close();
  delete f;
  printf( "2D cut saved to file %s\n", fname );
}
