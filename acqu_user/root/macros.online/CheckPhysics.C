void CheckPhysics(){
  Char_t* hname[] = {
    "1",  "EtaM_Nphoton",   "log",  "",
    "2",  "EtaM_Nproton",   "log",  "",
    //    "3",  "EtaM_Nneutron",  "log",  "",
    "3",  "EtaM_Npi0",      "log",  "",
    //    "5",  "EtaM_Npiplus",   "log",  "",
    "4",  "EtaM_Neta",      "log",  "",
    "5",  "EtaM_M2g",      "log",  "",
    "6",  "EtaM_M6g",      "lin",  "",
    "7",  "EtaM_EmEtaP",   "lin",  "",
    "7",  "EtaM_EmEtaR",   "lin",  "Same",
    "8",  "EtaM_EmEtapP",   "lin",  "",
    "8",  "EtaM_EmEtapR",   "lin",  "Same",
    //    "8",  "EtaM_EmComptonpP",   "lin",  "",
    //    "8",  "EtaM_EmComptonpR",   "lin",  "Same",
    "9",  "EtaM_M2gCBTAPS","lin",  "",
    "10", "EtaM_M2g_v_CB_NaI_ClNhitsOR",        "log",  "CONTZ",
    "11", "EtaM_M2g_v_TAGG_FPD_HitsPrompt",     "log",  "COLZ",
    "12", "EtaM_M6g_v_TAGG_FPD_HitsPrompt",       "log",  "COLZ",
  };
  Int_t col[] = { 2,2,2,2, 3,3,3,5,3,5,3,2,2,2,1,1 };
  Int_t xxx;
  Char_t* xname[] = {
    "Number Reconstructed #gamma per event",
    "Number Reconstructed p per event",
    //    "Number Reconstructed n per event",
    "Number Reconstructed #pi^{0} per event",
    //    "Number Reconstructed #pi^{+} per event",
    "Number Reconstructed #eta per event",
    "2-#gamma Invariant Mass",
    "6-#gamma Invariant Mass",
    "Missing Mass if QF #eta (Gr-prompt, Yell-rand) ",
    "",
    "Missing Energy if QF #eta-p (Gr-prompt, Yell-rand)",
    "",
    "2-#gamma Invariant Mass (1-#gamma CB, 1-#gamma TAPS)",
    "2-#gamma Invariant Mass vs CB Cluster Size",
    "2-#gamma Invariant Mass vs Prompt Tagger Hits",
    "6-#gamma Invariant Mass vs Prompt Tagger Hits",
  };
  Char_t* cloneName[] = {"rand1", "rand2", "rand3", "rand4"};
  TH1F* h1r;
  TH1F** ph1r[4];
  for( Int_t i=0; i<2; i++ ){
    h1r = (TH1F*)gROOT->FindObject(cloneName[i]);
    if( !h1r ) h1r = new TH1F("cloneName[i]","Photo-Hist",1000,0,1000);
    ph1r[i] = h1r;
  }
  //
  TH1F* h1;
  TH2F* h2;
  //  TH1F* h3 = new TH1F("h3x","Photo-Hist",10000,-5000,5000);
  //  TH1F* h4 = new TH1F("h4x","Photo-Hist",10000,-5000,5000);
  //  TH1F* h5 = new TH1F("h5x","Photo-Hist",10000,-5000,5000);
  TCanvas* canv = new TCanvas("Physics-Spectra","Physics-Online",240,180,1240,890);
  canv->Divide(4,3);
  Int_t i = 0;
  Int_t k = 0;
  Int_t j;
  Int_t pad_no;
  for( j=0; j<44; j+=4){
    sscanf(hname[j], "%d", &pad_no);           // pad_no is 1st arg.
    canv->cd(pad_no);
    //    printf("%d %d %d %s %s %s %s\n",pad_no,j,i,hname[j],hname[j+1],hname[j+2],hname[j+3]);
    h1 = (TH1F*)gROOT->FindObjectAny(hname[j+1]); // name is 2nd argument
    if( !h1 ){
      printf("No object named '%s' was found.\n Check your macro!\n", hname[j+1]);
    }
    else{
      if(strcmp(hname[j+3], "Same") == 0){
	h1r = ph1r[k];
	h1->Copy(*h1r);
	*h1r = *h1r * 0.5;
	h1 = h1r;
	k++;
      }
      if(strcmp(hname[j+2], "log") == 0)         // lin/log 3rd argument
	canv->GetPad(pad_no)->SetLogy();
      h1->GetXaxis()->SetTitle(xname[i]);
      h1->SetLineColor( 1 );
      h1->SetFillColor( col[i] );
      h1->Draw(hname[j+3]);               // parameter 4th arg.
      gPad->RedrawAxis();
    }
    i++;
  }
  for( j=44; j<56; j+=4){
    sscanf(hname[j], "%d", &pad_no);           // pad_no is 1st arg.
    canv->cd(pad_no);
    //    printf("%d %d %d\n",pad_no,j,i);
    h2 = (TH2F*)gROOT->FindObjectAny(hname[j+1]); // name is 2nd argument
    if( !h2 ){
      printf("No object named '%s' was found.\n Check your macro!\n", hname[j+1]);
    }
    else{
      if(strcmp(hname[j+2], "log") == 0)         // lin/log 3rd argument
	canv->GetPad(pad_no)->SetLogz();
      h2->GetXaxis()->SetTitle(xname[i]);
      h2->Draw(hname[j+3]);               // parameter 4th arg.
      gPad->RedrawAxis();
    }
    i++;
  }
  return;
}

PhysicsClear(){
  // Does not work properly. (Reset()!)
  //  TA2Apparatus* t = (TA2Apparatus*)gROOT->FindObject("EtaM");
  //  if(t)  t->ZeroAll();
  //  else printf("EtaM Physics not found\n");
  //  Scalers364to491->Reset();
  printf("Clear Physics spectra not implemented...try Clear Everything\n");
}
