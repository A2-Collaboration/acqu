//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of TAPS Veto spectra
//

CheckLiveTimes(){
  Char_t* hname[] = {
    "LiveTotal",
    "LiveTrig",
    "LiveTaggA",
    "LiveTaggB",
    "LiveBeamPol",
    "LiveADC1A",
    "LiveADC1B",
    "LiveADC2A",
    "LiveADC2B",
    "LiveTDCA",
    "LiveTDCB",
    "LiveMWPCA",
    "LiveMWPCB",
    "LiveTrigTAPS",
    "LiveBafA",
    "LiveBafB",
    "LiveBafC",
    "LiveBafD",
    "LiveBafE",
    "LiveBafF",
    "LiveVetoAB",
    "LiveVetoCD",
    "LiveVetoEF",
    "LivePbWO4",
 };
  Char_t* xname[] = {
    "Livetime Total",
    "Livetime Trig",
    "Livetime TaggA",
    "Livetime TaggB",
    "Livetime Beam-Pol-Mon",
    "Livetime ADC1A",
    "Livetime ADC1B",
    "Livetime ADC2A",
    "Livetime ADC2B",
    "Livetime TDCA",
    "Livetime TDCB",
    "Livetime MWPCA",
    "Livetime MWPCB",
    "Livetime TrigTAPS",
    "Livetime BafA",
    "Livetime BafB",
    "Livetime BafC",
    "Livetime BafD",
    "Livetime BafE",
    "Livetime BafF",
    "Livetime VetoAB",
    "Livetime VetoCD",
    "Livetime VetoEF",
    "Livetime PbWO4",
};
  Int_t log[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  Int_t col[] = { 4,4,5,5,2,2,2,2,2,2,5,5,3,3,3,4,4,4,4,4,4,2,2,2,5 };
 //
  TH1F* h1;
  TCanvas* canv;
  canv = new TCanvas("Livetimes-Spectra","LiveTime Spectra",240,180,1240,890);
  canv->SetFillStyle(4000);
  canv->Divide(4,6,0.01,0.01);
  for( Int_t i=0; i<24; i++ ){
      h1 = (TH1F*)(gROOT->FindObject(hname[i]));
      if( !h1 ){ printf("No root histogram %s\n",hname[i]); continue; }
      h1->SetLineColor(1);
      h1->SetFillColor( col[i] );
      canv->cd(i+1);
      if( log[i]==1 ) canv->GetPad(i+1)->SetLogy();
      h1->GetXaxis()->SetTitle(xname[i]);
      h1->Draw();
  }
}


void LiveTimesClear(){
  printf("Not yet implemented\n");
}

void SaveLiveTimes(){
  printf("Not yet implemented\n");
}
