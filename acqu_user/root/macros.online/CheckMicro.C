//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of Microscope spectra
//

CheckMicro(){
  Char_t* hname[] = {
    "Micro_Nhits",
    "Micro_NhitsPrompt",
    "Micro_NhitsRand",  
    "Micro_Hits",
    "Micro_ScalerCurr",
    "Micro_ScalerAcc",
    "Micro_EnergyOR",
    "Micro_TimeOR",
    "Micro_NMuHits", 
    "Micro_MuHits",
    "Micro_Hits_v_EnergyOR",
    "Micro_Hits_v_TimeOR",
  };
  Int_t log1D[] = { 0,0,0,0,1,1,0,0,0,0 };
  Int_t col1D[] = { 2,2,2,4,5,5,3,3,2,4 };
  Char_t* xname[] = {
    "Number of Microscope Hits per Event",
    "Number of Prompt Microscope Hits per Event",
    "Number of Random Microscope Hits per Event",
    "Microscope Hits distribution",
    "Microscope Current Scaler Read",
    "Microscope Accumulated Scalers",
    "Microscope Pulse-Height (MeV)",
    "Microscope Time (ns)",
    "Number of Microscope Overlap Hits per Event",
    "Microscope Overlap Hits Distribution",
    "Microscope Pulse-Height vs Hits",
    "Microscope Time vs Hits",
  };
  TH1F* h1;
  TCanvas* canv;
  //
  canv = new TCanvas("Micro-Spectra","Micro-Online",240,180,1240,890);
  canv->SetFillStyle(4000);
  canv->Divide(3,4,0.01,0.01);
  for( Int_t i=0; i<10; i++ ){
      h1 = (TH1F*)(gROOT->FindObject(hname[i]));
      if( !h1 ){
	printf("No root histogram %s\n",hname[i]);
	continue;
      }
      h1->SetLineColor( 1 );
      h1->SetFillColor( col1D[i] );
      canv->cd(i+1);
      if( log1D[i] ) canv->GetPad(i+1)->SetLogy();
      h1->GetXaxis()->SetTitle(xname[i]);
      h1->Draw();
  }
  //  return;
  TH2F* h2;
  for( Int_t i=10; i<12; i++ ){
    h2 = (TH2F*)(gROOT->FindObject(hname[i]));
    if( !h2 ){
      printf("No root histogram %s\n",hname[i]);
      continue;
    }
    canv->cd(i+1);
    h2->GetXaxis()->SetTitle(xname[i]);
    h2->Draw("COLZ");
  }
  return;
}


void MicroClear(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("Micro"));
  t->ZeroAll();
}

void SaveMicro(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("Micro"));
  t->SaveHist();
}
