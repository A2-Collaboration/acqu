//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of Tagger spectra
//

CheckTagger(){
  Char_t* hname[] = {
    "FPD_Nhits", "FPD_NhitsPrompt", "FPD_NhitsRand",  
    "FPD_Hits", "FPD_ScalerCurr", "FPD_ScalerAcc",
    "FPD_EnergyOR","FPD_TimeOR","TAGG_PhotonEnergy",
    "FPD_Hits_v_EnergyOR",
    "FPD_Hits_v_TimeOR",
    "TAGG_Micro_Hits_v_FPD_Hits",
  };
  Int_t log[] = { 0,0,0,0,1,1,0,0,0 };
  Int_t col[] = { 2,2,2,4,5,5,3,3,3 };
  Char_t* xname[] = {
    "Number of Tagger Hits per Event",
    "Number of Prompt Tagger Hits per Event",
    "Number of Random Tagger Hits per Event",
    "Tagger Hits distribution",
    "Tagger Current Scaler Read",
    "Tagger Accumulated Scalers",
    "Tagger Pulse-Height (MeV)",
    "Tagger Time (ns)",
    "Tagger Photon Energy",
    "Tagger Pulse-Height vs Hits",
    "Tagger Time vs Hits",
    "Microscope Hits vs Tagger Hits",
  };
  TH1F* h1;
  TCanvas* canv;
  //
  canv = new TCanvas("Tagger-Spectra","Tagger-Online",240,180,1240,890);
  canv->SetFillStyle(4000);
  canv->Divide(3,4,0.01,0.01);
  for( Int_t i=0; i<9; i++ ){
      h1 = (TH1F*)(gROOT->FindObject(hname[i]));
      if( !h1 ){
	printf("No root histogram %s\n",hname[i]);
	continue;
      }
      h1->SetLineColor( 1 );
      h1->SetFillColor( col[i] );
      canv->cd(i+1);
      if( log[i] ) canv->GetPad(i+1)->SetLogy();
      h1->GetXaxis()->SetTitle(xname[i]);
      h1->Draw();
  }
  //  return;
  TH2F* h2;
  for( Int_t i=9; i<12; i++ ){
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


void TaggerClear(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("FPD"));
  t->ZeroAll();
}

void SaveTagger(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("FPD"));
  t->SaveHist();
}
