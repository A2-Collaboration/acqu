//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of TAPS Veto spectra
//

CheckTAPSVeto(){
  Char_t* hname[] = {
    "VetoBaF2_NADChits",
    "VetoBaF2_NTDChits",
    "VetoBaF2_Nhits",
    "VetoBaF2_ADCHits",
    "VetoBaF2_TDCHits",
    "VetoBaF2_Hits",
    "VetoBaF2_TimeOR",
    "VetoBaF2_Hits_v_EnergyOR",
    "VetoBaF2_Hits_v_TimeOR",
  };
  Int_t log[] = { 0,0,0,1,1,1,0,1,0 };
  Int_t col[] = { 2,2,2,4,4,4,3,1,1 };
  Char_t* xname[] = {
    "Number of TAPS-Veto ADC Hits per Event",
    "Number of TAPS-Veto TDC Hits per Event",
    "Number of TAPS-Veto Hits per Event",
    "TAPS-Veto ADC Hits distribution",
    "TAPS-Veto TDC Hits distribution",
    "TAPS-Veto Hits distribution",
    "TAPS-Veto Time OR (ns)",
    "TAPS-Veto Pulse-Height vs Hits",
    "TAPS-Veto Time vs Hits",
  };
  TH1F* h1;
  TCanvas* canv;
  //
  canv = new TCanvas("TAPS-Veto-Spectra","TAPS-Veto-Online",240,180,1240,890);
  canv->SetFillStyle(4000);
  canv->Divide(3,3,0.01,0.01);
  for( Int_t i=0; i<7; i++ ){
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
  for( Int_t i=7; i<9; i++ ){
    h2 = (TH2F*)(gROOT->FindObject(hname[i]));
    if( !h2 ){
      printf("No root histogram %s\n",hname[i]);
      continue;
    }
    canv->cd(i+1);
    if( log[i] ) canv->GetPad(i+1)->SetLogz();
    h2->GetXaxis()->SetTitle(xname[i]);
    h2->Draw("COLZ");
  }
  return;
}


void TAPSVetoClear(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("VetoBaF2"));
  t->ZeroAll();
}

void SaveTAPSVeto(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("VetoBaF2"));
  t->SaveHist();
}
