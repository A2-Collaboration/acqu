// dk - 28th Mai 04
// Short macro to check/clear some TAPS BaF2 spectra
// Feel free to share your share your wisdom and improve it...
// dk - 22th june 04
// More specrtra...

void CheckTAPSBaF2(){
  Char_t* hname[] = {
    "BaF2_ADCHits",
    "BaF2_TDCHits",
    "BaF2_Hits",
    "BaF2_ClNhits",
    "BaF2_ClNhitsOR",
    "BaF2_Nhits",
    "BaF2_Etot",
    "BaF2_TimeOR",
    "BaF2_ClNhitsOR_v_ClEnergyOR",
    "BaF2_Hits_v_EnergyOR",
    "BaF2_Hits_v_TimeOR",
    "BaF2_ClPhi_v_ClTheta",
  };
  Char_t* xname[] = {
    "TAPS Distribution of ADC Hits",
    "TAPS Distribution of TDC Hits",
    "TAPS Distribution of Hits",
    "TAPS Number Hits per Event",
    "TAPS Number of Clusters per Event",
    "TAPS Number of Hits per Cluster",
    "TAPS Total Energy in Calorimeter",
    "TAPS Time OR",
    "TAPS Hits per Cluster vs Cluster Energy",
    "TAPS Energy vs Hits",
    "TAPS Time vs Hits",
    "TAPS #theta - #phi distribution of Cluster Hits",
  };
  Int_t log[] = { 0,0,0,0,0,0,1,0,1,1,0,0 };
  Int_t col[] = { 4,4,4,2,2,2,3,3,1,1,1,1 };
  Char_t* opt2d[] = { "cont0", "colz", "colz", "colz" };
  TH1F* h1;
  TH2F* h2;
  TCanvas* canv = new TCanvas("TAPS-BaF2-Spectra","TAPS-BaF2-Online",240,180,1240,890);
  canv->Divide(3,4);
  for( Int_t i=0; i<8; i++ ){
      h1 = (TH1F*)(gROOT->FindObjectAny(hname[i]));
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
  Int_t j=0;
  for( Int_t i=8; i<12; i++ ){
    h2 = (TH2F*)(gROOT->FindObjectAny(hname[i]));
    if( !h2 ){
      printf("No root histogram %s\n",hname[i]);
      continue;
    }
    canv->cd(i+1);
    h2->GetXaxis()->SetTitle(xname[i]);
    if( log[i] ) canv->GetPad(i+1)->SetLogz();
    //    canv->GetPad(i+1)->SetStyle(col[i]);
    h2->Draw(opt2d[j]);
    j++;
  }
  return;
}

TAPSBaF2Clear(){
  // Does not work properly. (Reset()!)
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("BaF2"));
  t->ZeroAll();
}
