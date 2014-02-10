//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of Microscope spectra
//

CheckPID(){
  Char_t* hname[] = {
    "PID_Nhits", "PID_Hits", "PID_Hits_v_EnergyOR", "PID_Hits_v_TimeOR",
    "CB_DeltaE0_v_Echarged0",
    "CB_DeltaE1_v_Echarged1",
    "CB_DeltaE2_v_Echarged2",
    "CB_DeltaE3_v_Echarged3",
    "CB_DeltaE4_v_Echarged4",
    "CB_DeltaE5_v_Echarged5",
    "CB_DeltaE6_v_Echarged6",
    "CB_DeltaE7_v_Echarged7",
    "CB_DeltaE8_v_Echarged8",
    "CB_DeltaE9_v_Echarged9",
    "CB_DeltaE10_v_Echarged10",
    "CB_DeltaE11_v_Echarged11",
    "CB_DeltaE12_v_Echarged12",
    "CB_DeltaE13_v_Echarged13",
    "CB_DeltaE14_v_Echarged14",
    "CB_DeltaE15_v_Echarged15",
    "CB_DeltaE16_v_Echarged16",
    "CB_DeltaE17_v_Echarged17",
    "CB_DeltaE18_v_Echarged18",
    "CB_DeltaE19_v_Echarged19",
    "CB_DeltaE20_v_Echarged20",
    "CB_DeltaE21_v_Echarged21",
    "CB_DeltaE22_v_Echarged22",
    "CB_DeltaE23_v_Echarged23",
  };
  Int_t log[] = { 0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  Int_t col[] = { 2,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,};
  Char_t* xname[] = {
    "Number of PID Hits per Event",
    "PID Hits Distribution",
    "PID Pulse-Height vs Hits",
    "PID Time vs Hits",
    "PID-CB 0 #deltaE-E Plot",
    "PID-CB 1 #deltaE-E Plot",
    "PID-CB 2 #deltaE-E Plot",
    "PID-CB 3 #deltaE-E Plot",
    "PID-CB 4 #deltaE-E Plot",
    "PID-CB 5 #deltaE-E Plot",
    "PID-CB 6 #deltaE-E Plot",
    "PID-CB 7 #deltaE-E Plot",
    "PID-CB 8 #deltaE-E Plot",
    "PID-CB 9 #deltaE-E Plot",
    "PID-CB 10 #deltaE-E Plot",
    "PID-CB 11 #deltaE-E Plot",
    "PID-CB 12 #deltaE-E Plot",
    "PID-CB 13 #deltaE-E Plot",
    "PID-CB 14 #deltaE-E Plot",
    "PID-CB 15 #deltaE-E Plot",
    "PID-CB 16 #deltaE-E Plot",
    "PID-CB 17 #deltaE-E Plot",
    "PID-CB 18 #deltaE-E Plot",
    "PID-CB 19 #deltaE-E Plot",
    "PID-CB 20 #deltaE-E Plot",
    "PID-CB 21 #deltaE-E Plot",
    "PID-CB 22 #deltaE-E Plot",
    "PID-CB 23 #deltaE-E Plot",
  };
  TH1F* h1;
  TCanvas* canv;
  //
  canv = new TCanvas("PID-Spectra","PID-Online",240,180,1240,890);
  canv->SetFillStyle(4000);
  canv->Divide(7,4,0.01,0.01);
  for( Int_t i=0; i<2; i++ ){
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
  for( Int_t i=2; i<28; i++ ){
    h2 = (TH2F*)(gROOT->FindObjectAny(hname[i]));
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


void PIDClear(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("PID"));
  TA2Apparatus* a = (TA2Apparatus*)(gAN->GetChild("CB"));
  if( t ) t->ZeroAll();
  else printf("PID detector class not found\n");
  if( a ) a->ZeroAll();
  else printf("CB apparatus class not found\n");
}

void SavePID(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("PID"));
  t->SaveHist();
}
