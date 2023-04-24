//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of CATS detector spectra
//

void CATSShieldClear(){
	TA2Detector* d = (TA2Detector*)(gAN->GetGrandChild("CATS_Shield"));
	if( d ) d->ZeroAll();
	else printf("CATS Shield detector class not found\n");
}

CheckCATSShield(TCanvas* canv){
	if(canv==NULL) {
		CATSShieldClear();
		return;
	}
 
  Char_t* hname[] = {
    "CATS_Shield_Nhits",
    "CATS_Shield_Energy0",
    "CATS_Shield_Energy1",
    "CATS_Shield_Energy2",
    "CATS_Shield_Energy3",
    "CATS_Shield_Energy4",
  };
  Int_t log[] = { 1,1,1,1,1,1 };
  Int_t col[] = { 2,4,4,4,4,4 };
  Char_t* xname[] = {
    "Shield N",
    "Shield 0",
    "Shield 1",
    "Shield 2",
    "Shield 3",
    "Shield 4",
  };
  TH1F* h1;
  canv->SetFillStyle(4000);
  canv->Divide(3,2,0.01,0.01);
  for( Int_t i=0; i<6; i++ ){
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
}

