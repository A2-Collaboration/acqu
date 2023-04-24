//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of CATS detector spectra
//

void CATSAnnulusClear(){
	TA2Detector* d = (TA2Detector*)(gAN->GetGrandChild("CATS_Annulus"));
	if( d ) d->ZeroAll();
	else printf("CATS Annulus detector class not found\n");
}

CheckCATSAnnulus(TCanvas* canv){
	if(canv==NULL) {
		CATSAnnulusClear();
		return;
	}
 
  Char_t* hname[] = {
    "CATS_Annulus_Energy0",
    "CATS_Annulus_Energy1",
    "CATS_Annulus_Energy2",
    "CATS_Annulus_Energy3",
    "CATS_Annulus_Energy4",
    "CATS_Annulus_Energy5",
  };
  Int_t log[] = { 1,1,1,1,1,1 };
  Int_t col[] = { 2,2,2,2,2,2 };
  Char_t* xname[] = {
    "Annulus 0",
    "Annulus 1",
    "Annulus 2",
    "Annulus 3",
    "Annulus 4",
    "Annulus 5",
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

