//--Author	JRM Annand   13th Jan 2007
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// Check plots of CATS detector spectra
//

void CATSCoreClear(){
	TA2Detector* d = (TA2Detector*)(gAN->GetGrandChild("CATS_Core"));
	if( d ) d->ZeroAll();
	else printf("CATS Core detector class not found\n");
}

CheckCATSCore(TCanvas* canv){
	if(canv==NULL) {
		CATSCoreClear();
		return;
	}
 
  Char_t* hname[] = {
    "CATS_Core_Nhits",
    "CATS_Core_Etot",
    "CATS_Core_Energy0",
    "CATS_Core_Energy1",
    "CATS_Core_Energy2",
    "CATS_Core_Energy3",
    "CATS_Core_Energy4",
    "CATS_Core_Energy5",
    "CATS_Core_Energy6",
  };
  Int_t log[] = { 0,0,0,0,0,0,0,0,0 };
  Int_t col[] = { 2,3,4,4,4,4,4,4,4 };
  Char_t* xname[] = {
    "Core N",
    "Core T",
    "Core 0",
    "Core 1",
    "Core 2",
    "Core 3",
    "Core 4",
    "Core 5",
    "Core 6",
  };
  TH1F* h1;
  canv->SetFillStyle(4000);
  canv->Divide(3,3,0.01,0.01);
  for( Int_t i=0; i<9; i++ ){
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

