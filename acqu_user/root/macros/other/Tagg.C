// Tagger Spectra
// JRM Annand 06/02/03
// DK         01/15/03
//                      some improvements...

TCanvas *CTaggerRawTDC = 0;

void TaggRawTDCplot(UInt_t i){
  // Plot 16 Tagger TDCs beginning @ ADC[i]
  char hist_name[10];
  TH1F *hist;
  UInt_t j, k;
  CTaggerRawTDC = new TCanvas("CTaggerRawTDC", "Tagger: RAW TDCs", 
			      200, 60, 700, 900);  
  CTaggerRawTDC->Divide(4, 4);
  for( j = i, k = 1; k<17; k++, ++j){
    CTaggerRawTDC->cd(k);
    sprintf(hist_name, "ADC%d", j);
    hist = (TH1F)gROOT->Get(hist_name);
    hist->Draw();
  }
  CTaggerRawTDC->Update();
  CTaggerRawTDC->Draw();
}



void TaggRawTDCfit(UInt_t i){
  // under construction. don't use!
  // it's only to stydy some root-features!
  char hist_name[10];
  TH1F *hist; 
  TF1  *fit;
  TLine *line;
  Double_t mu;
  CTaggerRawTDC = new TCanvas("CTaggerRawTDC", "Tagger: RAW TDCs", 
			      200, 60, 700, 900);  
  sprintf(hist_name, "ADC%d", i);
  hist = (TH1F)gROOT->Get(hist_name);
  hist->Fit("gaus", "Q"); // "q"
  fit = hist->GetFunction("gaus");
  mu = (fit->GetParameter(1)); // µ of the fit
  gPad->Update();
  line = new TLine(mu, 0, mu,  gPad->GetUymax());
  line->SetLineColor(kRed);
  hist->Draw();
  line->Draw();  
}

void TaggRawTDC(){
  // Create a nice menubar to check all the Tagger TDCs
  // ADC-channels of Tagger raw TDCs: 1200 - 1551.
  TControlBar* bar = new TControlBar("vertical", "Tagger - Raw TDCs");
  char bar_title[50];
  char bar_action[50];
  char c = 'A';
  for(UInt_t i = 1200; i<1552; i+=16){
    sprintf(bar_title,  "&%c: TDCs %d - %d", c++, i, i+15);
    sprintf(bar_action, "TaggRawTDCplot(%d)", i);
    bar->AddButton(bar_title, bar_action, " ");
  }
  bar->Show();
}

// -----------------------------------------------------------------------------

void TaggClear(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("Ladder"));
  t->ZeroAll();
}

void TaggSave(){
  TA2Detector* t = (TA2Detector*)(gAN->GetGrandChild("Ladder"));
  t->SaveHist();
}
// -----------------------------------------------------------------------------

Tagg(){
  TA2Detector* det = (TA2Detector*)(gAN->GetGrandChild("Ladder"));
  TCanvas* canv = new TCanvas("Tagg", "Tagger stuff", 200, 60, 700, 900);
  TObject *the_hist;
  int i, j;
  char* hnames[] = {"Ladder_TimeOR", "",
		    "Ladder_Hits", "", 
		    "Ladder_Nhits", "", 
		    "Ladder_Hits_v_TimeOR", "COLZ"};
  char *or_pos = 0;
  canv->SetFillStyle(4000);
  canv->Divide(2,2,0.01,0.01);
  for( j=0, i=1; i<5; j+=2, i++){
    canv->cd(i);
    //  printf("processing %s\n", hnames[j]);
    the_hist = gROOT->FindObject(hnames[j]);
    if (the_hist == 0){
      or_pos = strstr(hnames[j], "OR");
      if (or_pos > 0){
	// printf("OR %d  %s", or_pos, hnames[j]);
	*or_pos = '\0';
	or_pos = 0;
	//printf(" --> %s\n", hnames[j]);
	the_hist = gROOT->FindObject(hnames[j]);
      }
    }
    if(the_hist) the_hist->Draw(hnames[j+1]);
 }
}

