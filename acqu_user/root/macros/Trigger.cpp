
//-----------------------------------------------------------------------------

void DiffEsum(Char_t* Meson, Int_t nRebin, Double_t Low, Double_t High)
{
  Char_t Buff[256];
  TH1F* Temp;
  TH1F* ADC;
  TH1F* DADC;
  TF1* f;
  Int_t nBins;
  Int_t nLow, nHigh;
  Double_t y1, y2;
  Double_t peak;
  Double_t par[3];

  ADC =  (TH1F*)gROOT->FindObject("ESum");
  DADC = (TH1F*)gROOT->FindObject("DESum");
  if(ADC)  ADC->Delete();
  if(DADC) DADC->Delete();

  sprintf(Buff, "%s_ESum", Meson);
  Temp = (TH1F*)gROOT->FindObject(Buff);
  ADC = (TH1F*)Temp->Clone("ESum");
  if(nRebin) ADC->Rebin(nRebin);
  DADC = (TH1F*)ADC->Clone("DESum");

  nBins = ADC->GetXaxis()->GetNbins();
  for(Int_t t=2; t<nBins+1; t++)
  {
    y1 = ADC->GetBinContent(t-1);
    y2 = ADC->GetBinContent(t);
    DADC->SetBinContent(t, (y2-y1));
  }

  nLow = Low/DADC->GetBinWidth(1);
  for(Int_t t=0; t<nLow+1; t++)
    DADC->SetBinContent(t, 0);
  nHigh = High/DADC->GetBinWidth(1);
  //for(Int_t t=nHigh+1; t<999; t++)
    //DADC->SetBinContent(t, 0);

  peak = (DADC->GetMaximumBin())*DADC->GetBinWidth(1);
  f = new TF1("f", "gaus", peak-150.0, High);//peak+100.0);
  f->SetLineWidth(1);
  f->SetLineColor(2);
  DADC->Fit(f,"RQ+");
  f->GetParameters(&par[0]);

  if(par[1]==0.0) par[1]=999.999999;
  printf("%f +- %f", par[1], par[2]);
  if(par[2] < 0.5) printf("!!!");
  printf("\n");
  DADC->Draw();
}

//-----------------------------------------------------------------------------

void DiffThresNaI(Char_t* Trigger, Int_t nChan, Int_t nRebin, Double_t Low, Double_t High)
{
  Char_t Buff[256];
  TH1F* Temp;
  TH1F* ADC;
  TH1F* DADC;
  TF1* f;
  Int_t nBins;
  Int_t nLow, nHigh;
  Double_t y1, y2;
  Double_t peak;
  Double_t par[3];

  ADC =  (TH1F*)gROOT->FindObject("NaIThres");
  DADC = (TH1F*)gROOT->FindObject("DNaIThres");
  if(ADC)  ADC->Delete();
  if(DADC) DADC->Delete();
  
  sprintf(Buff, "%s_NaIThres%d", Trigger, nChan);
  Temp = (TH1F*)gROOT->FindObject(Buff);
  ADC = (TH1F*)Temp->Clone("NaIThres");
  if(nRebin) ADC->Rebin(nRebin);
  DADC = (TH1F*)ADC->Clone("DNaIThres");

  nBins = ADC->GetXaxis()->GetNbins();
  for(Int_t t=2; t<nBins+1; t++)
  {
    y1 = ADC->GetBinContent(t-1);
    y2 = ADC->GetBinContent(t);
    DADC->SetBinContent(t, (y2-y1));
  }

  nLow = Low/DADC->GetBinWidth(1);
  for(Int_t t=0; t<nLow+1; t++)
    DADC->SetBinContent(t, 0);
  nHigh = High/DADC->GetBinWidth(1);
  for(Int_t t=nHigh+1; t<999; t++)
    DADC->SetBinContent(t, 0);

  peak = (DADC->GetMaximumBin())*DADC->GetBinWidth(1);
  f = new TF1("f", "gaus", peak-2.0, peak+2.0);
  f->SetLineWidth(1);
  f->SetLineColor(2);
  DADC->Fit(f,"RQ+");
  f->GetParameters(&par[0]);

  if(par[1]==0.0) par[1] = 9999.999;
  printf("%3d  %8.3f  %8.3f\n", nChan, par[1], par[2]);
  DADC->Draw();
}

//-----------------------------------------------------------------------------

void DiffThresBaF2(Char_t* Trigger, Int_t nChan, Int_t nRebin, Double_t Low, Double_t High)
{
  Char_t Buff[256];
  TH1F* Temp;
  TH1F* ADC;
  TH1F* DADC;
  TF1* f;
  Int_t nBins;
  Int_t nLow, nHigh;
  Double_t y1, y2;
  Double_t peak;
  Double_t par[3];

  ADC =  (TH1F*)gROOT->FindObject("BaF2Thres");
  DADC = (TH1F*)gROOT->FindObject("DBaF2Thres");
  if(ADC)  ADC->Delete();
  if(DADC) DADC->Delete();

  sprintf(Buff, "%s_BaF2Thres%d", Trigger, nChan);
  Temp = (TH1F*)gROOT->FindObject(Buff);
  ADC = (TH1F*)Temp->Clone("BaF2Thres");
  if(nRebin) ADC->Rebin(nRebin);
  DADC = (TH1F*)ADC->Clone("DBaF2Thres");

  nBins = ADC->GetXaxis()->GetNbins();
  for(Int_t t=2; t<nBins+1; t++)
  {
    y1 = ADC->GetBinContent(t-1);
    y2 = ADC->GetBinContent(t);
    DADC->SetBinContent(t, (y2-y1));
  }

  nLow = Low/DADC->GetBinWidth(1);
  for(Int_t t=0; t<nLow+1; t++)
    DADC->SetBinContent(t, 0);
  nHigh = High/DADC->GetBinWidth(1);
  for(Int_t t=nHigh+1; t<999; t++)
    DADC->SetBinContent(t, 0);

  peak = (DADC->GetMaximumBin())*DADC->GetBinWidth(1);
  f = new TF1("f", "gaus", peak-2.0, peak+2.0);
  f->SetLineWidth(1);
  f->SetLineColor(2);
  DADC->Fit(f,"RQ+");
  f->GetParameters(&par[0]);

  if(par[1]==0.0) par[1] = 9999.999;
  printf("%3d  %8.3f  %8.3f\n", nChan, par[1], par[2]);
  DADC->Draw();
}

//-----------------------------------------------------------------------------

void CreateSumModel(Char_t* Filename1, Char_t* Filename2, Char_t* Meson, Int_t Low)
{
  Char_t Buffer[256];
  TFile* File1 = new TFile(Filename1);
  TFile* File2 = new TFile(Filename2);
  FILE* Output = fopen("Trigger_SumModel.txt","w");
  
  sprintf(Buffer,"%s_ESum", Meson);
  
  TH1F* Esum1 = (TH1F*)File1->Get(Buffer);
  TH1F* Esum2 = (TH1F*)File2->Get(Buffer);
  
  TH1F* Model = (TH1F*)Esum1->Clone("SumModel");
  Model->Divide(Esum2);
  for(Int_t t=0; t<(Low/2); t++)
    Model->SetBinContent(t+1, 0.0);
  Model->Scale(1.0/Model->GetBinContent(Model->GetMaximumBin()));

  for(Int_t t=0; t<750; t++)
  {    
    printf("%d  %f\n", t*2,   Model->GetBinContent(t+1));
    printf("%d  %f\n", t*2+1, Model->GetBinContent(t+1));
    fprintf(Output, "%d  %f\n", t*2,   Model->GetBinContent(t+1));
    fprintf(Output, "%d  %f\n", t*2+1, Model->GetBinContent(t+1));
  }
  for(Int_t t=1500; t<2001; t++)
  {
    printf("%d  %f\n", t, 1.0);
    fprintf(Output, "%d  %f\n", t, 1.0);
  }
  fclose(Output);

  Model->SetMinimum(0.00);
  Model->SetMaximum(1.05);
  Model->Draw();
}

//-----------------------------------------------------------------------------

void DiffEOrNaI(Int_t nRebin, Double_t Low, Double_t High)
{
  Char_t Buff[256];
  TH1F* Temp;
  TH1F* ADC;
  TH1F* DADC;
  TF1* f;
  Int_t nBins;
  Int_t nLow, nHigh;
  Double_t y1, y2;
  Double_t peak;
  Double_t par[3];

  ADC =  (TH1F*)gROOT->FindObject("ESum");
  DADC = (TH1F*)gROOT->FindObject("DESum");
  if(ADC)  ADC->Delete();
  if(DADC) DADC->Delete();

  sprintf(Buff, "NaI_EnergyOR");
  Temp = (TH1F*)gROOT->FindObject(Buff);
  ADC = (TH1F*)Temp->Clone("ESum");
  if(nRebin) ADC->Rebin(nRebin);
  DADC = (TH1F*)ADC->Clone("DESum");

  nBins = ADC->GetXaxis()->GetNbins();
  for(Int_t t=2; t<nBins+1; t++)
  {
    y1 = ADC->GetBinContent(t-1);
    y2 = ADC->GetBinContent(t);
    DADC->SetBinContent(t, (y2-y1));
  }

  nLow = Low/DADC->GetBinWidth(1);
  for(Int_t t=0; t<nLow+1; t++)
    DADC->SetBinContent(t, 0);
  nHigh = High/DADC->GetBinWidth(1);
  for(Int_t t=nHigh+1; t<999; t++)
    DADC->SetBinContent(t, 0);

  peak = (DADC->GetMaximumBin())*DADC->GetBinWidth(1);
  f = new TF1("f", "gaus", peak-150.0, peak+100.0);
  f->SetLineWidth(1);
  f->SetLineColor(2);
  DADC->Fit(f,"RQ+");
  f->GetParameters(&par[0]);

  if(par[1]==0.0) par[1]=999.999999;
  printf("%f +- %f", par[1], par[2]);
  if(par[2] < 0.5) printf("!!!");
  printf("\n");
  DADC->Draw();
}

//-----------------------------------------------------------------------------

void DiffEOrBaF2(Int_t nRebin, Double_t Low, Double_t High)
{
  Char_t Buff[256];
  TH1F* Temp;
  TH1F* ADC;
  TH1F* DADC;
  TF1* f;
  Int_t nBins;
  Int_t nLow, nHigh;
  Double_t y1, y2;
  Double_t peak;
  Double_t par[3];

  ADC =  (TH1F*)gROOT->FindObject("ESum");
  DADC = (TH1F*)gROOT->FindObject("DESum");
  if(ADC)  ADC->Delete();
  if(DADC) DADC->Delete();

  sprintf(Buff, "BaF2_EnergyOR");
  Temp = (TH1F*)gROOT->FindObject(Buff);
  ADC = (TH1F*)Temp->Clone("ESum");
  if(nRebin) ADC->Rebin(nRebin);
  DADC = (TH1F*)ADC->Clone("DESum");

  nBins = ADC->GetXaxis()->GetNbins();
  for(Int_t t=2; t<nBins+1; t++)
  {
    y1 = ADC->GetBinContent(t-1);
    y2 = ADC->GetBinContent(t);
    DADC->SetBinContent(t, (y2-y1));
  }

  nLow = Low/DADC->GetBinWidth(1);
  for(Int_t t=0; t<nLow+1; t++)
    DADC->SetBinContent(t, 0);
  nHigh = High/DADC->GetBinWidth(1);
  for(Int_t t=nHigh+1; t<999; t++)
    DADC->SetBinContent(t, 0);

  peak = (DADC->GetMaximumBin())*DADC->GetBinWidth(1);
  f = new TF1("f", "gaus", peak-150.0, peak+100.0);
  f->SetLineWidth(1);
  f->SetLineColor(2);
  DADC->Fit(f,"RQ+");
  f->GetParameters(&par[0]);

  if(par[1]==0.0) par[1]=999.999999;
  printf("%f +- %f", par[1], par[2]);
  if(par[2] < 0.5) printf("!!!");
  printf("\n");
  DADC->Draw();
}

//-----------------------------------------------------------------------------
