
//-----------------------------------------------------------------------------

void LadderTimes()
{
  TH1F* Time;
  Double_t par[4];
  Double_t Max;
  Char_t Buff[256];

  for(Int_t ch=0; ch<352; ch++)
  {
    sprintf(Buff, "Ladder_Time%d", ch);
    Time = (TH1F*)gROOT->FindObject(Buff);
    Max = Time->GetMaximumBin()*Time->GetBinWidth(Time->GetMaximumBin()) - 1000.0;
   
    TF1* gauss = new TF1("gauss", "gaus", Max-8, Max+8);
    Time->Fit(gauss, "RQ+");
    gauss->GetParameters(&par[0]);
    
    printf("%5.2f\n", par[1]/0.117710);
  }
}

//-----------------------------------------------------------------------------

void BaF2Times()
{
  TH1F* Time;
  Char_t Buff[256];

  for(Int_t ch=0; ch<438; ch++)
  {
    sprintf(Buff, "BaF2_Time%d", ch);
    Time = (TH1F*)gROOT->FindObject(Buff);
    printf("%5.1f\n", Time->GetMaximumBin()*Time->GetBinWidth(Time->GetMaximumBin()) - 300.0);
  }
}

//-----------------------------------------------------------------------------

void NaITimes()
{
  TH1F* Time;
  Char_t Buff[256];

  for(Int_t ch=0; ch<720; ch++)
  {
    sprintf(Buff, "NaI_Time%d", ch);
    Time = (TH1F*)gROOT->FindObject(Buff);
    printf("%5.1f\n", Time->GetMaximumBin()*Time->GetBinWidth(Time->GetMaximumBin()) - 200.0);
  }
}

//-----------------------------------------------------------------------------

void PIDTimes()
{
  TH1F* Time;
  Char_t Buff[256];

  for(Int_t ch=0; ch<24; ch++)
  {
    sprintf(Buff, "PID_Time%d", ch);
    Time = (TH1F*)gROOT->FindObject(Buff);
    printf("%5.1f\n", Time->GetMaximumBin()*Time->GetBinWidth(Time->GetMaximumBin()) - 200.0);
  }
}

void PIDPedes()
{
  TH1F* Energy;
  Char_t Buff[256];

  for(Int_t ch=0; ch<24; ch++)
  {
    sprintf(Buff, "PID_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    printf("%4.0f\n", Energy->GetMaximumBin()*Energy->GetBinWidth(Energy->GetMaximumBin()));
  }
}

void PIDEnergy()
{
  TH1F* Energy;
  Char_t Buff[256];
  Double_t MaxVal;
  Int_t MaxBin;

  for(Int_t ch=0; ch<24; ch++)
  {
    sprintf(Buff, "PID_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    MaxVal = 0.0;
    MaxBin = 0;
    for(Int_t bn=60; bn<4096; bn++)
      if(Energy->GetBinContent(bn) > MaxVal)
      {
        MaxVal = Energy->GetBinContent(bn);
        MaxBin = bn;
      }
    printf("%4d\n", MaxBin);
  }
}

void PIDThres()
{
  TH1F* Energy;
  Char_t Buff[256];
  Double_t MinVal;
  Int_t MinBin;

  for(Int_t ch=0; ch<24; ch++)
  {
    sprintf(Buff, "PID_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    MinVal = 1e38;
    MinBin = 0;
    for(Int_t bn=1; bn<200; bn++)
      if(Energy->GetBinContent(bn) < MinVal)
      {
        MinVal = Energy->GetBinContent(bn);
        MinBin = bn;
      }
    printf("%4d\n", MinBin);
  }
}

//-----------------------------------------------------------------------------

void VetoTimes()
{
  TH1F* Time;
  Char_t Buff[256];

  for(Int_t ch=0; ch<384; ch++)
  {
    sprintf(Buff, "Veto_Time%d", ch);
    Time = (TH1F*)gROOT->FindObject(Buff);
    printf("%6.1f\n", (Time->GetMaximumBin()*Time->GetBinWidth(Time->GetMaximumBin()) - 500.0)/0.05000);
  }
}

void VetoPedes()
{
  TH1F* Energy;
  Char_t Buff[256];

  for(Int_t ch=0; ch<438; ch++)
  {
    sprintf(Buff, "Veto_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    printf("%4.0f\n", Energy->GetMaximumBin()*Energy->GetBinWidth(Energy->GetMaximumBin()));
  }
}

void VetoEnergy()
{
  TH1F* Energy;
  Char_t Buff[256];
  Double_t MaxVal;
  Int_t MaxBin;

  for(Int_t ch=0; ch<438; ch++)
  {
    sprintf(Buff, "Veto_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    MaxVal = 0.0;
    MaxBin = 0;
    for(Int_t bn=60; bn<4096; bn++)
      if(Energy->GetBinContent(bn) > MaxVal)
      {
        MaxVal = Energy->GetBinContent(bn);
        MaxBin = bn;
      }
    printf("%4d\n", MaxBin);
  }
}

void VetoEnergy2(Int_t Lo, Int_t Hi)
{
  TCanvas* Window;
  TTimer* Refresh;
  TH1F* Energy;
  Char_t Buff[256];
  Char_t* Keyb;
  Double_t Pos[438];
  Double_t MaxVal;
  Int_t MaxBin;
  Int_t Start;

  Window = new TCanvas();
  Window->SetCrosshair();
  Window->ToggleEventStatus();
  Refresh = new TTimer("Flush()", 50, kFALSE);
  for(Int_t ch=Lo; ch<Hi+1; ch++)
  {
    sprintf(Buff, "Veto_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    Energy->GetXaxis()->SetRange(0, 500);
    Energy->Draw();
    Window->Update();
    sprintf(Buff, "Lower bound for %d: ", ch);
    Refresh->TurnOn();
    Refresh->Reset();
    Keyb = Getline(Buff);
    Refresh->TurnOff();
    Start = atoi(Keyb);
    MaxVal = 0.0;
    Pos[ch] = 0.0;
    for(Int_t bn=Start; bn<2048; bn++)
      if(Energy->GetBinContent(bn) > MaxVal)
      {
        MaxVal = Energy->GetBinContent(bn);
        Pos[ch] = 1.0*bn;
      }
  }
  for(Int_t ch=Lo; ch<Hi+1; ch++)
    printf("%4.0f\n", Pos[ch]);
}

void VetoEnergy3(Int_t Lo, Int_t Hi)
{
  TCanvas* Window;
  TTimer* Refresh;
  TH1F* Energy;
  TF1* Fit;
  Char_t Buff[256];
  Char_t* Keyb;
  Double_t Pos[438];
  Double_t Lower;
  Double_t Param[3];
  Int_t OK;

  Window = new TCanvas();
  Window->SetCrosshair();
  Window->ToggleEventStatus();
  Refresh = new TTimer("Flush()", 50, kFALSE);
  for(Int_t ch=Lo; ch<Hi+1; ch++)
  {
    sprintf(Buff, "Veto_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    Energy->Draw();
    Window->Update();

    Fit = new TF1("F1", "gaus", 0.7, 2.0);
    Energy->Fit("F1", "RQ");
    Fit->Draw("same");
    Window->Update();
    Fit->GetParameters(Param);

    Refresh->TurnOn();
    Refresh->Reset();
    sprintf(Buff, "OK for %d? ", ch);
    Keyb = Getline(Buff);
    Refresh->TurnOff();
    OK = atoi(Keyb);
    if(OK) goto Done;

Refit:
    Refresh->TurnOn();
    Refresh->Reset();
    sprintf(Buff, "Lower bound for %d: ", ch);
    Keyb = Getline(Buff);
    Refresh->TurnOff();
    Lower = atof(Keyb);

    Fit = new TF1("F1", "gaus", Lower, 2.0);
    Energy->Fit("F1", "RQ");
    Fit->Draw("same");
    Window->Update();
    Fit->GetParameters(Param);

    Refresh->TurnOn();
    Refresh->Reset();
    sprintf(Buff, "OK for %d? ", ch);
    Keyb = Getline(Buff);
    Refresh->TurnOff();
    OK = atoi(Keyb);
    if(!OK) goto Refit;

Done:
    Pos[ch] = Param[1];
  }
  for(Int_t ch=Lo; ch<Hi+1; ch++)
    printf("%8.6f\n", Pos[ch]);
}

void VetoThres()
{
  TH1F* Energy;
  Char_t Buff[256];
  Double_t MinVal;
  Int_t MinBin;

  for(Int_t ch=0; ch<438; ch++)
  {
    sprintf(Buff, "Veto_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    MinVal = 1e38;
    MinBin = 0;
    for(Int_t bn=1; bn<200; bn++)
      if(Energy->GetBinContent(bn) < MinVal)
      {
        MinVal = Energy->GetBinContent(bn);
        MinBin = bn;
      }
    printf("%4d\n", MinBin);
  }
}

void VetoThres2(Int_t Lo, Int_t Hi)
{
  TCanvas* Window;
  TTimer* Refresh;
  TH1F* Energy;
  Char_t Buff[256];
  Char_t* Keyb;
  Double_t Pos[438];

  Window = new TCanvas();
  Window->SetCrosshair();
  Window->ToggleEventStatus();
  Refresh = new TTimer("Flush()", 50, kFALSE);
  for(Int_t ch=Lo; ch<Hi+1; ch++)
  {
    sprintf(Buff, "Veto_Energy%d", ch);
    Energy = (TH1F*)gROOT->FindObject(Buff);
    Energy->GetXaxis()->SetRange(0, 250);
    Energy->Draw();
    Window->Update();
    sprintf(Buff, "Threshold for %d: ", ch);
    Refresh->TurnOn();
    Refresh->Reset();
    Keyb = Getline(Buff);
    Refresh->TurnOff();
    Pos[ch] = atof(Keyb);
  }
  for(Int_t ch=Lo; ch<Hi+1; ch++)
    printf("%5.1f\n", Pos[ch]);
}

//-----------------------------------------------------------------------------

void Flush()
{
  gSystem->ProcessEvents();
}
//-----------------------------------------------------------------------------

