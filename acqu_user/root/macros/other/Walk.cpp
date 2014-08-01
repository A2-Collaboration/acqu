Double_t fitt1(Double_t *x, Double_t *par)
{
  return (par[0] + par[1] * TMath::Power(x[0], par[2]));
}

//-----------------------------------------------------------------------------

void WalkFit1(int Crys)
{
  char Name[255];
  TH1F* Spect1;
  TH2F* Spect2;

  Double_t x;
  Double_t y1;
  Double_t y2;

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = +1.8;

  sprintf(&Name[0], "NaI_Time%d_v_Energy%d\0", Crys, Crys);
  Spect2 = (TH2F)gROOT->FindObject(Name);
  Spect1 = (TH1F)Spect2->ProfileX();

  TF1* f1 = new TF1("f1", fitt1, 0, 800, 3);
  par1[0] = 35.2398;
  par1[1] = 75.3237;
  par1[2] = -0.571325;

  f1->SetParameters(&par1[0]);
  f1->SetLineWidth(1);
  f1->SetLineColor(2);
  Spect1->Fit(f1,"RQ");

  TF1* f2 = new TF1("f2", "pol0", 150, 800);
  par2 = 39.1234;

  f2->SetParameters(&par2);
  f2->SetLineWidth(1);
  f2->SetLineColor(8);
  Spect1->Fit(f2,"RQ+");

  x = 0.0;
  y1 = f1->Eval(x);
  y2 = f2->Eval(x);;
  while(y1>y2)
  {
    x+=0.1;
    y1 = f1->Eval(x);
    y2 = f2->Eval(x);
    if(x>400.0) break;
  }

  cutoff = x;
  cutoff = 400.0;
  f1->GetParameters(&par1[0]);
  f2->GetParameters(&par2);

  printf("TimeWalk: %3d %12.6f %12.6f %12.6f %7.2f %12.6f %7.2f\n", Crys, par1[0], par1[1], par1[2], cutoff, par2, shift);
}

//-----------------------------------------------------------------------------

void WalkFit2()
{
  char Name[256];
  TH1F* Spect1;
  TH2F* Spect2;

  Double_t x;
  Double_t y1;
  Double_t y2;

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = 75.0;

  FILE* File1;
  FILE* File2;
  char Line[1024];
  File1 = fopen("/home/schumans/3v16/acqu/root/macros/NaITemplate.txt","r");
  File2 = fopen("/home/schumans/3v16/acqu/NaI.txt","w");
  fclose(File2);

  for(Int_t Crys=0; Crys<720; Crys++)
  {
    sprintf(&Name[0], "NaI_Time%d_v_Energy%d\0", Crys, Crys);
    Spect2 = (TH2F)gROOT->FindObject(Name);
    Spect1 = (TH1F)Spect2->ProfileX();

    TF1* f1 = new TF1("f1", fitt1, 0, 400, 3);
    par1[0] = 35.2398;
    par1[1] = 75.3237;
    par1[2] = -0.571325;

    f1->SetParameters(&par1[0]);
    f1->SetLineWidth(1);
    f1->SetLineColor(2);
    Spect1->Fit(f1,"RQ");

    TF1* f2 = new TF1("f2", "pol0", 150, 400);
    par2 = 39.1234;

    f2->SetParameters(&par2);
    f2->SetLineWidth(1);
    f2->SetLineColor(8);
    Spect1->Fit(f2,"RQ+");

    x = 0.0;
    y1 = f1->Eval(x);
    y2 = f2->Eval(x);;
    while(y1>y2)
    {
      x+=0.1;
      y1 = f1->Eval(x);
      y2 = f2->Eval(x);
      if(x>400.0) break;
    }

    cutoff = x;//188.6;
    f1->GetParameters(&par1[0]);
    f2->GetParameters(&par2);

    fgets(&Line[0], 1023, File1);
    File2 = fopen("/home/schumans/3v16/acqu/NaI.txt","a");
    fprintf(File2, "%sTimeWalk:  %f %f %f %f %f %f\n", Line, par1[0], par1[1], par1[2], cutoff, par2, shift);
    fclose(File2);
  }

  fclose(File1);
}

//-----------------------------------------------------------------------------

void Schieb()
{
  char Name[256];
  TH1F* Spect1;
  Double_t max;
  Double_t max2;
  Double_t par[4];

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = 75.0;

  FILE* File1;
  FILE* File2;
  char Line0[1024];
  char Line1[1024];
  char Out[1024];
  File1 = fopen("/home/schumans/3v16a/acqu/times.txt","r");
  File2 = fopen("/home/schumans/3v16a/acqu/times2.txt","w");
  fclose(File2);


  for(Int_t t=0; t<720; t++)
  {
    sprintf(&Name[0], "Sven_CBTagg%d\0", t);
    Spect1 = (TH1F)gROOT->FindObject(Name);
    max2 = Spect1->GetMaximumBin()*Spect1->GetBinWidth(10) - 100.0;

    max = Spect1->GetMaximumBin()*Spect1->GetBinWidth(10)- 100.0;
    TF1* f2 = new TF1("f2", "pol0", max-35, max+70);
    Spect1->Fit(f2,"RQ");
    TF1* f3 = new TF1("f3", "gaus", max-6, max+6);
    Spect1->Fit(f3,"RQ+");
    f2->GetParameters(&par[0]);
    f3->GetParameters(&par[1]);

    TF1 *to= new TF1("to", "pol0(0)+gaus(1)", max-40, max+90);
    to->SetParameters(par);
    Spect1->Fit(to,"RQ");
    to->GetParameters(&par[0]);

    max = par[2];

    fgets(&Line0[0], 1023, File1);
    fgets(&Line1[0], 1023, File1);

    File2 = fopen("/home/schumans/3v16a/acqu/times2.txt","a");
    sscanf(Line1, "%s %lf %lf %lf %lf %lf %lf", &Out, &par1[0], &par1[1], &par1[2], &cutoff, &par2, &shift);
    shift = shift - max2;
    fprintf(File2,"%s", Line0);
    fprintf(File2,"TimeWalk:  %f %f %f %f %f %f\n", par1[0], par1[1], par1[2], cutoff, par2, shift);
    fclose(File2);

    cout << t << ": " << -max << " <-> " << -max2  << endl;
  }

}

//-----------------------------------------------------------------------------

void WalkFit3()
{
  char Name[256];
  TH1F* Spect1;
  TH2F* Spect2;

  Double_t x;
  Double_t y1;
  Double_t y2;

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = -70.6;

  FILE* File1;
  FILE* File2;
  char Line[1024];
  File1 = fopen("/home/schumans/3v17/data/jul2walk/nowalk","r");
  File2 = fopen("/home/schumans/3v17/data/jul2walk/walk","w");
printf("%d\n",File2);
  fclose(File2);

  for(Int_t Crys=0; Crys<720; Crys++)
  {
    sprintf(&Name[0], "Sven_CBTagg%d_v_CB_NaI_Energy%d\0", Crys, Crys);
    Spect2 = (TH2F)gROOT->FindObject(Name);
    Spect1 = (TH1F)Spect2->ProfileX();

    TF1* f1 = new TF1("f1", fitt1, 0.0, 400.0, 3);
    par1[0] = 35.2398;
    par1[1] = 75.3237;
    par1[2] = -0.571325;

    f1->SetParameters(&par1[0]);
    f1->SetLineWidth(1);
    f1->SetLineColor(2);
    Spect1->Fit(f1,"RQ");

    /*TF1* f2 = new TF1("f2", "pol0", 150.0, 400.0);
    par2 = 39.1234;

    f2->SetParameters(&par2);
    f2->SetLineWidth(1);
    f2->SetLineColor(8);
    Spect1->Fit(f2,"RQ+");

    x = 0.0;
    y1 = f1->Eval(x);
    y2 = f2->Eval(x);;
    while(y1>y2)
    {
      x+=0.1;
      y1 = f1->Eval(x);
      y2 = f2->Eval(x);
      if(x>400.0) break;
    }*/

    //cutoff = x;
    cutoff = 800.0;
    f1->GetParameters(&par1[0]);
    //f2->GetParameters(&par2);
    par2 = 0.0;

    fgets(&Line[0], 1023, File1);
    File2 = fopen("/home/schumans/3v17/data/jul2walk/walk","a");
    fprintf(File2, "%sTimeWalk:  %f %f %f %f %f %f\n", Line, par1[0], par1[1], par1[2], cutoff, par2, shift);
    fclose(File2);
    printf("TimeWalk:  %f %f %f %f %f %f\n", par1[0], par1[1], par1[2], cutoff, par2, shift);
  }

  fclose(File1);
}

//-----------------------------------------------------------------------------

void WalkFit4(int Crys, double a, double b)
{
  char Name[255];
  TH1F* Spect1;
  TH2F* Spect2;

  Double_t x;
  Double_t y1;
  Double_t y2;

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = -72.2;

  sprintf(&Name[0], "Sven_CBTagg%d_v_CB_NaI_Energy%d\0", Crys, Crys);
  Spect2 = (TH2F)gROOT->FindObject(Name);
  Spect1 = (TH1F)Spect2->ProfileX();

  TF1* f1 = new TF1("f1", fitt1, a, b, 3);
  par1[0] = 35.2398;
  par1[1] = 75.3237;
  par1[2] = -0.571325;

  f1->SetParameters(&par1[0]);
  f1->SetLineWidth(1);
  f1->SetLineColor(2);
  Spect1->Fit(f1,"RQ");

  /*TF1* f2 = new TF1("f2", "pol0", 150.0, 400.0);
  par2 = 39.1234;

  f2->SetParameters(&par2);
  f2->SetLineWidth(1);
  f2->SetLineColor(8);
  Spect1->Fit(f2,"RQ+");

  x = 0.0;
  y1 = f1->Eval(x);
  y2 = f2->Eval(x);;
  while(y1>y2)
  {
    x+=0.1;
    y1 = f1->Eval(x);
    y2 = f2->Eval(x);
    if(x>400.0) break;
  }*/

  //cutoff = x;
  cutoff = 800.0;
  f1->GetParameters(&par1[0]);
  //f2->GetParameters(&par2);
  par2 = 0.0;

  printf("TimeWalk:  %f %f %f %f %f %f\n", par1[0], par1[1], par1[2], cutoff, par2, shift);
}

//-----------------------------------------------------------------------------

void machhin()
{
  char Name[256];
  TH1F* Spect1;
  Double_t max;
  Double_t max2;
  Double_t par[4];

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = 75.0;

  FILE* File1;
  FILE* File2;
  char Line0[1024];
  char Line1[1024];
  char Out[1024];
  File1 = fopen("/home/schumans/3v17/acqu/cal.txt","r");
  File2 = fopen("/home/schumans/3v17/acqu/cal2.txt","w");
  fclose(File2);

  for(Int_t t=0; t<720; t++)
  {
    fgets(&Line0[0], 1023, File1);
    fgets(&Line1[0], 1023, File1);

    File2 = fopen("/home/schumans/3v17/acqu/cal2.txt","a");
    fprintf(File2,"%s", Line0);
    fclose(File2);
  }

}


//-----------------------------------------------------------------------------

void Hui()
{
  Double_t par1[3];
  Double_t shift = -71.9;

  TF1* f1 = new TF1("f1", fitt1, 0.0, 400.0, 3);
  par1[0] = 24.034073;
  par1[1] = 54.729258;
  par1[2] = -0.306611;

  f1->SetParameters(&par1[0]);
  f1->SetLineWidth(4);
  f1->SetLineColor(1);
  f1->Draw("same");

  TF1* f2 = new TF1("f2", fitt1, 0.0, 400.0, 3);

  f2->SetParameters(&par1[0]);
  f2->SetLineWidth(2);
  f2->SetLineColor(10);
  f2->Draw("same");
}

//-----------------------------------------------------------------------------

void Hui2(Int_t Bin, Double_t a, Double_t b)
{
  Double_t par2[4];
  Char_t Filename[255];
  sprintf(Filename, "/home/schumann/3v17/acqu/cuts/PCorrCB_04Jul_LP883_2.txt");
  FILE* Params = fopen(Filename, "r");
  Int_t n;
  while(!feof(Params))
  {
    fscanf(Params, "%d %lf %lf %lf %lf",&n, &par2[0], &par2[1], &par2[2], &par2[3]);
    if(n==Bin) break;
  }

  TF1* f1 = new TF1("f1", "pol3", a, b);

  f1->SetParameters(&par2[0]);
  f1->SetLineWidth(4);
  f1->SetLineColor(1);
  f1->Draw("same");

  TF1* f2 = new TF1("f2", "pol3", a, b);

  f2->SetParameters(&par2[0]);
  f2->SetLineWidth(2);
  f2->SetLineColor(10);
  f2->Draw("same");
}

//-----------------------------------------------------------------------------

void Hui3(Int_t Bin, Double_t a, Double_t b)
{
  Double_t par2[4];
  Char_t Filename[255];
  sprintf(Filename, "/home/schumann/3v17/acqu/cuts/PCorrTAPS_All.txt");
  FILE* Params = fopen(Filename, "r");
  Int_t n;
  while(!feof(Params))
  {
    fscanf(Params, "%d %lf %lf %lf %lf",&n, &par2[0], &par2[1], &par2[2], &par2[3]);
    if(n==Bin) break;
  }

  TF1* f1 = new TF1("f1", "pol3", a, b);

  f1->SetParameters(&par2[0]);
  f1->SetLineWidth(4);
  f1->SetLineColor(1);
  f1->Draw("same");

  TF1* f2 = new TF1("f2", "pol3", a, b);

  f2->SetParameters(&par2[0]);
  f2->SetLineWidth(2);
  f2->SetLineColor(10);
  f2->Draw("same");
}

//-----------------------------------------------------------------------------

void Ch17()
{
  char Name[255];
  int Crys = 17;
  TH1F* Spect1;
  TH2F* Spect2;

  Double_t x;
  Double_t y1;
  Double_t y2;

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = +1.8;

  sprintf(&Name[0], "NaI_Time%d_v_Energy%d\0", Crys, Crys);
  Spect2 = (TH2F)gROOT->FindObject(Name);
  Spect1 = (TH1F)Spect2->ProfileX();

  //TF1* f1 = new TF1("f1", fitt1, 0, 800, 3);
  TF1* f1 = new TF1("f1", fitt1, 0, 45, 3);
  par1[0] = 35.2398;
  par1[1] = 75.3237;
  par1[2] = -0.571325;

  f1->SetParameters(&par1[0]);
  f1->SetLineWidth(1);
  f1->SetLineColor(2);
  Spect1->Fit(f1,"RQ");

  //TF1* f2 = new TF1("f2", "pol0", 150, 800);
  TF1* f2 = new TF1("f2", "pol0", 60, 800);
  par2 = 39.1234;

  f2->SetParameters(&par2);
  f2->SetLineWidth(1);
  f2->SetLineColor(8);
  Spect1->Fit(f2,"RQ+");

  x = 0.0;
  y1 = f1->Eval(x);
  y2 = f2->Eval(x);;
  while(y1>y2)
  {
    x+=0.1;
    y1 = f1->Eval(x);
    y2 = f2->Eval(x);
    if(x>400.0) break;
  }

  cutoff = x;
  cutoff = 45.0;
  f1->GetParameters(&par1[0]);
  f2->GetParameters(&par2);

  printf("TimeWalk: %3d %12.6f %12.6f %12.6f %7.2f %12.6f %7.2f\n", Crys, par1[0], par1[1], par1[2], cutoff, par2, shift);
}

//-----------------------------------------------------------------------------

void Ch672()
{
  char Name[255];
  int Crys = 672;
  TH1F* Spect1;
  TH2F* Spect2;

  Double_t x;
  Double_t y1;
  Double_t y2;

  Double_t par1[3];
  Double_t par2;
  Double_t cutoff;
  Double_t shift = +1.8;

  sprintf(&Name[0], "NaI_Time%d_v_Energy%d\0", Crys, Crys);
  Spect2 = (TH2F)gROOT->FindObject(Name);
  Spect1 = (TH1F)Spect2->ProfileX();

  //TF1* f1 = new TF1("f1", fitt1, 0, 800, 3);
  TF1* f1 = new TF1("f1", fitt1, 6, 45, 3);
  par1[0] = 35.2398;
  par1[1] = 75.3237;
  par1[2] = 0.571325;

  f1->SetParameters(&par1[0]);
  f1->SetLineWidth(1);
  f1->SetLineColor(2);
  Spect1->Fit(f1,"RQ");

  //TF1* f2 = new TF1("f2", "pol0", 150, 800);
  TF1* f2 = new TF1("f2", "pol0", 60, 800);
  par2 = 39.1234;

  f2->SetParameters(&par2);
  f2->SetLineWidth(1);
  f2->SetLineColor(8);
  Spect1->Fit(f2,"RQ+");

  x = 0.0;
  y1 = f1->Eval(x);
  y2 = f2->Eval(x);;
  while(y1>y2)
  {
    x+=0.1;
    y1 = f1->Eval(x);
    y2 = f2->Eval(x);
    if(x>400.0) break;
  }

  cutoff = x;
  cutoff = 45.0;
  f1->GetParameters(&par1[0]);
  f2->GetParameters(&par2);

  printf("TimeWalk: %3d %12.6f %12.6f %12.6f %7.2f %12.6f %7.2f\n", Crys, par1[0], par1[1], par1[2], cutoff, par2, shift);
}

//-----------------------------------------------------------------------------

void CheckCoinc(Char_t* Meson)
{
  TH1F* Spect;
  Char_t Histo[256];
  Double_t par[4];
  Double_t Max;

  sprintf(Histo, "%s_EventTime", Meson);
  Spect = (TH1F)gROOT->FindObject(Histo);
  Max = Spect->GetMaximumBin() * Spect->GetBinWidth(Spect->GetMaximumBin()) - 500.0;
  TF1* f2 = new TF1("f2", "pol0", Max-50.0, Max+50.0);
  Spect->Fit(f2,"RQ");
  TF1* f3 = new TF1("f3", "gaus", Max-5.0, Max+5.0);
  Spect->Fit(f3,"RQ+");
  f2->GetParameters(&par[0]);
  f3->GetParameters(&par[1]);

  TF1 *to= new TF1("to", "pol0(0)+gaus(1)", Max-50.0, Max+50.0);
  to->SetParameters(par);
  Spect->Fit(to,"RQ");
  to->GetParameters(&par[0]);
  printf("%5.3f =- %5.3f\n", par[2], par[3]);
}

//-----------------------------------------------------------------------------
