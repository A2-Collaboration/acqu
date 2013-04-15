void DblGauss(Char_t* Meson)
{
  Char_t Name[256];
  TF1* G = new TF1("Gauss","gaus", 497.5, 597.5);
  TF1* D = new TF1("DblGauss","gaus(0)+gaus(3)", 497.5, 597.5);
  sprintf(Name,"%s_MesonMinvW1", Meson);
  TH1F* M = (TH1F*)gROOT->FindObject(Name);
  M->Fit(G, "R");
  D->SetParameter(0, G->GetParameter(0));
  D->SetParameter(1, G->GetParameter(1));
  D->SetParameter(2, G->GetParameter(2));
  D->SetParameter(3, G->GetParameter(0) * 0.5);
  D->SetParameter(4, G->GetParameter(1) * 0.9);
  D->SetParameter(5, G->GetParameter(2) * 1.5);
  M->Fit(D, "R");
  printf("%f\n", D->GetMaximumX(497.5, 597.5));
}

void DblGaussExp(Char_t* Meson)
{
  Char_t Name[256];
  TF1* G = new TF1("Gauss","gaus", 497.5, 597.5);
  TF1* E = new TF1("Exp","expo", 200.0, 400.0);
  TF1* D = new TF1("DblGaussExp","expo(0)+gaus(2)+gaus(5)", 497.5, 597.5);
  sprintf(Name,"%s_MesonMinvW1", Meson);
  TH1F* M = (TH1F*)gROOT->FindObject(Name);
  M->Fit(G, "R");
  M->Fit(E, "R");
  D->SetParameter(0, E->GetParameter(0));
  D->SetParameter(1, E->GetParameter(1));
  D->SetParameter(2, G->GetParameter(0));
  D->SetParameter(3, G->GetParameter(1));
  D->SetParameter(4, G->GetParameter(2));
  D->SetParameter(5, G->GetParameter(0) * 0.5);
  D->SetParameter(6, G->GetParameter(1) * 0.9);
  D->SetParameter(7, G->GetParameter(2) * 1.5);
  M->Fit(D, "R");
  //printf("%f\n", D->GetMaximumX(497.5, 597.5));

  TF1* N = new TF1("DblGauss","gaus(0)+gaus(3)", 497.5, 597.5);
  N->SetParameter(0, D->GetParameter(2));
  N->SetParameter(1, D->GetParameter(3));
  N->SetParameter(2, D->GetParameter(4));
  N->SetParameter(3, D->GetParameter(5));
  N->SetParameter(4, D->GetParameter(6));
  N->SetParameter(5, D->GetParameter(7));
  printf("%f\n", N->GetMaximumX(497.5, 597.5));
}

void GaussExp(Char_t* Meson)
{
  Char_t Name[256];
  TF1* G = new TF1("Gauss","gaus", 110.0, 160.0);
  TF1* E = new TF1("Exp","expo", 70.0, 100.0);
  TF1* D = new TF1("GaussExp","expo(0)+gaus(2)", 110.0, 160.0);
  sprintf(Name,"%s_MesonMinvW1", Meson);
  TH1F* M = (TH1F*)gROOT->FindObject(Name);
  M->Fit(G, "R");
  M->Fit(E, "R");
  D->SetParameter(0, E->GetParameter(0));
  D->SetParameter(1, E->GetParameter(1));
  D->SetParameter(2, G->GetParameter(0));
  D->SetParameter(3, G->GetParameter(1));
  D->SetParameter(4, G->GetParameter(2));
  M->Fit(D, "R");
  //printf("%f\n", D->GetMaximumX(110.0, 160.0));

  TF1* N = new TF1("SglGauss","gaus(0)", 110.0, 160.0);
  N->SetParameter(0, D->GetParameter(2));
  N->SetParameter(1, D->GetParameter(3));
  N->SetParameter(2, D->GetParameter(4));
  printf("%f\n", N->GetMaximumX(110.0, 160.0));
}

void Gauss(Char_t* Meson)
{
  Char_t Name[256];
  TF1* G = new TF1("Gauss","gaus", 110.0, 160.0);
  sprintf(Name,"%s_MesonMinvW1", Meson);
  TH1F* M = (TH1F*)gROOT->FindObject(Name);
  M->Fit(G, "R");
  printf("%f\n", G->GetMaximumX(110.0, 160.0));
}
