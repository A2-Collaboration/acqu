void Pi0Minv(Char_t* Spect)
{
  TH1F* Hist = gROOT->FindObject(Spect);
  Double_t par[5];

  TF1* f2 = new TF1("f2", "expo", 70, 400);
  TF1* f3 = new TF1("f3", "gaus", 100, 150);

  f2->SetLineWidth(1);
  f2->SetLineColor(8);
  Hist->Fit(f2,"R0");
  f3->SetLineWidth(1);
  f3->SetLineColor(6);
  Hist->Fit(f3,"R0+");

  TF1 *to= new TF1("to", "expo(0)+gaus(2)", 70, 400);

  f2->GetParameters(&par[0]);
  f3->GetParameters(&par[2]);

  to->SetParameters(par);
  to->SetLineWidth(2);
  to->SetLineColor(2);
  Hist->Fit(to,"R+");
}

void EtaMinv(Char_t* Spect)
{
  TH1F* Hist = gROOT->FindObject(Spect);
  Double_t par[5];

  TF1* f2 = new TF1("f2", "expo", 200, 700);
  TF1* f3 = new TF1("f3", "gaus", 500, 600);

  f2->SetLineWidth(1);
  f2->SetLineColor(8);
  Hist->Fit(f2,"R0");
  f3->SetLineWidth(1);
  f3->SetLineColor(6);
  Hist->Fit(f3,"R0+");

  TF1 *to= new TF1("to", "expo(0)+gaus(2)", 200, 700);

  f2->GetParameters(&par[0]);
  f3->GetParameters(&par[2]);

  to->SetParameters(par);
  to->SetLineWidth(2);
  to->SetLineColor(2);
  Hist->Fit(to,"R+");
}
