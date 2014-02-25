void CheckLinPol(){
TString name("LinPolSpectra");
TCanvas *canvLinPol = new TCanvas(name, name, 200, 60, 700, 1200);
//gStyle->SetOptFit(1);
//gStyle->SetOptStat(0);
   canvLinPol->Range(-79.0074,-0.373242,79.3233,0.359577);
   canvLinPol->SetFillColor(10);
   canvLinPol->SetFillStyle(4000);
   canvLinPol->SetBorderSize(0);
   canvLinPol->SetFrameFillColor(10);
   canvLinPol->SetFrameFillStyle(0);
   canvLinPol->SetBottomMargin(0.12);
   canvLinPol->SetTopMargin(0.03);
//canvLinPol->SetStatFormat("5.3g");

TA2LinearPol* pLP = (TA2LinearPol*)(gAN->GetChild("LinPol"));

char* hnames[] = {"LinPol_Incoherent","LinPol_Coherent","LinPol_CoherentPara","LinPol_CoherentPerp","LinPol_EnhancementPara","LinPol_EnhancementPerp","LinPol_CohEdgeEnergy"};

TLatex *texPARA= new TLatex(0.63,0.85,"PARA");;
texPARA->SetNDC();
texPARA->SetTextColor(2);
TLatex *texPERP= new TLatex(0.63,0.8,"PERP");;
texPERP->SetNDC();
texPERP->SetTextColor(4);

//canvLinPol->SetFillStyle(4000);
canvLinPol->Divide(2,3,0.01,0.01);
TH1F* his;

/////////LinPol_Incoherent
canvLinPol->cd(1);
his=(TH1F*)gROOT->Get(hnames[0]);
his->Draw();

/////////LinPol_Coherent
canvLinPol->cd(2);
his=(TH1F*)gROOT->Get(hnames[1]);
his->Draw();
Int_t mode=pLP->GetGoniMode();//Para or perp?
//cout<<"Current mode "<<mode<<endl;
if(mode==1)
{
texPARA->Draw();
his->SetLineColor(2);
}
else if (mode==2)
{
texPERP->Draw();
his->SetLineColor(4);
}

/////////LinPol_CoherentPara
canvLinPol->cd(3);
his=(TH1F*)gROOT->Get(hnames[2]);
his->SetLineColor(2);
his->Draw();

/////////LinPol_CoherentPerp
canvLinPol->cd(4);
his=(TH1F*)gROOT->Get(hnames[3]);
his->SetLineColor(4);
his->Draw("");

/////////LinPol_EnhancementPara
canvLinPol->cd(5);
his=(TH1F*)gROOT->Get(hnames[4]);
 his->SetLineColor(2);

if(his->Integral()){
  his->Draw();
  texPARA->Draw();
}
//////////LinPol_CoherentPerp
//canvLinPol->cd(6);
if(his->Integral()){
  his=(TH1F*)gROOT->Get(hnames[5]);
  his->SetLineColor(4);
  his->Draw("same");
  texPERP->Draw();
}
else{
  his=(TH1F*)gROOT->Get(hnames[5]);
  if(his->Integral()){
  his->SetLineColor(4);
  his->Draw("");
  texPERP->Draw();
 }
}

//////////LinPol_CoherentPerp
canvLinPol->cd(6);
his=(TH1F*)gROOT->Get(hnames[6]);
his->Draw("");

}
