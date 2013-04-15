TGraphErrors* TaggEffPlot;
TH1F* histEelec;             // electron energy
TH1F* histHits;              // hits histogram
TH1D* histSc;                // scaler histogram
TH1D* histScB;               // background scaler histogram
TH1D* histS;                 // trigger scalers histogram
TH1D* histSB;                // background trigger scalers histogram

TH1D* histScB2;               // background scaler histogram
TH1D* histSB2;                // background trigger scalers histogram

Double_t* Eelec;             // electron energy
Double_t* Hits;              // hits array
Double_t* Scalers;           // accumulated scalers
Double_t* ScalersB;          // accumulated scalers background run
Double_t* TaggEff;           // tagging efficiency
Double_t* dTaggEff;          // uncertainty
Double_t* TaggEffB;          // background corrected tagging efficiency
Double_t* dTaggEffB;         // uncertainty

Double_t meanEff, meanEffB, dmeanEff, dmeanEffB;

FILE* TaggEffLog;
FILE* SvenLog;
Char_t Buffer[256];


//void TaggEff(Char_t* ffore, Char_t* fback,
//	     Int_t nchan=352, Int_t jclk=16, Int_t jinh=17)
void TaggEff(Char_t* fback1, Char_t* ffore, Char_t* fback2,
         Int_t nchan=352, Int_t jclk=1, Int_t jinh=1)
{
  // calculate tagging efficiency
  Eelec = new Double_t[nchan];
  Hits = new Double_t[nchan];
  Scalers = new Double_t[nchan];
  ScalersB = new Double_t[nchan];
  TaggEff = new Double_t[nchan];
  dTaggEff = new Double_t[nchan];
  TaggEffB = new Double_t[nchan];
  dTaggEffB = new Double_t[nchan];

  TFile* fore = new TFile(ffore);
  if( !fore ){ printf("File %s not found\n", ffore) return; }
  TFile* back1 = new TFile(fback1);
  if( !back1 ){ printf("File %s not found\n", fback1) return; }
  TFile* back2 = new TFile(fback2);
  if( !back2 ){ printf("File %s not found\n", fback2) return; }
  if( !(TaggEffLog = fopen("./TaggEff.log","w+")) ){
    printf("Unable to open TaggEff.log\n");
    return;
  }
  sprintf(Buffer, "%s.txt", ffore);
  SvenLog = fopen(Buffer, "w");

  histEelec =  (TH1F*)fore->Get("Ladder_Eelect");
  histHits = (TH1F*)fore->Get("Ladder_Hits");
  histSc = (TH1D*)fore->Get("Ladder_ScalerAcc");
  histScB = (TH1D*)back1->Get("Ladder_ScalerAcc");
  //histS = (TH1D*)fore->Get("SumScalers0to43");
  //histSB = (TH1D*)back1->Get("SumScalers0to43");
  histS = (TH1D*)fore->Get("SumScalers0to11");
  histSB = (TH1D*)back1->Get("SumScalers0to11");

  histScB2 = (TH1D*)back2->Get("Ladder_ScalerAcc");
  histSB2 = (TH1D*)back2->Get("SumScalers0to11");

  if( !histEelec ) printf("Electron energy hist not found\n");
  if( !histHits ) printf("Hits hist not found\n");
  if( !histSc ) printf("Scalers hist not found\n");
  if( !histS ) printf("Trig-Scalers hist not found\n");
  if( !histScB ) printf("Background Scalers hist not found\n");
  if( !histSB ) printf("Background Trig-Scalers hist not found\n");

  Double_t ratio12 = histSB->GetBinContent(jclk)/histSB2->GetBinContent(jclk);
  histScB2->Scale(ratio12);
  histSB2->Scale(ratio12);
  histScB->Add(histScB2);
  histSB->Add(histSB2);
  
  Double_t clk = (Double_t)histS->GetBinContent(jclk);
  Double_t inh = (Double_t)histS->GetBinContent(jinh);
  Double_t live = clk/inh;
  live = 1.0;
  //TH1F* HistLive = (TH1F*)fore->Get("Ladder_Livetime");
  //live = HistLive->GetMean();
  Double_t clkB = (Double_t)histSB->GetBinContent(jclk);
  Double_t bscale = clk/clkB;
  Double_t ScB, dScB;

  fprintf(TaggEffLog,
	 " Tagging efficiency calculation\n J.R.M.Annand 19/11/06\n\n");

  fprintf(TaggEffLog,"Clock = %10g   Inhibited Clock = %10g   Live = %8g\n\
Background clock = %10g    background scale factor = %8g\n\n",
	  clk,inh,1/live,clkB,bscale);

  fprintf(TaggEffLog,
"   Gi  Hits   Scaler  Background    TaggEff           TaggEff(back-sub)\n");

  for(UInt_t i = 0; i< nchan; i++){
    //    Eelec[i] = (Double_t)histEelec->GetBinContent(i);
    Eelec[i] = (Double_t)i;
    Hits[i]  = (Double_t)histHits->GetBinContent(i);
    Scalers[i] = (Double_t)histSc->GetBinContent(i);
    ScalersB[i] = (Double_t)histScB->GetBinContent(i);

    if( (Scalers[i] == 0) || (Hits[i] == 0) ){
      TaggEff[i] = TaggEffB[i] = dTaggEff[i] = dTaggEffB[i] = 0.0;
    }
    else{
      ScB = Scalers[i] - bscale*ScalersB[i];
      dScB = TMath::Sqrt(ScalersB[i]) * bscale;
      dScB = TMath::Sqrt( Scalers[i] + dScB*dScB )/ScB;
      
      TaggEff[i] = live * Hits[i]/Scalers[i];
      TaggEffB[i] = live * Hits[i]/( Scalers[i] - bscale*ScalersB[i] );
      dTaggEff[i] = TaggEff[i] * TMath::Sqrt( 1./Hits[i] + 1./Scalers[i] );
      dTaggEffB[i] = TaggEffB[i] * TMath::Sqrt( 1./Hits[i] + dScB*dScB );
    }
    fprintf(TaggEffLog,
	    "%5g %5g %8g %8g   %7.3g +-%6.2g  %7.3g +-%6.2g\n",
	   Eelec[i], Hits[i], Scalers[i], ScalersB[i],
	   TaggEff[i], dTaggEff[i], TaggEffB[i], dTaggEffB[i]);
    fprintf(SvenLog, "%3d    %10.1f  %10.4f    %10.2f  %10.3f    %9.6f  %9.6f\n",
            i,
            live*Hits[i], live*TMath::Sqrt(Hits[i]),
            Scalers[i]-bscale*ScalersB[i], TMath::Sqrt(Scalers[i]+bscale*ScalersB[i]),
            TaggEffB[i], dTaggEffB[i]);
  }
  TaggEffDraw(nchan,ffore);
  fprintf(TaggEffLog,
   "\nMean tagging efficiency %7.3g +- %6.2g\n\
                        %7.3g +- %6.2g without background subtraction\n",
	  meanEffB, dmeanEffB, meanEff, dmeanEff );
  //fprintf(SvenLog, "999   %8.6f   %8.6f\n", meanEffB, dmeanEffB);

  fclose(TaggEffLog);
  fclose(SvenLog);
  
  fore->Close();
  back1->Close();
  back2->Close();
}

void TaggEffDraw(Int_t nchan, Char_t* run = ""){
  // refills the TeffData into gr_data avoiding "bad" channels and
  // plots the graph.
  TCanvas* canv;
  if( !( canv = (TCanvas*)(gROOT->FindObject("Teff")) ) ){
    TCanvas* canv = new TCanvas("Teff","Tagging Efficiency",240,180,700,900);
    canv->SetFillStyle(4000);
    canv->Divide(1,2,0.01,0.01);
    canv->cd(1);
  }
  else canv->cd(1);
  char title[255];
  UInt_t j = 0;
  for(UInt_t i=0; i<nchan; i++){ 
    if( TaggEff[i] > 0 ) j++;
  }

  Double_t* ee = new Double_t[j];
  Double_t* tf = new Double_t[j];
  Double_t* dtf = new Double_t[j];
  j = 0;
  meanEff = dmeanEff = 0.0;
  for(UInt_t i=0; i<nchan; i++){
    if( TaggEff[i] > 0 ){
      ee[j] = Eelec[i];
      tf[j] = TaggEff[i];
      meanEff = meanEff + tf[j];
      dtf[j] = dTaggEff[i];
      dmeanEff = dmeanEff + dtf[j]*dtf[j];
      j++;
    }
  }
  meanEff = meanEff/j;
  dmeanEff = TMath::Sqrt(dmeanEff)/j;

  TaggEffPlot = new TGraphErrors(j, ee, tf, NULL, dtf);
  strcpy(title,"Tagging-Efficiency");
  if (strlen(run)){ 
    strcat(title, " of ");
    strcat(title, run);
  }
  strcat(title, " (Preliminary)");
  TaggEffPlot->SetTitle(title);
  TaggEffPlot->GetXaxis()->SetTitle("Tagger Channel (Goosey)");
  TaggEffPlot->GetYaxis()->SetTitle("#epsilon_{tagg}");
  TaggEffPlot->GetXaxis()->CenterTitle();
  TaggEffPlot->GetYaxis()->CenterTitle();
  TaggEffPlot->SetMarkerColor(kBlue);
  TaggEffPlot->SetMarkerStyle(21);
  TaggEffPlot->SetMarkerSize(0.7);
  TaggEffPlot->Draw("AP");

  //
  canv->cd(2);
  j = 0;
  meanEffB = 0.0;
  for(UInt_t i=0; i<nchan; i++){
    if( TaggEff[i] > 0 ){
      ee[j] = Eelec[i];
      tf[j] = TaggEffB[i];
      meanEffB = meanEffB + tf[j];
      dtf[j] = dTaggEffB[i];
      dmeanEffB = dmeanEffB + dtf[j]*dtf[j];
      j++;
    }
  }
  meanEffB = meanEffB/j;
  dmeanEffB = TMath::Sqrt(dmeanEffB)/j;

  TaggEffPlotB = new TGraphErrors(j, ee, tf, NULL, dtf);
  strcpy(title,"Background-Corrected Tagging-Efficiency");
  if (strlen(run)){ 
    strcat(title, " of ");
    strcat(title, run);
  }
  strcat(title, " (Preliminary)");
  TaggEffPlotB->SetTitle(title);
  TaggEffPlotB->GetXaxis()->SetTitle("Tagger Channel (Goosey)");
  TaggEffPlotB->GetYaxis()->SetTitle("#varepsilon_{tagg}");
  TaggEffPlotB->GetXaxis()->CenterTitle();
  TaggEffPlotB->GetYaxis()->CenterTitle();
  TaggEffPlotB->SetMarkerColor(kRed);
  TaggEffPlotB->SetMarkerStyle(20);
  TaggEffPlotB->SetMarkerSize(0.7);
  TaggEffPlotB->Draw("AP");
}
