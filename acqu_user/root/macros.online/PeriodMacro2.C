// This code gets executed by the online analysis
// every Nth events (see online config Period: statement)
int x=5;
int  bufferCounter = 0;

void PeriodMacro() {
  // Dump the current FP scalers to a nfs exported file system
  // if(gROOT->FindObject("FPD_ScalerCurr")){
  //   if((FPD_ScalerCurr->Integral()) > 0){
  //     FILE  *fp=fopen("scratch/tagscaler.txt","w");
  //     //FILE  *fp=fopen("/exports/a2raid5/tagscaler.txt","w");
  //     for(int n=1;n<=352;n++){
  // 	fprintf(fp,"%d\n", FPD_ScalerCurr->GetBinContent(n));
  //     }
  //     fclose(fp);
  //   }
  // }

  //replaced writing out the tagger scalers to file, with writing to EPICS.
   // write the tagger scalers to data

  cout << "Period Mac 1 " << endl;


  if(gROOT->FindObject("FPD_ScalerCurr")){
    if((FPD_ScalerCurr->Integral())>0) {
      stringstream cmd;
      cmd << "caput -a TAGGER:RAW_SCALER 352";
      for(int n=352; n>=1; n--) {
	cmd << FPD_ScalerCurr->GetBinContent(n) << " ";
      } 
      cmd << " > /dev/null";
      system(cmd.str().c_str());
    }
  }
  
  // calculate the background-subtracted pairspec data
  if(gROOT->FindObject("PairSpec_Open")){
    if((PairSpec_Open->Integral())>0) {
      stringstream cmd;
      cmd << "caput -a BEAM:PAIRSPEC:TaggEff 352";
      for(int n=1; n<=352; n++) {
	Double_t open = PairSpec_Open->GetBinContent(n);
	Double_t gated = PairSpec_Gated->GetBinContent(n);
	Double_t gated_dly = PairSpec_GatedDly->GetBinContent(n);
	Double_t taggeff = (gated-gated_dly)/open;
	if(!TMath::Finite(taggeff))
          taggeff = 0;
	cmd << taggeff << " ";
      }
      cmd << " > /dev/null";
      system(cmd.str().c_str());
    }
  }

  Double_t dError = 0;

  // look for hole in MWPC
  if(gROOT->FindObject("MWPC_Wires_Hits")){
    Int_t iNBins = MWPC_Wires_Hits->GetNbinsX();
    if((MWPC_Wires_Hits->Integral()) > (100*iNBins)){
      Int_t iPrev, iThis;
      Double_t dDiff;
      Int_t iProb = 0;
      iThis = MWPC_Wires_Hits->GetBinContent(1);
      for(Int_t i=1; i<iNBins; i++){
	iPrev = iThis;
	iThis = MWPC_Wires_Hits->GetBinContent(i+1);
	dDiff = (TMath::Abs((iThis-iPrev)/(1.*iPrev)));
	if(dDiff > 0.5) iProb++;
      }
      if(iProb > 16){
	printf("Possible problem in MWPC Wires - Event %d\n\t\t\t%d jumps found\n\n",gAN->GetNDAQEvent(),iProb);
	//dError += 1000;
      }
    }
    if((MWPC_Wires_Hits->Integral()) > (400*iNBins)) MWPC_Wires_Hits->Reset();
  }


  // look for shift in FPD
  if(gROOT->FindObject("FPD_TimeOR")){
    Int_t iNBins = FPD_TimeOR->GetNbinsX();
    if((FPD_TimeOR->Integral()) > (100*iNBins)){
      TH1D *Temp_FPD = (TH1D*)FPD_TimeOR->Clone("Temp_FPD");
      if(gROOT->FindObject("Prev_FPD")){
	Temp_FPD->Add(Prev_FPD,-1);
	if((Temp_FPD->Integral()) > (100*iNBins)){
	  delete Prev_FPD;
	  TH1D *Prev_FPD = (TH1D*)FPD_TimeOR->Clone("Prev_FPD");
	}
      }
      else TH1D *Prev_FPD = (TH1D*)FPD_TimeOR->Clone("Prev_FPD");
      
      Int_t iPeak = Temp_FPD->GetMaximumBin();
      Double_t dTime = Temp_FPD->GetBinCenter(iPeak);
      //Double_t dTime = Temp_FPD->GetMean();
      if(dTime < 5 || dTime > 25){
	printf("Possible problem in FPD_TimeOR - Event %d\n\t\t\tPeak at %f ns\n\n",gAN->GetNDAQEvent(),dTime);
	dError += 2000;
      }	  
      delete Temp_FPD;
    }
  }
  /*
  // look for shift in NaI
  Bool_t bShift = false;
  if(gROOT->FindObject("NaI_Hits_v_TimeOR")){
    Int_t iNBins = ((NaI_Hits_v_TimeOR->GetNbinsX())*(NaI_Hits_v_TimeOR->GetNbinsX()));
    if((NaI_Hits_v_TimeOR->Integral()) > (10000*720/8)){
      TH2D *Temp_NaI = (TH2D*)NaI_Hits_v_TimeOR->Clone("Temp_NaI");
      if(gROOT->FindObject("Prev_NaI")){
	Temp_NaI->Add(Prev_NaI,-1);  
	if((Temp_NaI->Integral()) > (10000*720/8)){
	  delete Prev_NaI;
	  TH2D *Prev_NaI = (TH2D*)NaI_Hits_v_TimeOR->Clone("Prev_NaI");
	}
      }
      else TH2D *Prev_NaI = (TH2D*)NaI_Hits_v_TimeOR->Clone("Prev_NaI");
	
      TH1D *Proj_NaI;
      for(Int_t i=0; i<720; i+=8){
	if((i==32) || (i==352) || (i==360) || (i==680)) continue;
	Proj_NaI = (TH1D*)Temp_NaI->ProjectionX("Proj_NaI",i+1,i+8);
	
	Double_t dTime = Proj_NaI->GetMean();
	if((dTime < 0) || (dTime > 40)){
	  if(!bShift) printf("Possible problem in NaI_Hits_v_TimeOR - Event %d\n",gAN->GetNDAQEvent());
	  printf("\t\t\tPeak at %f ns (Channels %3d-%3d)\n",dTime,i,i+7);
	  if(!bShift) dError += 4000;
	  bShift = true;
	}
      }
      if(gROOT->FindObject("Error_NaI")) delete Error_NaI;
      if(bShift){
	TH2D *Error_NaI = (TH2D*)Temp_NaI->Clone("Error_NaI");
	printf("\n");
      }
      delete Proj_NaI;
      delete Temp_NaI;
    }
  }
  */
  // look for shift in NaI
  Bool_t bShift = false;
  if(gROOT->FindObject("NaI_Hits_v_TimeOR")){
    Int_t iNBins = ((NaI_Hits_v_TimeOR->GetNbinsX())*(NaI_Hits_v_TimeOR->GetNbinsX()));
    if((NaI_Hits_v_TimeOR->Integral()) > (10000*720/8)){	
      TH1D *Proj_NaI;
      for(Int_t i=0; i<720; i+=8){
	if((i==32) || (i==352) || (i==360) || (i==680)) continue;
	Proj_NaI = (TH1D*)NaI_Hits_v_TimeOR->ProjectionX("Proj_NaI",i+1,i+8);
	
	Double_t dTime = Proj_NaI->GetMean();
	if((dTime < 30) || (dTime > 60)){
	  if(!bShift) printf("Possible problem in NaI_Hits_v_TimeOR - Event %d\n",gAN->GetNDAQEvent());
	  printf("\t\t\tPeak at %f ns (Channels %3d-%3d)\n",dTime,i,i+7);
	  if(!bShift) dError += 4000;
	  bShift = true;
	}
      }
      if(bShift) printf("\n");
      else NaI_Hits_v_TimeOR->Reset();
    }
  }

  // determine number of events with a hardware error
  if((gROOT->FindObject("NHardwareError")) && (dError == 0)){
    if((gAN->GetNDAQEvent()) < 3000) NHardwareError->Reset();
    dError = (NHardwareError->Integral(2,-1));
  }

  // now write the resulting errors to epics
  stringstream cmd;
  cmd << "caput GEN:MON:EventsWithError.A " << dError << " > /dev/null";
  system(cmd.str().c_str());

  // read out pi0 counts from 
  if(gROOT->FindObject("PHYS_IM_2g")){
    Double_t dNPi0 = 0;
    if(gROOT->FindObject("Prev_IM")){
      TH1D *Temp_IM = (TH1D*)PHYS_IM_2g->Clone("Temp_IM");
      Temp_IM->Add(Prev_IM,-1);
      if((Temp_IM->Integral()) > 100){
	dNPi0 = Temp_IM->Integral();
	  
	delete Temp_IM;
	delete Prev_IM;
	  
	TH1D *Prev_IM = (TH1D*)PHYS_IM_2g->Clone("Prev_IM");
      }
    }
    else{
      dNPi0 =PHYS_IM_2g->Integral();
      TH1D *Prev_IM = (TH1D*)PHYS_IM_2g->Clone("Prev_IM");
    }

    stringstream cmd;
    cmd << "caput GEN:MON:Pi0PerScRead.A " << dNPi0 << " > /dev/null";
    system(cmd.str().c_str());
  }
  
  if(bufferCounter>10){
    // calculate the background-subtracted pairspec data
    if(gROOT->FindObject("PairSpec_Open")){
      if((PairSpec_Open->Integral())>0) {
	stringstream cmd;
	cmd << "caput -a BEAM:PAIRSPEC:TaggEff 352";
	for(int n=1; n<=352; n++) {
	  Double_t open = PairSpec_Open->GetBinContent(n);
	  Double_t gated = PairSpec_Gated->GetBinContent(n);
	  Double_t gated_dly = PairSpec_GatedDly->GetBinContent(n);
	  Double_t taggeff = (gated-gated_dly)/open;
	  if(!TMath::Finite(taggeff))
	    taggeff = 0;
	  cmd << taggeff << " ";
	}
	cmd << " > /dev/null";
	system(cmd.str().c_str());
      }
    }
    cout << "Period Mac Buffer Couned " << endl;
    bufferCounter=0;
  }
  bufferCounter++;
  



}
