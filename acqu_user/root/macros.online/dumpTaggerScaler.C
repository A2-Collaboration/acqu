//Dump the current FP scalers to a nfs exported file system

void dumpAcc(){
  FILE  *fp=fopen("/exports/a2raid5/scaler_acc.dat","w");
  for(int n=1;n<=352;n++){
    //    fprintf(fp,"%d\n", FPD_ScalerCurr->GetBinContent(353-n));
    //   fprintf(fp,"%d\n", FPD_ScalerCurr->GetBinContent(n));
    fprintf(fp,"%d\n", FPD_ScalerAcc->GetBinContent(n));
  }
  fclose(fp);
}
