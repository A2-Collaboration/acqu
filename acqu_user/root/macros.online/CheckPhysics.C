void CheckPhysics()
{
	Char_t* hname[] = {  // pad number, histogram name, lin or log y axis, draw option
		"1", "Online_invMass2g", "lin", "",
		"2", "Online_invMass2g1p", "lin", "",
		"3", "Online_invMass2CB", "lin", "",
		"4", "Online_invMass2CB1TAPS", "lin", "",
		"5", "Online_invMass3g", "lin", "",
		"6", "Online_invMass3g1p", "lin", "",
		"7", "Online_invMass3CB", "lin", "",
		"8", "Online_invMass3CB1TAPS", "lin", "",
		"9", "Online_invMass6g", "lin", "",
		"10", "Online_invMass6g1p", "lin", "",
		"11", "Online_invMass6CB", "lin", "",
		"12", "Online_invMass6CB1TAPS", "lin", ""
	};
/*  Char_t* hname[] = {
    "1",  "EtaM_Nphoton",   "log",  "",
    "2",  "EtaM_Nproton",   "log",  "",
    //    "3",  "EtaM_Nneutron",  "log",  "",
    "3",  "EtaM_Npi0",      "log",  "",
    //    "5",  "EtaM_Npiplus",   "log",  "",
    "4",  "EtaM_Neta",      "log",  "",
    "5",  "EtaM_M2g",      "log",  "",
    "6",  "EtaM_M6g",      "lin",  "",
    "7",  "EtaM_EmEtaP",   "lin",  "",
    "7",  "EtaM_EmEtaR",   "lin",  "Same",
    "8",  "EtaM_EmEtapP",   "lin",  "",
    "8",  "EtaM_EmEtapR",   "lin",  "Same",
    //    "8",  "EtaM_EmComptonpP",   "lin",  "",
    //    "8",  "EtaM_EmComptonpR",   "lin",  "Same",
    "9",  "EtaM_M2gCBTAPS","lin",  "",
    "10", "EtaM_M2g_v_CB_NaI_ClNhitsOR",        "log",  "CONTZ",
    "11", "EtaM_M2g_v_TAGG_FPD_HitsPrompt",     "log",  "COLZ",
    "12", "EtaM_M6g_v_TAGG_FPD_HitsPrompt",       "log",  "COLZ",
  };*/

	Char_t* xname[] = {
		"invariant mass of 2 photons",
		"invariant mass of 2 photons + 1 proton",
		"invariant mass of 2 CB clusters",
		"invariant mass of 2 CB clusters + 1 TAPS cluster",
		"invariant mass of 3 photons",
		"invariant mass of 3 photons + 1 proton",
		"invariant mass of 3 CB clusters",
		"invariant mass of 3 CB clusters + 1 TAPS cluster",
		"invariant mass of 6 photons",
		"invariant mass of 6 photons + 1 proton",
		"invariant mass of 6 CB clusters",
		"invariant mass of 6 CB clusters + 1 TAPS cluster"
	};

	Int_t col[] = {2, 3, 4, 5, 2, 3, 4, 5, 2, 3, 4, 5};

	UInt_t i, j, pad_no;
	TH1F *h1;

	TCanvas* canv = new TCanvas("Physics-Spectra", "Physics-Online", 240, 180, 1240, 890);
	canv->Divide(4, 3);
	for (i = j = 0; j < 48; i++, j += 4) {
		sscanf(hname[j], "%d", &pad_no);  // pad_no is 1st argument
		canv->cd(pad_no);
		h1 = (TH1F*)gROOT->FindObjectAny(hname[j+1]);  // name is 2nd argument
		if (strcmp(hname[j+2], "log") == 0)  // lin/log is 3rd argument
			canv->GetPad(pad_no)->SetLogy();
		if (h1) {
			h1->Draw(hname[j+3]);  // parameter is 4th argument
			h1->GetXaxis()->SetTitle(xname[i]);
			h1->SetLineColor(1);
			h1->SetFillColor(col[i]);
		} else
			printf("No histogram named \"%s\" found.\nCheck macro!\n", hname[j+1]);
	}

  return;
}

PhysicsClear(){
  // Does not work properly. (Reset()!)
  //  TA2Apparatus* t = (TA2Apparatus*)gROOT->FindObject("EtaM");
  //  if(t)  t->ZeroAll();
  //  else printf("EtaM Physics not found\n");
  //  Scalers364to491->Reset();
  printf("Clear Physics spectra not implemented...try Clear Everything\n");
}
