
void InitToolbar(){
  // create a toolbar and connect macros to the buttons.
  // the nacros need to be loaded before.

  char *button_names[] =    {"Check CB",   "Clear CB",
			     "Check TAPS BaF2", "Clear TAPS BaF2",
			     "Check TAPS PWO", "Clear TAPS PWO",
			     "Check TAPS Veto", "Clear TAPS Veto",
			     "Check Tagger", "Clear Tagger",
			     "Check Microscope", "Clear Microscope",
			     // "Check EPT", "Clear EPT",
	      		     "Check LinPol",
			     "Check PID", "Clear PID",
			     "Check MWPC", "Clear MWPC",
			     "Check Trigger", "Clear Trigger",
			     "Check LiveTimes", "Clear LiveTimes",                        
			     "Check Physics", "Clear Physics",
              		     "Save All as PNG",
			     "Clear Everything",
			     0};
  
  char* button_actions[] =  {"CheckCB()", "CBClear()",
			     "CheckTAPSBaF2()","TAPSBaF2Clear()", 
			     "CheckTAPSPbWO4()", "TAPSPbWO4Clear()",
			     "CheckTAPSVeto()", "TAPSVetoClear()",
			     "CheckTagger2()", "TaggerClear()",
			     "CheckMicro()", "MicroClear()",
                             //"CheckEPT()", "EPTClear()",
	      		     "CheckLinPolDave()",
			     "CheckPID()", "PIDClear()",
			     "CheckMWPC()", "MWPCClear()",
			     "CheckTrigger()", "TriggerClear()",
                             "CheckLiveTimes()","LiveTimesClear()",                             
			     "CheckPhysics()", "PhysicsClear()",
		             "SaveOpenSpectras()",
			     "gAN->ZeroAll()",
			     0};
  char* button_comments[] = {"Representative CB spectra", 
			     "Clear CB spectra",
			     "Representative TAPS BaF2 spectra",
			     "Clear TAPS BaF2 spectra",
			     "Representative TAPS-PbWO4 spectra",
			     "Clear TAPS-PbWO4 spectra",
			     "Representative TAPS Veto-Time spectra",
			     "Clear TAPS Veto-Time spectra", 
			     "Representative Tagger spectra", 
			     "Clear Tagger spectra",
			     "Representative Microscope spectra", 
			     "Clear Microscope spectra",
                             //"Representative EPT spectra",
	    		     //"Clear EPT spectra",
	      		     "Representative coherent bremsstrahlung spectra", 
			     "Representative Particle ID (CB) spectra",
			     "Clear Particle ID spectra",
			     "Representative MWPC (CB) tracker spectra",
			     "Clear MWPC spectra",
			     "Representative trigger and CB-Tagger-TAPS synchronisation",
			     "Clear trigger spectra", 
                             "LiveTime spectra",
                             "Clear livetime spectra"
			     "Online reconstruction of detected particles",
			     "Clear particle reconstruction", 
			     "Save all open spectras",
			     "Clear ALL ONLINE Spectra", 
			     0};

  TControlBar* bar = new TControlBar("vertical", "Check CB-TAPS-Tagger Setup");
  UInt_t i = 0;
  char bar_title[50];
  char bar_action[50];
  char bar_comment[255];

  while (button_names[i]){
    sprintf(bar_title,   "%s", button_names[i]);
    sprintf(bar_action,  "%s", button_actions[i]);
    sprintf(bar_comment, "%s", button_comments[i]);
    bar->AddButton(bar_title, bar_action, bar_comment);
    i++;
  }
  bar->Show();

}

void InitPackages(){
  // load some macro packages
  // at the moment it is not checked it they really exist.
  char macro_path[] = "/home/wagners/acqu/acqu_user/root/macros.online/";
  char *macro_names[] =   {"CheckCB.C",
			   "CheckTAPSBaF2.C", 
			   "CheckTAPSPbWO4.C",
			   "CheckTAPSVeto.C",
			   "CheckTagger2.C",
			   "CheckMicro.C",
 			   //"CheckEPT.C",
			   "CheckLinPol.C",
			   "CheckPID.C",
			   "CheckMWPC.C",
			   "CheckTrigger.C",
                           "CheckLiveTimes.C",
			   "CheckPhysics.C",
			   "PrintThesis.C",
			   "SaveSpectras.C",
			    0};
  UInt_t i = 0;
  char cmd[255];
  while (macro_names[i]){
    strcpy(cmd, ".L ");
    strcat(cmd, macro_path);
    strcat(cmd, macro_names[i]);
    printf("%s\n", cmd);
    gROOT->ProcessLine(cmd);
    i++;
  }
}

CheckSystem(){
  InitPackages();
  InitToolbar();
 gStyle->SetPalette(1);
}
