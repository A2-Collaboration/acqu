// Ken Livingston 2nd Oct 2012
//macro to test the EPICS readout defined in an AcquDAQ setup file
//run like this:
//root [0] .L testEpics.C 
//root [1] testEpics("AcquRoot/data/EpicsEvent_0.dat")
//
Bool_t loaded=kFALSE;          //flag that these libs get loaded
Char_t fEventBuff[4000];       //event buffer

void* out = fEventBuff;   // output buffer pointer
Short_t index;
Short_t id;
Char_t EventName[50];;
Short_t nchan;
Float_t outbuff[200];
Int_t type;
Int_t nElem;

const char *epicsTypeName[] = {"epicsBYTE",   "epicsSTRING",    "epicsSHORT",   "epicsLONG",   "epicsFLOAT",   "epicsDOUBLE",      NULL};


void loadLibs(){
  gSystem->Load( "libGui.so" );
  gSystem->Load( "/opt/acqu/build/lib/libezcaRoot.so" );
  gSystem->Load( "/opt/acqu/build/lib/libAcquRoot.so" );
  gSystem->Load( "/opt/acqu/build/lib/libAcquDAQ.so" );  
  loaded=kTRUE;
}

void testEpics(char *epicsSetupFile, char* outfile=NULL, int singleOnly=0){
  //epicsSetupFile is the name of the setup file for the EPICS module to be loaded in AcquDAQ
  
  Char_t* name="epicsTest"; //name for the module
  FILE* log=NULL;           //no logging here
  Char_t line[120];         //line for AcquDAQ to set up the module
  void *out = fEventBuff;   //pointer to the event buffer

  if(!loaded)loadLibs();    //load the libs if required

  TEPICSmodule *epics;      //Epics module
 
  fprintf(stderr, "\n\n------- TESTING SINGE READOUT SHOT-----------\n");
  
  sprintf(line,"dummy dummy dummy timer 0");                             //create the setup line 0 period
  epics = new TEPICSmodule(name,epicsSetupFile,log,line);                //create a module
  //  epics->SetVerbose(kTRUE);
  epics->FileConfig(epicsSetupFile);                                     //configure with the setup file
  epics->PostInit();                                                     //do the post init
  epics->WriteEPICS(&out);                                               //do the actual readout
  epics->DumpBuffer(fEventBuff,outfile);                                 //Check that the buffer is good
 
  if(singleOnly){
    delete epics;
    return;
  }


  fprintf(stderr, "\n\n------- TESTING TIMED PERIODIC READOUT MODE-----------\n");
  epics->Start();                                       //start the timer
  for(int n=0;n<2;n++){                                 //do 2 events
    fprintf(stderr,"Event %d\n",n);
    if(epics->IsTimedOut()){                            //check to see if timed out
      out = fEventBuff;                                 //retes out every time since WriteEPICS(&out) moves it along the buffer
      epics->WriteEPICS(&out);                          //do the actual readout
      epics->DumpBuffer(fEventBuff,outfile);            //Check that the buffer is good
      fprintf(stderr,"----- Read %d bytes from the EPICS buffer ----\n\n",out-fEventBuff); //count the bytes of the ptr in the buffer.  
       epics->Start();                                   //restart the timer   
    }
    gSystem->Sleep(500);
  }
  delete epics;
  
  
  sprintf(line,"dummy dummy dummy scaler 2");                            //create the setup line mode scaler  period = 2 scaler events
  epics = new TEPICSmodule(name,epicsSetupFile,log,line);                //create a module
  //  epics->SetVerbose(kTRUE);
  epics->FileConfig(epicsSetupFile);                                     //configure with the setup file
  epics->PostInit();                                                     //do the post init

  fprintf(stderr, "\n\n------- #SCALER EVENTS PERIODIC READOUT MODE-----------\n");
  epics->Start();                                       //start the counter
  for(int n=0;n<5;n++){                                 //do 5 events
    fprintf(stderr,"Scaler Event %d\n",n);
    epics->Count();                                     //force a dummy scaler count
    if(epics->IsCountedOut()){                          //check to see if count is reached
      out = fEventBuff;                                 //reset out every time since WriteEPICS(&out) moves it along the buffer
      epics->WriteEPICS(&out);                          //do the actual readout
      epics->DumpBuffer(fEventBuff,outfile);            //Check that the buffer is good
      fprintf(stderr,"----- Read %d bytes from the EPICS buffer ----\n\n",out-fEventBuff); //count the bytes of the ptr in the buffer.  
      epics->Start();                                   //restart the counter 
    }
    gSystem->Sleep(100);
  }
  delete epics;
}

void getBufferInfo(char *epicsSetupFile){
  Char_t* name="epicsTest"; //name for the module
  FILE* log=NULL;           //no logging here
  Char_t line[120];         //line for AcquDAQ to set up the module
  void *out = fEventBuff;   //pointer to the event buffer

  if(!loaded)loadLibs();    //load the libs if required

  TEPICSmodule *epics;      //Epics module


  sprintf(line,"dummy dummy dummy timer 370");                           //create the setup line mode timer 0.37s period
  epics = new TEPICSmodule(name,epicsSetupFile,log,line);                //create a module
  //  epics->SetVerbose(kTRUE);
  epics->FileConfig(epicsSetupFile);                                     //configure with the setup file
  epics->PostInit();                                                     //do the post init
  epics->WriteEPICS(&out);                                               //do the actual readout

  epics->GetBufferInfo(&index, &id, EventName, &nchan,fEventBuff, &timestamp);
  fprintf(stderr, "index: %hd, id: %hd, EventName: %s, nchan: %hd, time: %s\n",index,id,EventName,nchan,asctime(localtime(&timestamp)));
  delete epics;

}

void getChannel(char *epicsSetupFile,Char_t *PVName){
  Char_t* name="epicsTest"; //name for the module
  FILE* log=NULL;           //no logging here
  Char_t line[120];         //line for AcquDAQ to set up the module
  void *out = fEventBuff;   //pointer to the event buffer

  if(!loaded)loadLibs();    //load the libs if required

  TEPICSmodule *epics;      //Epics module
  sprintf(line,"dummy dummy dummy timer 370");                           //create the setup line mode timer 0.37s period

  epics = new TEPICSmodule(name,epicsSetupFile,log,line);                //create a module
  epics->FileConfig(epicsSetupFile);                                     //configure with the setup file
  epics->PostInit();                                                     //do the post init
  epics->WriteEPICS(&out);                                               //do the actual readout

  epics->GetChannel(PVName, &type, fEventBuff, outbuff, &nElem);
  
  fprintf(stdout,"\nPV: %s, type: %s, nElement: %d \n", PVName,epicsTypeName[type],nElem);
  
  for(int n=0;n<nElem;n++){
    fprintf(stdout,"%s[%d] = %f\n",PVName,n,outbuff[n]);
  }
  delete epics;
}
