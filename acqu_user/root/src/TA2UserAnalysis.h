//--Author	JRM Annand   27th Jan 2004...Apapt A.Starostin code
//--Rev         JRM Annand   12th May 2004...TA2Calorimeter etc. added
//--Rev         JRM Annand   21st Nov 2004...TA2CosmicCal added
//--Rev         JRM Annand   15th Jul 2005...TA2CrystalBall, TA2TAPS
//--Update      JRM Annand   15th Jul 2005...Use TA2Analysis::LoadVariable
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2UserAnalysis
//
// User analysis recognises user-written apparatus and physics classes
// New apparati and physics should be entered in Map_t kKnownChild
// and also in the switch statement of CreateChild() where the apparatus
// constructers are called

#ifndef __TA2UserAnalysis_h__
#define __TA2UserAnalysis_h__

#include "TA2Analysis.h"
#include <time.h>

class TA2UserAnalysis : public TA2Analysis
{
 protected:
  Char_t* fFileName; //Name of the file currently analysed
  Bool_t  fFileDone; //Flag, if analysed file is finished
 public:
  TA2UserAnalysis(const char*);
  virtual ~TA2UserAnalysis();
  virtual TA2DataManager* CreateChild(const char*, Int_t);
  virtual void ParseDisplay(char*);
  virtual void LoadVariable();
  virtual void InitSaveTree(Char_t*);
  virtual void ChangeTreeFile(Char_t*);
  virtual void CloseEvent();
  virtual void Periodic();
  virtual void EndFile();
  virtual void Finish();
  virtual void Decode();
  virtual void Process();
  virtual Char_t* GetFileName() const { return fFileName; };
  void ReadRunName(Char_t*);
  void Speed();
  void ZeroAll();
  void ZeroHist(Char_t*);
  void SaveAll();
  void SaveAll(Char_t*);
  void PrivateInit();

  ClassDef(TA2UserAnalysis,1)
};

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::ReadRunName(Char_t* RunName)
{
  Char_t* Ptr;

  //Get filename of analysed .dat or .rd0 file
  if(strlen(fFileName)) //Online, i.e. .dat files:
  {
    strcpy(RunName, fFileName);
  }
  else //Offline, i.e. .rd0 files
  {
    strcpy(RunName, gROOT->GetFile()->GetName());
    //Cut directory names
    Ptr = strrchr(RunName, '/');
    strcpy(RunName, Ptr+1);
  }
  //Cut file extension (.dat or .rd0)
  Ptr = strrchr(RunName, '.');
  *Ptr = '\0';
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::Speed()
{
  static Int_t iStart = time(NULL);
  static Int_t nEvents = 0;

  nEvents++;
  if((time(NULL) - iStart)>=60)
  {
    iStart = time(NULL);
    printf("Processed %d Events/sec\n", nEvents/60);
    nEvents = 0;
  }
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::LoadVariable( )
{
  //Input name - variable pointer associations for any subsequent
  //cut or histogram setup.
  //LoadVariable( "name", pointer-to-variable, type-spec );
  //NB scaler variable pointers need the preceeding &
  //   array variable pointers do not.
  //type-spec ED prefix for a Double_t variable
  //          EI prefix for an Int_t variable
  //type-spec SingleX for a single-valued variable
  //          MultiX  for a multi-valued variable

  //                            name     pointer          type-spec
  //  TA2DataManager::LoadVariable("Mmiss",  &fMmiss,         EDSingleX);

  TA2Analysis::LoadVariable();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::ParseDisplay(char* line)
{
  TA2Analysis::ParseDisplay(line);            //Standard parse
  //if(!fIsConfigPass) return;                //Was it a standard command?
  //TA2HistManager::ParseDisplay(line);
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::InitSaveTree(Char_t* fname)
{
  if(fPhysics)fPhysics->InitSaveTree();
  TA2Analysis::InitSaveTree(fname);
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::ChangeTreeFile(Char_t* fname)
{
  CloseEvent();        //Close tree & file
  InitSaveTree(fname); //New tree & file
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::CloseEvent()
{
  if(fPhysics) fPhysics->CloseEvent();
  TA2Analysis::CloseEvent();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::Periodic()
{
  if(fPhysics) fPhysics->Periodic();
  TA2Analysis::Periodic();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::EndFile()
{
  if(fPhysics) fPhysics->EndFile();
  TA2Analysis::EndFile();
  fFileDone = true;
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::Finish()
{
  if(IsEndOfFile()) EndFile();
  if(fPhysics) fPhysics->Finish();
  TA2Analysis::Finish();
  fFileDone = true;
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::Decode()
{
  //Set correct file name. This code tries to get the new file name as long as
  //it differs from the previous one, in order to get along with possible
  //race conditions in file name handling when changing files
  if(gAR->IsOnline() && fFileDone) //Previous file was finished and file name has really changed
  {
    Char_t* TempName = strrchr(gAR->GetFileName(), '/') + 1;
    //If new file name is different from previous one...
    if(strcmp(fFileName, TempName))
    {
      delete[] fFileName;
      fFileName = new Char_t[strlen(TempName) + 1];
      strcpy(fFileName, TempName); //...get new name...
      fFileDone = false;           //...and don't try again
    }
  }

  //Call standard decode
  TA2Analysis::Decode();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::Process()
{
  //Measure processed events per second
  // This was a most irritating feature JRMA
  // Speed();

  //Call standard decode
  TA2Analysis::Process();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::ZeroHist(Char_t* name)
{
  //Reset histogram, find histogram by name

  //Check there is something to display
  if(!fIsDisplay)
  {
    printf(" Sorry no histograms enabled for this %s\n", this->ClassName());
    return;
  }

  TH1* h = (TH1*)(gROOT->FindObject(name)); //All histograms inherit from TH1
  if(h)
    h->Reset();
  else
    //Histogram not found...print message
    printf("Histogram %s not found for class %s\n", name, this->ClassName());
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::ZeroAll()
{
  if(fPhysics) fPhysics->ZeroAll();
  TA2Analysis::ZeroAll();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::SaveAll(Char_t* name)
{
  TFile AcquHist(name, "RECREATE");
  gROOT->GetList()->Write();
  AcquHist.Close();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::SaveAll()
{
  TFile AcquHist("ARHistograms.root", "RECREATE");
  gROOT->GetList()->Write();
  AcquHist.Close();
}

//-----------------------------------------------------------------------------

inline void TA2UserAnalysis::PrivateInit()
{
  //Initialise file name handling
  fFileDone = true;          //File name should be updated in Decode()
  fFileName = new Char_t[1]; //Allocate some memory and...
  fFileName[0] = '\0';       //...create ('nonsense') file name
}

//-----------------------------------------------------------------------------

#endif


