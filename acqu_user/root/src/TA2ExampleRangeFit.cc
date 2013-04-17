#include <iostream>
#include <fstream>

// AcquROOT includes
#include "TA2UserAnalysis.h"
#include "TA2CentralApparatus.h"
#include "TA2CalArray.h"
#include "TA2PlasticPID.h"
#include "TA2CylMwpc.h"
#include "TA2CentralTrack.h"
#include "TA2ExampleRangeFit.h"
#include "TA2Tagger.h"
#include "TA2RangeFit.h"

using namespace std;

static const Double_t mproton  = 0.938272013; // GeV
static const Double_t mneutron = 0.939565346; // GeV
static const Double_t mpi0     = 0.1349766;   // GeV
static const Double_t mpi      = 0.139570;    // GeV
static const Double_t RadToDeg = TMath::RadToDeg();
static const Double_t DegToRad = TMath::DegToRad();

ClassImp(TA2ExampleRangeFit)


// ================================================================
  TA2ExampleRangeFit::TA2ExampleRangeFit( const char* name, TA2Analysis* analysis ) : TA2Physics( name, analysis )
{
  // Initialise Physics variables here
  // Default null pointers, zeroed variables
  
  // Apparatuses & Detectrors
  fCB   = NULL;
  fNaI  = NULL;
  fMwpc = NULL;
  fPid  = NULL;
  fNintersTrue = NULL;
  fFPDTimeOR   = new Double_t[1056];
  fFPDHits     = new Int_t[1056];
  
  // Neutral
  fThetaNe = NULL;
  fPhiNe   = NULL;
  fEne     = NULL;
  fMclNe   = NULL;
  fTne     = NULL;
  fCentralIndexNaINe = NULL;

  // Charged
  fThetaCh   = NULL;
  fThetaChCB = NULL;
  fPhiCh     = NULL;
  fEch       = NULL;
  fMclCh     = NULL;
  fTch       = NULL;
  fEpid      = NULL;
  fIhitPid   = NULL;
  fPsVertexX = NULL;
  fPsVertexY = NULL;
  fPsVertexZ = NULL;
  fPsVertexGoodZ = NULL;
  fChTrackType   = NULL;
  fCentralIndexNaICh = NULL;

  for(Int_t e=0;e<24;e++) {
    fPID_dE[e]=NULL;
    fPID_dE_prot[e] = NULL;
    fPID_dE_pion[e] = NULL;
    fPID_dE_raw[e] = NULL;
  }
  fEmwpc    = NULL ;
  fDEmwpc   = NULL ;
  EventTime = NULL; 
  fTotal_Energy      = NULL;
  fTotal_Energy_prot = NULL;
  fTotal_Energy_pion = NULL;
}


// ================================================================
TA2ExampleRangeFit::~TA2ExampleRangeFit()
{
  // Free up allocated memory...after checking its allocated
  // detector and cuts lists

  // Neutral
  delete [] fThetaNe;
  delete [] fPhiNe;
  delete [] fEne;
  delete [] fMclNe;
  delete [] fTne;
  delete [] fCentralIndexNaINe;

  // Charged
  delete [] fThetaCh;
  delete [] fThetaChCB;
  delete [] fPhiCh;
  delete [] fEch;
  delete [] fMclCh;
  delete [] fTch;
  delete [] fCentralIndexNaICh;
  delete [] fEpid;
  delete [] fIhitPid;
  delete [] fPsVertexX;
  delete [] fPsVertexY;
  delete [] fPsVertexZ;
  delete [] fPsVertexGoodZ;
  delete [] fChTrackType;
  delete [] fDE_MWPC;
  delete [] EventTime;
  for(Int_t e=0;e<24;e++) {
    delete [] fPID_dE[e];
    delete [] fPID_dE_prot[e];
    delete [] fPID_dE_pion[e];
    delete [] fPID_dE_raw[e];
  }
  delete [] fTotal_Energy;
  delete [] fTotal_Energy_prot;
  delete [] fTotal_Energy_pion;
}


// ================================================================
void TA2ExampleRangeFit::PostInit()
{
  // Initialise arrays to contain 4 momenta and plotable scaler variables
  // Missing mass, missing energy, cm momentum, energies, angles
  // Initialisation will abort if CB or Tagger not initialised
  // TAPS is optional

  nBeam = 1000; // Tagger

  // CB
  fCB = (TA2CentralApparatus*)((TA2Analysis*)fParent)->GetChild("CB");
  if (!fCB) PrintError("","<No Central Apparatus class found in annalysis>",EErrFatal);

  // NaI
  fNaI = (TA2CalArray*)((TA2Analysis*)fParent)->GetGrandChild("NaI");
  if (!fNaI) PrintError("Warning!","<No NaI class found in annalysis>");
  
  // Mwpc
  fMwpc = (TA2CylMwpc*)((TA2Analysis*)fParent)->GetGrandChild("CylMWPC");
  if (!fMwpc)
    PrintError("Warning!","<No Mwpc class found in annalysis>");
  else
    fNintersTrue = fMwpc->GetNintersTrue();
      
  // Pid
  fPid = (TA2PlasticPID*)((TA2Analysis*)fParent)->GetGrandChild("PID");
  if (!fPid) PrintError("Warning!","<No Pid class found in annalysis>");

  // TAGGER
  fTagger = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
  if (!fTagger)
    PrintError("","<No TAGGER class found in analysis>",EErrFatal);
  else {
    fTAGGp4 = fTagger->GetP4();
    // FPD
    fFPD = (TA2Ladder*)fTagger->GetChild("FPD", "TA2Detector");
    if (!fFPD)
      PrintError("","<No FPD class found in analysis>",EErrFatal);
  }

  fRangeFit = new TA2RangeFit("RangeFit", (TA2Analysis*)fParent);               // gets pointer to RangeFit
  if(!fRangeFit) 
    PrintError("","<No RangeFit class found in analysis>",EErrFatal);

  // Neutral
  Int_t nMaxTracks = fCB->GetMaxTrack()+1;
  fThetaNe = new Double_t[nMaxTracks];
  fPhiNe   = new Double_t[nMaxTracks];
  fEne     = new Double_t[nMaxTracks];
  fMclNe   = new Double_t[nMaxTracks];
  fTne     = new Double_t[nMaxTracks];
  fCentralIndexNaINe = new Double_t[nMaxTracks];
  // Charged
  fThetaCh   = new Double_t[nMaxTracks];
  fThetaChCB = new Double_t[nMaxTracks];
  fPhiCh     = new Double_t[nMaxTracks];
  fEch       = new Double_t[nMaxTracks];
  fMclCh     = new Double_t[nMaxTracks];
  fTch       = new Double_t[nMaxTracks];
  fEpid      = new Double_t[nMaxTracks];
  fIhitPid   = new Double_t[nMaxTracks];
  fPsVertexX = new Double_t[nMaxTracks];
  fPsVertexY = new Double_t[nMaxTracks];
  fPsVertexZ = new Double_t[nMaxTracks];
  fPsVertexGoodZ = new Double_t[nMaxTracks];
  fChTrackType   = new Double_t[nMaxTracks];
  fCentralIndexNaICh = new Double_t[nMaxTracks];

  for(Int_t e=0;e<24;e++) {
    fPID_dE[e]= new Double_t[nMaxTracks];
    fPID_dE_prot[e] = new Double_t[nMaxTracks];
    fPID_dE_pion[e] = new Double_t[nMaxTracks];
    fPID_dE_raw[e]  = new Double_t[nMaxTracks];
  }
  EventTime = new Double_t[1000];//[nBeam];        // Timing spectrum
  // These spectra are created for each defined time window (prompt/random)
  for(Int_t w=1; w<OPT_NWIND+1; w++) {
    BeamPol10[w]    = new Int_t[nBeam];
    Egamma[w]       = new Double_t[1000];
    PhiPi[w]        = new Double_t[1000];
    ThetaPi[w]      = new Double_t[1000];
    PhiPi0[w]       = new Double_t[1000];
    ThetaPi0[w]     = new Double_t[1000];
    ThetaPi0Kin[w]  = new Double_t[1000];
    ThetaNeuKin[w]  = new Double_t[1000];
    DiffThetaPi0[w] = new Double_t[1000];
    DiffThetaNeu[w] = new Double_t[1000];
    EnePi[w]        = new Double_t[1000];
    EnePiKin[w]     = new Double_t[1000];
    EnePi0[w]       = new Double_t[1000];
    EnePi0Kin[w]    = new Double_t[1000];
    EneNeuKin[w]    = new Double_t[1000];
    PhiProton[w]    = new Double_t[1000];
    ThetaProton[w]  = new Double_t[1000];
    EneProton[w]    = new Double_t[1000];
    EneProtonKin[w] = new Double_t[1000];
    deltaphi[w]     = new Double_t[1000];
    KinProton[w]    = new Double_t[1000];
    KinPi0[w]       = new Double_t[1000];
    KinNeu[w]       = new Double_t[1000];
    KinPi[w]        = new Double_t[1000];
    TrfBadChi2[w]   = new Double_t[1000];
    TkineBadChi2[w] = new Double_t[1000];
    TrfOkChi2[w]    = new Double_t[1000];
    TkineOkChi2[w]  = new Double_t[1000];
    TrfBadChi2Pion[w]   = new Double_t[1000];
    TkineBadChi2Pion[w] = new Double_t[1000];
    TrfOkChi2Pion[w]    = new Double_t[1000];
    TkineOkChi2Pion[w]  = new Double_t[1000];
    EneCBProt[w]      = new Double_t[1000];
    CheckEneCBProt[w] = new Double_t[1000];
    CheckEneCBPion[w] = new Double_t[1000];
    EneCBKinProt[w]   = new Double_t[1000];
    EnePIDProt[w]     = new Double_t[1000];
    EnePIDKinProt[w]  = new Double_t[1000];
    EnePIDNormProt[w] = new Double_t[1000];
    EneCBPion[w]      = new Double_t[1000];
    EneCBKinPion[w]   = new Double_t[1000];
    EnePIDPion[w]     = new Double_t[1000];
    EnePIDKinPion[w]  = new Double_t[1000];
    EnePIDNormPion[w] = new Double_t[1000];
    MissingEneProt[w] = new Double_t[1000];
    MissingEnePion[w] = new Double_t[1000];
    CheckMEProt[w]    = new Double_t[1000];
    CheckMEPion[w]    = new Double_t[1000];
    EneCBStopPion[w]  = new Double_t[1000];
    EnePIDStopPion[w] = new Double_t[1000];
    CheckEnePIDProt[w] = new Double_t[1000];
    CheckEnePIDPion[w] = new Double_t[1000];
    MissingEnePIDProt[w] = new Double_t[1000];
    MissingEnePIDPion[w] = new Double_t[1000];
    EneCBKinStopPion[w]  = new Double_t[1000];
    EneCBKinLowPion[w]   = new Double_t[1000];
    EnePIDKinStopPion[w] = new Double_t[1000];
    CheckEneCBStopPion[w]= new Double_t[1000];
    OkEnePIDNormProt[w] = new Double_t[1000];
    OkEneCBProt[w] = new Double_t[1000];
    OkEnePIDNormPion[w] = new Double_t[1000];
    OkEneCBPion[w] = new Double_t[1000];
  }
  fEmwpc  = new Double_t[nMaxTracks];
  fDEmwpc = new Double_t[nMaxTracks];
  fTotal_Energy = new Double_t[nMaxTracks];
  fTotal_Energy_prot = new Double_t[nMaxTracks];
  fTotal_Energy_pion = new Double_t[nMaxTracks];
  
  // Default physics initialisation
  TA2Physics::PostInit();
}


// ================================================================
void TA2ExampleRangeFit::LoadVariable()
{
  // Input name - variable pointer associations for any subsequent
  // Default physics initialisation cut or histogram setup
  // LoadVariable( "name", pointer-to-variable, type-spec );

  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  Char_t Name[255];
  Char_t* VarName;

  TA2Physics::LoadVariable();
  TA2DataManager::LoadVariable("Ntr", &fNtracks, EISingleX);

  // Neutral
  TA2DataManager::LoadVariable("Nne",		    &fNne, 	        EISingleX);
  TA2DataManager::LoadVariable("ThetaNe",	    fThetaNe,	        EDMultiX);
  TA2DataManager::LoadVariable("PhiNe",		    fPhiNe,             EDMultiX);
  TA2DataManager::LoadVariable("Ene",               fEne,               EDMultiX);
  TA2DataManager::LoadVariable("MclNe",             fMclNe,             EDMultiX);
  TA2DataManager::LoadVariable("Tne",               fTne,               EDMultiX);
  TA2DataManager::LoadVariable("M2g",               &fM2g,              EDSingleX);
  TA2DataManager::LoadVariable("CentralIndexNaINe", fCentralIndexNaINe, EDMultiX);

  // Charged
  TA2DataManager::LoadVariable("Nch",           &fNch,          EISingleX);
  TA2DataManager::LoadVariable("ThetaCh",       fThetaCh,       EDMultiX);
  TA2DataManager::LoadVariable("PhiCh",         fPhiCh,         EDMultiX);
  TA2DataManager::LoadVariable("Ech",           fEch,           EDMultiX);
  TA2DataManager::LoadVariable("MclCh",         fMclCh,         EDMultiX);
  TA2DataManager::LoadVariable("Tch",           fTch,           EDMultiX);
  TA2DataManager::LoadVariable("Epid",          fEpid,          EDMultiX);
  TA2DataManager::LoadVariable("IhitPid",       fIhitPid  ,     EDMultiX);
  TA2DataManager::LoadVariable("PsVertexX",     fPsVertexX,     EDMultiX);
  TA2DataManager::LoadVariable("PsVertexY",     fPsVertexY,     EDMultiX);
  TA2DataManager::LoadVariable("PsVertexZ",     fPsVertexZ,     EDMultiX);
  TA2DataManager::LoadVariable("ChTrackType",   fChTrackType,   EDMultiX);
  TA2DataManager::LoadVariable("PsVertexGoodZ", fPsVertexGoodZ, EDMultiX);
  TA2DataManager::LoadVariable("EventTime",     EventTime,      EDMultiX);

  TA2DataManager::LoadVariable("EnergyCB",      &fECB,          EDSingleX);
  TA2DataManager::LoadVariable("EnergyPID",     &fEPID,         EDSingleX);

  // plots for range fit
  TA2DataManager::LoadVariable("EneRFProt",     &EneRFProt,     EDSingleX);
  TA2DataManager::LoadVariable("EneRFPion",     &EneRFPion,     EDSingleX);
  TA2DataManager::LoadVariable("chi2Prot",      &chi2Prot,      EDSingleX);
  TA2DataManager::LoadVariable("chi2ProtPion",  &chi2ProtPion,  EDSingleX);
  TA2DataManager::LoadVariable("chi2Pion",      &chi2Pion,      EDSingleX);
  TA2DataManager::LoadVariable("chi2PionProt",  &chi2PionProt,  EDSingleX);
  TA2DataManager::LoadVariable("CumulProt",     &CumulProt,     EDSingleX);
  TA2DataManager::LoadVariable("CumulPion",     &CumulPion,     EDSingleX);
  TA2DataManager::LoadVariable("CumulAllProt",  &CumulAllProt,  EDSingleX);
  TA2DataManager::LoadVariable("CumulAllPion",  &CumulAllPion,  EDSingleX);
  TA2DataManager::LoadVariable("RanErrProt",    &RanErrProt,    EISingleX);
  TA2DataManager::LoadVariable("RanErrPion",    &RanErrPion,    EISingleX);
  TA2DataManager::LoadVariable("CentralIndexNaICh",fCentralIndexNaICh,EDMultiX);


  // histos for bananas (24 PID sectors)
  for(Int_t e=0;e<24;e++) {
    sprintf(Name, "PID_dE%d", e);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, fPID_dE[e], EDMultiX);

    sprintf(Name, "PID_dE_prot%d", e);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, fPID_dE_prot[e], EDMultiX);

    sprintf(Name, "PID_dE_pion%d", e);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, fPID_dE_pion[e], EDMultiX);

    sprintf(Name, "PID_dE_raw%d", e);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, fPID_dE_raw[e], EDMultiX);
  }

  for(Int_t w=1; w<OPT_NWIND+1; w++) {
    sprintf(Name, "BeamPol10W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, BeamPol10[w], EIMultiX);

    sprintf(Name, "EgammaW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, Egamma[w], EDMultiX);  
			 
    sprintf(Name, "PhiPiW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, PhiPi[w], EDMultiX);

    sprintf(Name, "ThetaPiW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, ThetaPi[w], EDMultiX);

    sprintf(Name, "PhiPi0W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, PhiPi0[w], EDMultiX);

    sprintf(Name, "ThetaPi0W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, ThetaPi0[w], EDMultiX);
    
    sprintf(Name, "ThetaNeuKinW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, ThetaNeuKin[w], EDMultiX);
    
    sprintf(Name, "ThetaPi0KinW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, ThetaPi0Kin[w], EDMultiX);

    sprintf(Name, "DiffThetaNeuW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, DiffThetaNeu[w], EDMultiX);

    sprintf(Name, "EnePiW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePi[w], EDMultiX);

    sprintf(Name, "EnePiKinW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePiKin[w], EDMultiX);

    sprintf(Name, "EnePi0W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePi0[w], EDMultiX);

    sprintf(Name, "EnePi0KinW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePi0Kin[w], EDMultiX);

    sprintf(Name, "EneNeuKinW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneNeuKin[w], EDMultiX);

    sprintf(Name, "ThetaProtonW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, ThetaProton[w], EDMultiX);

    sprintf(Name, "PhiProtonW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, PhiProton[w], EDMultiX);

    sprintf(Name, "EneProtonW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneProton[w], EDMultiX);

    sprintf(Name, "EneProtonKinW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneProtonKin[w], EDMultiX);

    sprintf(Name, "MissingEneProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, MissingEneProt[w], EDMultiX);

    sprintf(Name, "MissingEnePionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, MissingEnePion[w], EDMultiX);

    sprintf(Name, "MissingEnePIDProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, MissingEnePIDProt[w], EDMultiX);

    sprintf(Name, "MissingEnePIDPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, MissingEnePIDPion[w], EDMultiX);

    sprintf(Name, "CheckMEProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckMEProt[w], EDMultiX);

    sprintf(Name, "CheckMEPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckMEPion[w], EDMultiX);

    sprintf(Name, "KinProtonW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, KinProton[w], EDMultiX);

    sprintf(Name, "KinPiW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, KinPi[w], EDMultiX);

    sprintf(Name, "KinPi0W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, KinPi0[w], EDMultiX);

    sprintf(Name, "KinNeuW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, KinNeu[w], EDMultiX);

    sprintf(Name, "TrfBadChi2W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TrfBadChi2[w], EDMultiX);

    sprintf(Name, "TrfOkChi2W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TrfOkChi2[w], EDMultiX);

    sprintf(Name, "TkineBadChi2W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TkineBadChi2[w], EDMultiX);

    sprintf(Name, "TkineOkChi2W%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TkineOkChi2[w], EDMultiX);

    sprintf(Name, "TrfBadChi2PionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TrfBadChi2Pion[w], EDMultiX);

    sprintf(Name, "TrfOkChi2PionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TrfOkChi2Pion[w], EDMultiX);

    sprintf(Name, "TkineBadChi2PionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TkineBadChi2Pion[w], EDMultiX);

    sprintf(Name, "TkineOkChi2PionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, TkineOkChi2Pion[w], EDMultiX);

    sprintf(Name, "deltaphiW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, deltaphi[w], EDMultiX);

    sprintf(Name, "EneCBProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBProt[w], EDMultiX);

    sprintf(Name, "CheckEneCBProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckEneCBProt[w], EDMultiX);

    sprintf(Name, "CheckEneCBPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckEneCBPion[w], EDMultiX);

    sprintf(Name, "CheckEnePIDProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckEnePIDProt[w], EDMultiX);

    sprintf(Name, "CheckEnePIDPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckEnePIDPion[w], EDMultiX);

    sprintf(Name, "EneCBKinProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBKinProt[w], EDMultiX);

    sprintf(Name, "EnePIDProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDProt[w], EDMultiX);

    sprintf(Name, "EnePIDKinProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDKinProt[w], EDMultiX);

    sprintf(Name, "EnePIDNormProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDNormProt[w], EDMultiX);

    sprintf(Name, "EneCBPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBPion[w], EDMultiX);

    sprintf(Name, "EneCBKinPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBKinPion[w], EDMultiX);

    sprintf(Name, "EneCBStopPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBStopPion[w], EDMultiX);

    sprintf(Name, "EneCBKinStopPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBKinStopPion[w], EDMultiX);

    sprintf(Name, "EneCBKinLowPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EneCBKinLowPion[w], EDMultiX);

    sprintf(Name, "EnePIDPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDPion[w], EDMultiX);

    sprintf(Name, "EnePIDKinPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDKinPion[w], EDMultiX);

    sprintf(Name, "EnePIDStopPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDStopPion[w], EDMultiX);

    sprintf(Name, "EnePIDKinStopPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDKinStopPion[w], EDMultiX);

    sprintf(Name, "EnePIDNormPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, EnePIDNormPion[w], EDMultiX);

    sprintf(Name, "CheckEneCBStopPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, CheckEneCBStopPion[w], EDMultiX);

    sprintf(Name, "OkEnePIDNormProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, OkEnePIDNormProt[w], EDMultiX);

    sprintf(Name, "OkEneCBProtW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, OkEneCBProt[w], EDMultiX);

    sprintf(Name, "OkEnePIDNormPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, OkEnePIDNormPion[w], EDMultiX);

    sprintf(Name, "OkEneCBPionW%d", w);
    VarName = new Char_t[strlen(Name)+1];
    strcpy(VarName,Name);
    TA2DataManager::LoadVariable(VarName, OkEneCBPion[w], EDMultiX);
  }

  TA2DataManager::LoadVariable("E_MWPC", fEmwpc, EDMultiX);
  TA2DataManager::LoadVariable("DE_MWPC", fDEmwpc, EDMultiX);
  TA2DataManager::LoadVariable("Total_Energy", fTotal_Energy, EDMultiX);
  TA2DataManager::LoadVariable("Total_Energy_prot", fTotal_Energy_prot, EDMultiX);
  TA2DataManager::LoadVariable("Total_Energy_pion", fTotal_Energy_pion, EDMultiX);

  Read_Cal();
}


// ================================================================
void TA2ExampleRangeFit::Reconstruct()
{
  evtDaq = gAN->GetNDAQEvent();
  
  //   cout << endl;
  //       cout << "Event # " << evtDaq << endl;
  //   fMwpc->Example( kTRUE);
  //   if (fNch==1)
  //     fCB->Example(kTRUE);
  
  Int_t iWindow;
  Double_t Tvirtual;

  EneRFProt  = -1;
  EneRFPion  = -1;
  chi2Prot   = -1;
  chi2Pion   = -1;
  chi2ProtPion= -1;
  chi2PionProt = -1;
  RanErrProt = -1;
  RanErrPion = -1;
  CumulProt  = -1;
  CumulPion  = -1;
  CumulAllProt = -1;
  CumulAllPion = -1;

  // proton bananas 
  Double_t x[24][8];
  Double_t y[24][8];
 
  for (Int_t i=0; i<24 ;i++) {
    x[i][0] = 0;
    x[i][1] = 100;
    x[i][2] = 250;
    x[i][3] = 400;

    x[i][4] = 400;
    x[i][5] = 250;
    x[i][6] = 100;
    x[i][7] = 0; 

    y[i][4] = 10;
    y[i][5] = 10;
    y[i][6] = 10;
    y[i][7] = 10;
  }
  
  y[0][0]  = 2.00;   y[0][1]  = 0.80;   y[0][2]  = 0.65;   y[0][3]  = 0.55;
  y[1][0]  = 1.90;   y[1][1]  = 0.80;   y[1][2]  = 0.65;   y[1][3]  = 0.55;
  y[2][0]  = 2.10;   y[2][1]  = 1.00;   y[2][2]  = 0.70;   y[2][3]  = 0.65;
  y[3][0]  = 1.75;   y[3][1]  = 0.90;   y[3][2]  = 0.70;   y[3][3]  = 0.65;  
  y[4][0]  = 1.75;   y[4][1]  = 1.00;   y[4][2]  = 0.85;   y[4][3]  = 0.75;
  y[5][0]  = 1.40;   y[5][1]  = 0.80;   y[5][2]  = 0.60;   y[5][3]  = 0.55;
  y[6][0]  = 1.80;   y[6][1]  = 0.85;   y[6][2]  = 0.60;   y[6][3]  = 0.60;
  y[7][0]  = 1.80;   y[7][1]  = 0.80;   y[7][2]  = 0.60;   y[7][3]  = 0.55;
  y[8][0]  = 1.85;   y[8][1]  = 0.95;   y[8][2]  = 0.55;   y[8][3]  = 0.60;
  y[9][0]  = 3.20;   y[9][1]  = 1.20;   y[9][2]  = 0.75;   y[9][3]  = 0.70;
  y[10][0] = 2.20;   y[10][1] = 1.00;   y[10][2] = 0.75;   y[10][3] = 0.70;
  y[11][0] = 2.20;   y[11][1] = 0.90;   y[11][2] = 0.70;   y[11][3] = 0.65;
  y[12][0] = 2.40;   y[12][1] = 0.95;   y[12][2] = 0.70;   y[12][3] = 0.65;
  y[13][0] = 1.50;   y[13][1] = 0.75;   y[13][2] = 0.55;   y[13][3] = 0.50;
  y[14][0] = 1.50;   y[14][1] = 0.75;   y[14][2] = 0.55;   y[14][3] = 0.50;
  y[15][0] = 2.50;   y[15][1] = 1.00;   y[15][2] = 0.65;   y[15][3] = 0.60;
  y[16][0] = 3.10;   y[16][1] = 1.40;   y[16][2] = 0.90;   y[16][3] = 0.90;
  y[17][0] = 2.10;   y[17][1] = 1.10;   y[17][2] = 0.75;   y[17][3] = 0.75;
  y[18][0] = 2.10;   y[18][1] = 0.90;   y[18][2] = 0.60;   y[18][3] = 0.55;
  y[19][0] = 2.10;   y[19][1] = 0.90;   y[19][2] = 0.60;   y[19][3] = 0.55;
  y[20][0] = 2.10;   y[20][1] = 1.00;   y[20][2] = 0.65;   y[20][3] = 0.65; 
  y[21][0] = 2.90;   y[21][1] = 1.10;   y[21][2] = 0.70;   y[21][3] = 0.70;
  y[22][0] = 2.70;   y[22][1] = 1.00;   y[22][2] = 0.65;   y[22][3] = 0.65;
  y[23][0] = 2.65;   y[23][1] = 1.00;   y[23][2] = 0.65;   y[23][3] = 0.65;

  // Speed measurement
  Speed();

  ResetEvent();

  // Get number of tracks in CB
  fNtracks = fCB->GetNtracks();
  fNne = fCB->GetNne(); // # of "neutral" tracks in CB
  fNch = fCB->GetNch(); // # of charged tracks in CB
  if (fNch < 1) {
    MarkEndBuffers();
    return;
  }

  // Definition of some useful variables
  // Get beam polarization =10 anti-parallel/ =4 parallel
  Int_t nBeamPol=0;
  nBeamPol = (Int_t)fADC[6]&fgBeamPolConst;

  // Choice of the Trigger
  Int_t nCBSumTrig = 0 ; 
  Int_t nTAPSOrTrig = 0;
  Int_t nPulserTrig = 0;
  Int_t nNoTrig     = 0;

  // Trigger selection (Using trigger pattern --- Sasha)
  Int_t iTrig = fParent->GetBitPattern()->GetHits(0,0);
  switch(iTrig) {
  case 0:  nCBSumTrig  = 1; break;
  case 3:  nTAPSOrTrig = 1; break;
  case 4:  nPulserTrig = 1; break;
  default: nNoTrig     = 1;
  }
  if (!nCBSumTrig) return;

  // Get ptr to the CB tracks
  const TA2CentralTrack *tracks = fCB->GetTracks();
 
  fECB = -1000.0;
  fEPID = -1000.0;
  
  // Neutral
  const Int_t *iNe = fCB->GetIneTracks();
  const TA2CentralTrack *track;
  for (Int_t i=0; i<fNne; ++i) {
    track = tracks + iNe[i];
    fThetaNe[i] = track->GetTheta()*RadToDeg;            // Track theta (polar)   angle (deg)
    fPhiNe[i]   = track->GetPhi()*RadToDeg;              // Track phi (azimuthal) angle (deg)
    fEne[i]     = track->GetEclNaI();                    // Cluster Deposited Energy inside CB (MeV)
    fTne[i]     = track->GetTclNaI();                    // Cluster Time (ns) 
    fMclNe[i]   = track->GetMclNaI();                    // Cluster Size
    fCentralIndexNaINe[i] = track->GetCentralIndexNaI(); // Index of the CB central cluster crystal
  }

  // Charged
  const Int_t *iCh = fCB->GetIchTracks();
  for (Int_t i=0; i<fNch; ++i) {
    track = tracks + iCh[i];
    fThetaCh[i] = track->GetTheta()*RadToDeg;                  // Track theta (polar)   angle (deg)
    fPhiCh[i]   = track->GetPhi()*RadToDeg;                    // Track phi (azimuthal) angle (deg)
    fEch[i]     = track->GetEclNaI();                          // Cluster Deposited Energy inside CB (MeV)
    fTch[i]     = track->GetTclNaI();                          // Cluster Time (ns) 
    fMclCh[i]   = track->GetMclNaI();                          // Cluster Size
    fCentralIndexNaICh[i] = track->GetCentralIndexNaI();       // Index of the CB central cluster crystal
    fEpid[i] = track->GetEhitPid();                            // Energy (MeV) of the hit PID sector 
    fIhitPid[i] = *(fPid->GetHits() + track->GetIhitPid());    // Index of the hit PID sector --new Paolo
    fPsVertexX[i] = track->GetPsVertex().X();                  // Track Pseudo-vertex - X coordinate
    fPsVertexY[i] = track->GetPsVertex().Y();                  // Track Pseudo-vertex - Y coordinate
    fPsVertexZ[i] = track->GetPsVertex().Z();                  // Track Pseudo-vertex - Z coordinate
    fChTrackType[i] = track->GetType();                        // Track Type
    fEmwpc[i] = track->GetEtrackMwpc();                        // Strip Amplitude
    Int_t iClusterCB =  track->GetIclNaI();                    // new Theta from CB
    fThetaChCB[i] = fCB->GetPositionsNaI(iClusterCB)->Theta(); // in rad

//     cout << evtDaq << "\t" << i << "\t" << fThetaCh[i] << "\t" << fPhiCh[i] << "\t" << fEch[i] << "\t" << fEpid[i] << endl;
  }

  // Build 2g invariant mass
  TLorentzVector g[2];
  TVector3 p[2];
  Double_t massTmp, delta, deltaMin = 1000.;
  Double_t theta0, theta1, phi0, phi1;
  Double_t energy0, energy1;
  Double_t Pi0Theta, Pi0Phi, Pi0Energy;
  Bool_t GoodPi0 = kFALSE;

  // gamma per pi0 --> energy, theta, phi
  for (Int_t i=0; i<fNtracks-1; ++i) {
    if (!tracks[i].IsNeutral()) continue;
    p[0] = tracks[i].GetDirCos().Unit() * tracks[i].GetEclNaI();
    g[0].SetPxPyPzE(p[0].X(),p[0].Y(),p[0].Z(), tracks[i].GetEclNaI());
    theta0 = tracks[i].GetTheta();
    phi0 = tracks[i].GetPhi();
    energy0 = tracks[i].GetEclNaI();
    for (Int_t j=i+1; j<fNtracks; ++j) {
      if (!tracks[j].IsNeutral()) continue;
      p[1] = tracks[j].GetDirCos().Unit() * tracks[j].GetEclNaI();
      g[1].SetPxPyPzE(p[1].X(),p[1].Y(),p[1].Z(), tracks[j].GetEclNaI());
      theta1 = tracks[j].GetTheta();
      phi1 = tracks[j].GetPhi();
      energy1 = tracks[j].GetEclNaI();
      massTmp = (g[0] + g[1]).Mag()/1000.; // GeV
      if (massTmp < 0.) continue;
      delta = TMath::Abs(massTmp - mpi0);
      if (delta < deltaMin) {
	deltaMin = delta;
	fM2g = massTmp;
	Double_t num = energy0*TMath::Cos(theta0)+energy1*TMath::Cos(theta1);
	Double_t den = TMath::Sqrt(energy0*energy0 + energy1*energy1 + 2*energy0*energy1*TMath::Cos(phi0-phi1));
	if (den==0.) 
	  cout << "In pi0 theta calculation... den=0 " << den << "\t" << fM2g << "\t" <<gAN->GetNEvent() << "\t" << endl;
	// pi0 theta!!!!
	if (den!=0.) 
	  Pi0Theta = TMath::ACos(num/den)*RadToDeg;
	fM2g = massTmp;
	Pi0Energy = energy0 + energy1;
	Double_t num2 = p[0].Mag()*TMath::Sin(theta0)*TMath::Sin(phi0)+p[1].Mag()*TMath::Sin(theta1)*TMath::Sin(phi1);
	Double_t den2 = p[0].Mag()*TMath::Sin(theta0)*TMath::Cos(phi0)+p[1].Mag()*TMath::Sin(theta1)*TMath::Cos(phi1);
	Pi0Phi = TMath::ATan2(num2,den2)*RadToDeg;
	if (0.06 <= fM2g && fM2g <= 0.18)
	  GoodPi0 = kTRUE;
      }
    } 
  }


  Bool_t GoodProton = kFALSE;
  Bool_t GoodPion = kFALSE;
  Bool_t GoodDelta = kFALSE;
  Bool_t GoodMwpcProt = kFALSE;
  Bool_t GoodMwpcPion = kFALSE;
  Double_t ProtonEne, ProtonTheta, ProtonPhi, ProtonMom, DeltaPhi;
  Double_t PiEne, PiTheta, PiPhi, PiMom;
  Double_t NeutronTheta, NeutronPhi;
  Double_t ECB_prot = -999, ECB_pion = -999;
  Double_t CorrECB_prot = -999, CorrECB_pion = -999;
  Double_t EPID_prot = -999, EPID_norm_prot = -999;
  Double_t EPID_pion = -999, EPID_norm_pion = -999;
  Double_t ECB = -999, EPID = -999, EPID_norm = -999;
  Double_t *xp, *yp;
  Double_t fThetaPion = -999, fThetaMwpcPion = -999;
  Double_t fThetaProton = -999, fThetaMwpcProton = -999;
 
  Int_t HitPID = 0;
  for (Int_t ii = 0; ii<fNch; ii++) {
    track = tracks + iCh[ii];
    if (track->HasNaI() && track->HasPid()) {
      HitPID = fIhitPid[ii];
      xp = &x[HitPID][0];
      yp = &y[HitPID][0];
      fECB = fEch[ii];
      fEPID = fEpid[ii];
      fEPIDs =  fEPID*TMath::Sin(fThetaChCB[ii]);
      fPID_dE[HitPID][IndPID[HitPID]] = fEPIDs;
      IndPID[HitPID]++;
      fTotal_Energy[nCB] = fECB;
      nCB++;
  //     cout << evtDaq << "\t" << ii << "\t" << fECB << "\t" << fEPID << "\t" << fTotal_Energy[nCB-1] << endl;
      if (fECB==0)
	cout << "fECB = 0" << endl;
      if (fTotal_Energy[nCB-1]==0)
	cout << "fTotal = 0 " << endl;

      // if this is a proton, fill in the proton band
      if (TMath::IsInside(fECB, fEPIDs, 8, xp, yp)) {
 	fPID_dE_prot[HitPID][IndPID_prot[HitPID]] = fEPIDs;
 	IndPID_prot[HitPID]++;
	fTotal_Energy_prot[nCB_prot] = fECB;
	nCB_prot++;
	
 	if (fNch==1) {
 	  GoodProton = kTRUE;
 	  fThetaProton = fThetaChCB[ii]*RadToDeg;
	  
	  if (fChTrackType[ii]==11 || fChTrackType[ii]==13 || fChTrackType[ii]==15) {
	    GoodMwpcProt = kTRUE;
	    fThetaMwpcProton = fThetaCh[ii];
	    ProtonEne = (fECB + fEPIDs)*1e-03 + mproton;
	    ProtonMom = TMath::Sqrt(ProtonEne*ProtonEne - mproton*mproton);
	    ProtonTheta = fThetaCh[ii];
	    ProtonPhi = fPhiCh[ii];
	    ECB_prot = fECB;
	    EPID_prot = fEPID;
	    EPID_norm_prot = fEPIDs;
	  }
 	}
      }
      
      // if this is a pion, fill in the pion band
      else {
 	fPID_dE_pion[HitPID][IndPID_pion[HitPID]] = fEPIDs;
 	IndPID_pion[HitPID]++;
	fTotal_Energy_pion[nCB_pion] = fECB;
	nCB_pion++;
	if (fNch==1) {
	  GoodPion = kTRUE;
	  fThetaPion = fThetaChCB[ii]*RadToDeg;
	  if (fChTrackType[ii]==11 || fChTrackType[ii]==13 || fChTrackType[ii]==15) {
	    GoodMwpcPion = kTRUE;
	    fThetaMwpcPion = fThetaCh[ii];
	    PiEne = (fECB + fEPIDs)*1e-03 + mpi;
	    PiMom = TMath::Sqrt(PiEne*PiEne - mpi*mpi);
	    PiTheta = fThetaCh[ii];
	    PiPhi = fPhiCh[ii];
	    ECB_pion = fECB;
	    EPID_pion = fEPID;
	    EPID_norm_pion = fEPIDs;
	    
	    // look for 1 neutral cluster (neutron)
	    if (fNne==1) {
	      NeutronTheta = fThetaNe[0];
	      NeutronPhi = fPhiNe[0];
	      // check if pi+ and n are back to back
	      DeltaPhi = 180 - TMath::Abs(PiPhi - NeutronPhi);
	      if (DeltaPhi >=-30 && DeltaPhi <= 30)
		GoodDelta = kTRUE;
	    }
	  }
	}
      }
    }
  }

  //   if (fNch==3 && fChTrackType[0]==15 && fChTrackType[1]==15 && fChTrackType[2]==15)
//   if (fNch==3 && fChTrackType[0]>10 && fChTrackType[1]>10 && fChTrackType[2]>10)
//     fCB->Example(kTRUE);

  Bool_t GoodProtPi0 = kFALSE;
  Bool_t GoodNPiPlus = kFALSE;
  if (GoodProton && GoodMwpcProt && GoodPi0)
    GoodProtPi0 = kTRUE;
  if (GoodPion && GoodMwpcPion && GoodDelta)
    GoodNPiPlus = kTRUE;  

  fFPDNHits=0 ;
  for(UInt_t m=0; m<fFPD->GetNMultihit(); ++m) {
    for(UInt_t i=0; i<fFPD->GetNhitsM(m); ++i) {
      fFPDTimeOR[fFPDNHits] = fFPD->GetTimeORM(m)[i];
      fFPDHits[fFPDNHits]   = fFPD->GetHitsM(m)[i];
      ++fFPDNHits;
    }
  }
  
  Int_t nCBTaggerAll = 0;
  
  for(Int_t iTagg=0; iTagg<fFPDNHits; ++iTagg) {
    // CB-Tagger
    Tvirtual=1000. ;
    Int_t iprompt=0 ; Double_t pippo2 ; Int_t ientr=0 ;
    
    Double_t MinTime = 10000.0; 
    Double_t TShift  = 108.5 ;
    
    for(Int_t ichar=0; ichar<fNch; ++ichar) {
      iprompt=0;
//       cout << "fFPD " << fFPDTimeOR[iTagg] << "\t" << fTch[ichar] << endl;
      Double_t pippo = fTch[ichar] - fFPDTimeOR[iTagg] ;
      EventTime[nEventTime] = pippo;
      nEventTime++;
	   
      if( (pippo > -45) && (pippo < -27)) {
	iprompt=1 ; 
	Tvirtual=100; // fake prompt
      }
	   
      if( (Tvirtual !=100) && (pippo > -600.0) && (pippo < -45.0)) Tvirtual=30;   // fake random
      if( (Tvirtual !=100) && (pippo > -27.0)  && (pippo < 400.0)) Tvirtual=150;  // fake random
      pippo2 = pippo+TShift ;
      ientr=1 ;
      if (TMath::Abs(MinTime) > TMath::Abs(pippo2)) MinTime = pippo2 ;
      ++nCBTaggerAll;
    }
	 
    for(Int_t w = 1; w<OPT_NWIND+1; w++) {  // Compare with lower/upper bounds
	   
      if ((Tvirtual > Window[w][0]) && (Tvirtual < Window[w][1])) { // of all defined windows
	iWindow = w;  
	Double_t Egamma_val = EgammaCal[fFPDHits[iTagg]];
	Egamma[iWindow][nEgamma[iWindow]] = Egamma_val;
	nEgamma[iWindow]++;

	Bool_t GoodEgamma = kFALSE;
	// 	if (Egamma_val<0.400) 
	GoodEgamma = kTRUE;
	
	deltaphi[iWindow][ndeltaphi[iWindow]] = DeltaPhi;
	ndeltaphi[iWindow]++;
	
	if (GoodProtPi0) {
	  BeamPol10[iWindow][nBeamPol10[iWindow]] = fFPDHits[iTagg];
	  nBeamPol10[iWindow]++;
	  // proton
	  PhiProton[iWindow][nPhiProton[iWindow]]     = ProtonPhi;
	  nPhiProton[iWindow]++;
	  ThetaProton[iWindow][nThetaProton[iWindow]] = ProtonTheta;
	  nThetaProton[iWindow]++;
	  EneProton[iWindow][nEneProton[iWindow]] = ProtonEne;
	  nEneProton[iWindow]++;
	  // pi0
	  PhiPi0[iWindow][nPhiPi0[iWindow]] = Pi0Phi;
	  nPhiPi0[iWindow]++;
	  ThetaPi0[iWindow][nThetaPi0[iWindow]] = Pi0Theta;
	  nThetaPi0[iWindow]++;
	  EnePi0[iWindow][nEnePi0[iWindow]] = Pi0Energy;
	  nEnePi0[iWindow]++;

	  Double_t p_prot=0, E_prot=0, T_prot=0, th_pi0=0, p_pi0=0, E_pi0=0, T_pi0=0;
	  Int_t ierr;
	  ComputeKineticEne(Egamma_val, ProtonTheta, p_prot, E_prot, T_prot, th_pi0, p_pi0, E_pi0, T_pi0, ierr);

	  if (ierr!=1) {
	    EneProtonKin[iWindow][nEneProtonKin[iWindow]] = E_prot;
	    nEneProtonKin[iWindow]++;
	    KinProton[iWindow][nKinProton[iWindow]] = T_prot*1e3;
	    nKinProton[iWindow]++;
	    ThetaPi0Kin[iWindow][nThetaPi0Kin[iWindow]] = th_pi0;
	    nThetaPi0Kin[iWindow]++;
	    DiffThetaPi0[iWindow][nDiffThetaPi0[iWindow]] = th_pi0 - PiTheta;
	    nDiffThetaPi0[iWindow]++;
	    EnePi0Kin[iWindow][nEnePi0Kin[iWindow]] = E_pi0*1000;
	    nEnePi0Kin[iWindow]++;
	    KinPi0[iWindow][nKinPi0[iWindow]] = T_pi0*1e3;
	    nKinPi0[iWindow]++;

	    // *   ipart            particle type   (1 = proton / 2 = pion)
	    // *   itarget          = 1 : Hydrogen
	    // *                    = 2 : Deuterium
	    // *                    = 3 : He-3
	    // *                    = 4 : He-4
	    // *                    = 5 : Carbon
	    // *                    = 6 : CH2
	    // *                    = 7 : Butanol
	    
	    // ********** PROTONS ************
	    Int_t ipart = 1;
	    Int_t itarg = 1; 
	    Double_t kinetic_prot = T_prot*1000; 
	    Double_t output[3];
	    Double_t xcos[6];
	    xcos[0] = fPsVertexX[0];
	    xcos[1] = fPsVertexY[0];
	    xcos[2] = fPsVertexZ[0];
	    xcos[3] = TMath::Sin(ProtonTheta*DegToRad)*TMath::Cos(ProtonPhi*DegToRad);
	    xcos[4] = TMath::Sin(ProtonTheta*DegToRad)*TMath::Sin(ProtonPhi*DegToRad);
	    xcos[5] = TMath::Cos(ProtonTheta*DegToRad);

	    fRangeFit->range_los(&ipart, &kinetic_prot, xcos, &itarg, output);

 	    Int_t istop = 2; // CB
 	    Double_t CorrECB_prot = Read_CBProton(ECB_prot);
 	    Double_t CorrEPID_prot = Read_PIDProton(EPID_prot);

	    Double_t e[3];
	    e[0] = CorrEPID_prot;
	    e[1] = CorrECB_prot;
	    e[2] = 0;
	    Double_t ey[3];
	    ey[0] = -0.09036 + 0.1379*CorrEPID_prot;
	    ey[1] = -2.572 + 0.5054*CorrECB_prot - 0.003526*pow(CorrECB_prot,2) + 9.274e-06*pow(CorrECB_prot,3);
	    ey[2] = 0;

	    // test rangefit for proton
	    fRangeFit->range_fit(&ipart, &istop, xcos, e, ey, &itarg);

	    Double_t TRFprot = fRangeFit->GetEner();
	    RanErrProt = fRangeFit->GetIerr();
	    CumulAllProt = fRangeFit->GetCumulus();
	    if (fRangeFit->GetIerr() == 0) {
	      EneRFProt = TRFprot;
 	      CumulProt = fRangeFit->GetCumulus();
	    }
	    chi2Prot = fRangeFit->GetChi2();

	    EnePIDNormProt[iWindow][nEnePIDNormProt[iWindow]] = EPID_norm_prot; // MeV
	    nEnePIDNormProt[iWindow]++;
 	    if (output[1]>0.) {
	      // CB
 	      EneCBProt[iWindow][nEneCBProt[iWindow]] = ECB_prot; // MeV
 	      nEneCBProt[iWindow]++;
	      EneCBKinProt[iWindow][nEneCBKinProt[iWindow]] = output[1]; // MeV
 	      nEneCBKinProt[iWindow]++; 
 	      MissingEneProt[iWindow][nMissingEneProt[iWindow]] = output[1] - ECB_prot; // MeV
 	      nMissingEneProt[iWindow]++;
	      CheckEneCBProt[iWindow][nCheckEneCBProt[iWindow]] = CorrECB_prot;
 	      nCheckEneCBProt[iWindow]++;
	      CheckMEProt[iWindow][nCheckMEProt[iWindow]] = output[1] - CorrECB_prot;
	      nCheckMEProt[iWindow]++;
	      // PID
// 	      fPID_dE_raw[HitPID][iWindow][nfPID_dE_raw[HitPID][iWindow]] = (Int_t) fADC[100+HitPID];
	      if (HitPID==1) {
// 		EnePIDProt[iWindow][nEnePIDProt[iWindow]] = EPID_prot;
 		EnePIDProt[iWindow][nEnePIDProt[iWindow]] = (Int_t) fADC[100+HitPID];
		nEnePIDProt[iWindow]++;
	      }
	      EnePIDKinProt[iWindow][nEnePIDKinProt[iWindow]] = output[0]; // MeV
 	      nEnePIDKinProt[iWindow]++;
 	      MissingEnePIDProt[iWindow][nMissingEnePIDProt[iWindow]] = output[0] - EPID_prot;
 	      nMissingEnePIDProt[iWindow]++;
	      CheckEnePIDProt[iWindow][nCheckEnePIDProt[iWindow]] = CorrEPID_prot;
 	      nCheckEnePIDProt[iWindow]++;

	      if (chi2Prot>5) {
		TrfBadChi2[iWindow][nTrfBadChi2[iWindow]] = TRFprot;
		nTrfBadChi2[iWindow]++;
		TkineBadChi2[iWindow][nTkineBadChi2[iWindow]] = T_prot*1e3;
		nTkineBadChi2[iWindow]++;
	      }
	      if (chi2Prot<=5) {
		TrfOkChi2[iWindow][nTrfOkChi2[iWindow]] = TRFprot;
		nTrfOkChi2[iWindow]++;
		TkineOkChi2[iWindow][nTkineOkChi2[iWindow]] = T_prot*1e3;
		nTkineOkChi2[iWindow]++;
		OkEnePIDNormProt[iWindow][nOkEnePIDNormProt[iWindow]] = CorrEPID_prot*TMath::Sin(fThetaProton*DegToRad);
		nOkEnePIDNormProt[iWindow]++;
		OkEneCBProt[iWindow][nOkEneCBProt[iWindow]] = CorrECB_prot;
		nOkEneCBProt[iWindow]++;
	      }
 	    }
	    // Range Fit for protons as pions
	    ipart = 2;

	    CorrECB_prot = Read_CBPion(ECB_prot);
 	    CorrEPID_prot = Read_PIDPion(EPID_prot);

	    Double_t f[3];
	    f[0] = CorrEPID_prot;
	    f[1] = CorrECB_prot;
	    f[2] = 0;
	    Double_t fy[3];
// 	    fy[0] = 1.972 - 0.2122*CorrEPID_prot; 
// 	    fy[0] /= 2.;
	    fy[0] = -0.09036 + 0.1379*CorrEPID_prot; // stessa risoluz PID per prot e pion
	    fy[1] = 3.175 + 0.1348*CorrECB_prot; 
	    fy[2] = 0;

	    // test rangefit for proton
	    fRangeFit->range_fit(&ipart, &istop, xcos, f, fy, &itarg);
	    chi2ProtPion = fRangeFit->GetChi2();
	  }
	}

	if (GoodNPiPlus) {

	  // pi+
	  PhiPi[iWindow][nPhiPi[iWindow]] = PiPhi;
	  nPhiPi[iWindow]++;
	  ThetaPi[iWindow][nThetaPi[iWindow]] = PiTheta;
	  nThetaPi[iWindow]++;
	  EnePi[iWindow][nEnePi[iWindow]] = PiEnergy;
	  nEnePi[iWindow]++;

	  Double_t p_pi=0, E_pi=0, T_pi=0, th_n=0, p_n=0, E_n=0, T_n=0;
	  Int_t ierr2;
	  ComputeKineticEneNPi(Egamma_val, PiTheta, p_pi, E_pi, T_pi, th_n, p_n, E_n, T_n, ierr2);

 	  if (ierr2!=1) {
 	    EnePiKin[iWindow][nEnePiKin[iWindow]] = E_pi;
 	    nEnePiKin[iWindow]++;
 	    KinPi[iWindow][nKinPi[iWindow]] = T_pi*1e3;
 	    nKinPi[iWindow]++;
 	    ThetaNeuKin[iWindow][nThetaNeuKin[iWindow]] = th_n;
 	    nThetaNeuKin[iWindow]++;
 	    DiffThetaNeu[iWindow][nDiffThetaNeu[iWindow]] = th_n - NeutronTheta;
 	    nDiffThetaNeu[iWindow]++;
 	    EneNeuKin[iWindow][nEneNeuKin[iWindow]] = E_n*1000;
 	    nEnePi0Kin[iWindow]++;
 	    KinNeu[iWindow][nKinNeu[iWindow]] = T_n*1e3;
 	    nKinNeu[iWindow]++;
	    
	    // ********** PIONS ************
	    Int_t ipart = 2;
	    Int_t itarg = 1; 
	    Double_t kinetic_pion = T_pi*1000;
	    Double_t outputN[3];
	    Double_t ycos[6];
	    ycos[0] = fPsVertexX[0];
	    ycos[1] = fPsVertexY[0];
	    ycos[2] = fPsVertexZ[0];
	    ycos[3] = TMath::Sin(PiTheta*DegToRad)*TMath::Cos(PiPhi*DegToRad);
	    ycos[4] = TMath::Sin(PiTheta*DegToRad)*TMath::Sin(PiPhi*DegToRad);
	    ycos[5] = TMath::Cos(PiTheta*DegToRad);
	    
	    fRangeFit->range_los(&ipart, &kinetic_pion, ycos, &itarg, outputN);

 	    Double_t CorrECB_pion = Read_CBPion(ECB_pion);
 	    Double_t CorrEPID_pion = Read_PIDPion(EPID_pion);

	    Double_t ee[3];
	    ee[0] = CorrEPID_pion;
	    ee[1] = CorrECB_pion;
	    ee[2] = 0;
	    Double_t eey[3];
	    // 	    eey[0] = 1.972 - 0.2122*CorrEPID_pion; // vecchia, gia' brutta
	    eey[0] =  -0.09036 + 0.1379*CorrEPID_pion;// stessa risol del PID protoni
// 	    eey[0] /= 2.;
	    eey[1] = 3.175 + 0.1348*CorrECB_pion; 
	    eey[2] = 0;
	    Int_t istop = 2; // CB
	    // test rangefit for pion
	    fRangeFit->range_fit(&ipart, &istop, ycos, ee, eey, &itarg);



 	    EneCBPion[iWindow][nEneCBPion[iWindow]] = ECB_pion; // MeV
 	    nEneCBPion[iWindow]++;

 	    if (outputN[1]>0.) {
	      // CB
 	      CheckEneCBPion[iWindow][nCheckEneCBPion[iWindow]] = CorrECB_pion;
 	      nCheckEneCBPion[iWindow]++;
	      EneCBKinPion[iWindow][nEneCBKinPion[iWindow]] = outputN[1]; // MeV
	      nEneCBKinPion[iWindow]++; 
	      // PID
	      EnePIDKinPion[iWindow][nEnePIDKinPion[iWindow]] = outputN[0]; // MeV
	      nEnePIDKinPion[iWindow]++;
	      EnePIDPion[iWindow][nEnePIDPion[iWindow]] = EPID_pion; // MeV
	      nEnePIDPion[iWindow]++;
	      EnePIDNormPion[iWindow][nEnePIDNormPion[iWindow]] = EPID_norm_pion; // MeV
	      nEnePIDNormPion[iWindow]++;
 	      if (outputN[2]==0.) {
 		EneCBKinStopPion[iWindow][nEneCBKinStopPion[iWindow]] = outputN[1];
 		nEneCBKinStopPion[iWindow]++;
 		if (Egamma_val<0.4) {

		  Double_t TRFpion = fRangeFit->GetEner();
		  RanErrPion = fRangeFit->GetIerr();
		  CumulAllPion = fRangeFit->GetCumulus();
		  if (fRangeFit->GetIerr() == 0) {
		    EneRFPion = TRFpion;
		    CumulPion = fRangeFit->GetCumulus();
		  }
		  chi2Pion = fRangeFit->GetChi2();

		  // CB
 		  EneCBStopPion[iWindow][nEneCBStopPion[iWindow]] = ECB_pion;
 		  nEneCBStopPion[iWindow]++;
		  CheckEneCBStopPion[iWindow][nCheckEneCBStopPion[iWindow]] = CorrECB_pion;
		  nCheckEneCBStopPion[iWindow]++;
 		  EneCBKinLowPion[iWindow][nEneCBKinLowPion[iWindow]] = outputN[1];
 		  nEneCBKinLowPion[iWindow]++;
 		  MissingEnePion[iWindow][nMissingEnePion[iWindow]] = outputN[1] - ECB_pion; // MeV
 		  nMissingEnePion[iWindow]++;
		  CheckMEPion[iWindow][nCheckMEPion[iWindow]] = outputN[1] - CorrECB_pion;
		  nCheckMEPion[iWindow]++;
		  // PID
 		  EnePIDStopPion[iWindow][nEnePIDStopPion[iWindow]] = EPID_pion;
 		  nEnePIDStopPion[iWindow]++;
		  EnePIDKinStopPion[iWindow][nEnePIDKinStopPion[iWindow]] = outputN[0];
 		  nEnePIDKinStopPion[iWindow]++;
		  MissingEnePIDPion[iWindow][nMissingEnePIDPion[iWindow]] = outputN[0] - EPID_pion;
 		  nMissingEnePIDPion[iWindow]++;
		  CheckEnePIDPion[iWindow][nCheckEnePIDPion[iWindow]] = CorrEPID_pion;
		  nCheckEnePIDPion[iWindow]++;
		  if (chi2Pion>5) {
		    TrfBadChi2Pion[iWindow][nTrfBadChi2Pion[iWindow]] = TRFpion;
		    nTrfBadChi2Pion[iWindow]++;
		    TkineBadChi2Pion[iWindow][nTkineBadChi2Pion[iWindow]] = T_pi*1e3;
		    nTkineBadChi2Pion[iWindow]++;
		  }
		  if (chi2Pion<=5) {
		    TrfOkChi2Pion[iWindow][nTrfOkChi2Pion[iWindow]] = TRFpion;
		    nTrfOkChi2Pion[iWindow]++;
		    TkineOkChi2Pion[iWindow][nTkineOkChi2Pion[iWindow]] = T_pi*1e3;
		    nTkineOkChi2Pion[iWindow]++;
		    OkEnePIDNormPion[iWindow][nOkEnePIDNormPion[iWindow]] = CorrEPID_pion*TMath::Sin(fThetaPion*DegToRad);
		    nOkEnePIDNormPion[iWindow]++;
		    OkEneCBPion[iWindow][nOkEneCBPion[iWindow]] = CorrECB_pion;
		    nOkEneCBPion[iWindow]++;
		  }
		}
	      }
	    }
	    // Range Fit for pions as protons
	    ipart = 1;

	    CorrECB_pion = Read_CBProton(ECB_pion);
 	    CorrEPID_pion = Read_PIDProton(EPID_pion);

	    Double_t ff[3];
	    ff[0] = CorrEPID_pion;
	    ff[1] = CorrECB_pion;
	    ff[2] = 0;
	    Double_t ffy[3];
	    ffy[0] = -0.09036 + 0.1379*CorrEPID_pion;
	    ffy[1] = -2.572 + 0.5054*CorrECB_pion - 0.003526*pow(CorrECB_pion,2) + 9.274e-06*pow(CorrECB_pion,3);
	    ffy[2] = 0;

	    // test rangefit for proton
	    fRangeFit->range_fit(&ipart, &istop, ycos, ff, ffy, &itarg);
	    if (outputN[1]>0. && outputN[2]==0. && Egamma_val < 0.4)
	      chi2PionProt = fRangeFit->GetChi2();
 	  }
 	}
      } 
    } 
  }  

  MarkEndBuffers();
}


// ================================================================
void TA2ExampleRangeFit::ResetEvent()
{
  // Set default value for the variables being chenged during an event

  fNtracks = 0;
  fNch = 0;
  fNne = 0;
  fM2g = ENullFloat;
  fECB = ENullFloat;
  fEPID = ENullFloat;
  chi2Prot     = ENullFloat;
  chi2ProtPion = ENullFloat;
  chi2Pion     = ENullFloat;
  chi2PionProt = ENullFloat;
  fFPDTimeOR[0] = EBufferEnd;
  fFPDNHits     = 0 ;
  fFPDHits[0]   = EBufferEnd;
  
  nCB      = 0;
  nCB_prot = 0;
  nCB_pion = 0;
  nCBnoPID = 0;
  nEventTime = 0;
  for(Int_t w=1; w<OPT_NWIND+1; w++) {
      nBeamPol10[w]    = 0;
      nEgamma[w]       = 0;
      nPhiPi[w]        = 0;
      nThetaPi[w]      = 0;
      nPhiPi0[w]       = 0;
      nThetaPi0[w]     = 0;
      nThetaPi0Kin[w]  = 0;
      nThetaNeuKin[w]  = 0;
      nDiffThetaPi0[w] = 0;
      nDiffThetaNeu[w] = 0;
      nEnePi[w]        = 0;
      nEnePiKin[w]     = 0;
      nEnePi0[w]       = 0;
      nEnePi0Kin[w]    = 0;
      nEneNeuKin[w]    = 0;
      nPhiProton[w]    = 0;
      nThetaProton[w]  = 0;
      nEneProton[w]    = 0;
      nEneProtonKin[w] = 0;
      nKinProton[w]    = 0;
      nKinPi0[w]       = 0;
      nKinNeu[w]       = 0;
      nKinPi[w]        = 0;
      ndeltaphi[w]     = 0;
      nTrfBadChi2[w]   = 0;
      nTkineBadChi2[w] = 0;
      nTrfOkChi2[w]    = 0;
      nTkineOkChi2[w]  = 0;
      nTrfBadChi2Pion[w]   = 0;
      nTkineBadChi2Pion[w] = 0;
      nTrfOkChi2Pion[w]    = 0;
      nTkineOkChi2Pion[w]  = 0;
      nEneCBProt[w]      = 0;
      nCheckEneCBProt[w] = 0;
      nCheckEneCBPion[w] = 0;
      nEneCBKinProt[w]   = 0;
      nEnePIDProt[w]     = 0;
      nEnePIDKinProt[w]  = 0;
      nEnePIDNormProt[w] = 0;
      nMissingEneProt[w] = 0;
      nCheckMEProt[w]    = 0;
      nCheckMEPion[w]    = 0;
      nEneCBPion[w]      = 0;
      nEneCBKinPion[w]   = 0;
      nEnePIDPion[w]     = 0;
      nEnePIDKinPion[w]  = 0;
      nEnePIDNormPion[w] = 0;
      nMissingEnePion[w] = 0;
      nEneCBKinLowPion[w]= 0;
      nEneCBStopPion[w]  = 0;
      nEnePIDStopPion[w] = 0;
      nCheckEnePIDProt[w] = 0;
      nCheckEnePIDPion[w] = 0;
      nMissingEnePIDProt[w] = 0;
      nMissingEnePIDPion[w] = 0;
      nEnePIDKinStopPion[w] = 0;
      nEneCBKinStopPion[w]  = 0;
      nCheckEneCBStopPion[w]= 0;
      nOkEnePIDNormProt[w] = 0;
      nOkEneCBProt[w] = 0;
      nOkEnePIDNormPion[w] = 0;
      nOkEneCBPion[w] = 0;
    }
  
  for(Int_t e=0;e<24;e++) {
    IndPID[e]=0;
    IndPID_prot[e] = 0;
    IndPID_pion[e] = 0;
  }
}


// ================================================================
void TA2ExampleRangeFit::MarkEndBuffers()
{
  // Mark end buffers for arrays
  fFPDTimeOR[fFPDNHits] = EBufferEnd;
  fFPDHits[fFPDNHits]   = EBufferEnd;


  // Neutral
  fThetaNe[fNne] = EBufferEnd;
  fPhiNe[fNne]   = EBufferEnd;
  fEne[fNne]     = EBufferEnd;
  fMclNe[fNne]   = EBufferEnd;
  fTne[fNne]     = EBufferEnd;
  fCentralIndexNaINe[fNne] = EBufferEnd;
  // Charged
  fThetaCh[fNch]   = EBufferEnd;
  fThetaChCB[fNch] = EBufferEnd;
  fPhiCh[fNch]     = EBufferEnd;
  fEch[fNch]       = EBufferEnd;
  fMclCh[fNch]     = EBufferEnd;
  fTch[fNch]       = EBufferEnd;
  fCentralIndexNaICh[fNch] = EBufferEnd;
  fEpid[fNch]      = EBufferEnd;
  fIhitPid[fNch]   = EBufferEnd;
  fPsVertexX[fNch] = EBufferEnd;
  fPsVertexY[fNch] = EBufferEnd;
  fPsVertexZ[fNch] = EBufferEnd;
  fPsVertexGoodZ[fNch] = EBufferEnd;
  fChTrackType[fNch] = EBufferEnd;
  fEmwpc[fNch] = EBufferEnd;
  fTotal_Energy[nCB]       = EBufferEnd;
  fTotal_Energy_prot[nCB_prot] = EBufferEnd;
  fTotal_Energy_pion[nCB_pion] = EBufferEnd;
  //
  for(Int_t e=0;e<24;e++)
    {
      fPID_dE[e][IndPID[e]] = EBufferEnd;
      fPID_dE_prot[e][IndPID_prot[e]] = EBufferEnd;
      fPID_dE_pion[e][IndPID_pion[e]] = EBufferEnd;
      fPID_dE_raw[e][IndPID_pion[e]]  = EBufferEnd;
    }
  fDEmwpc[nCBnoPID]= EBufferEnd;
  EventTime[nEventTime] = EBufferEnd;
  for(Int_t w=1; w<OPT_NWIND+1; w++) {
    BeamPol10[w][nBeamPol10[w]]       = EBufferEnd;
    Egamma[w][nEgamma[w]]             = EBufferEnd;
    ThetaProton[w][nThetaProton[w]]   = EBufferEnd;
    PhiProton[w][nPhiProton[w]]       = EBufferEnd;
    EneProton[w][nEneProton[w]]       = EBufferEnd;
    EneProtonKin[w][nEneProtonKin[w]] = EBufferEnd;
    KinProton[w][nKinProton[w]]       = EBufferEnd;
    KinPi0[w][nKinPi0[w]]             = EBufferEnd;
    KinNeu[w][nKinNeu[w]]             = EBufferEnd;
    KinPi[w][nKinPi[w]]               = EBufferEnd;
    PhiPi[w][nPhiPi[w]]               = EBufferEnd;
    ThetaPi[w][nThetaPi[w]]           = EBufferEnd;
    PhiPi0[w][nPhiPi0[w]]             = EBufferEnd;
    ThetaPi0[w][nThetaPi0[w]]         = EBufferEnd;
    ThetaPi0Kin[w][nThetaPi0Kin[w]]   = EBufferEnd;
    ThetaNeuKin[w][nThetaNeuKin[w]]   = EBufferEnd;
    DiffThetaPi0[w][nDiffThetaPi0[w]] = EBufferEnd;
    DiffThetaNeu[w][nDiffThetaNeu[w]] = EBufferEnd;
    EnePi[w][nEnePi[w]]               = EBufferEnd;
    EnePiKin[w][nEnePiKin[w]]         = EBufferEnd;
    EnePi0[w][nEnePi0[w]]             = EBufferEnd;
    EnePi0Kin[w][nEnePi0Kin[w]]       = EBufferEnd;
    EneNeuKin[w][nEneNeuKin[w]]       = EBufferEnd;
    deltaphi[w][ndeltaphi[w]]         = EBufferEnd;
    TrfBadChi2[w][nTrfBadChi2[w]]     = EBufferEnd;
    TkineBadChi2[w][nTkineBadChi2[w]] = EBufferEnd;
    TrfOkChi2[w][nTrfOkChi2[w]]       = EBufferEnd;
    TkineOkChi2[w][nTkineOkChi2[w]]   = EBufferEnd;
    TrfBadChi2Pion[w][nTrfBadChi2Pion[w]]     = EBufferEnd;
    TkineBadChi2Pion[w][nTkineBadChi2Pion[w]] = EBufferEnd;
    TrfOkChi2Pion[w][nTrfOkChi2Pion[w]]       = EBufferEnd;
    TkineOkChi2Pion[w][nTkineOkChi2Pion[w]]   = EBufferEnd;
    EneCBProt[w][nEneCBProt[w]]           = EBufferEnd;
    CheckEneCBProt[w][nCheckEneCBProt[w]] = EBufferEnd;
    CheckEneCBPion[w][nCheckEneCBPion[w]] = EBufferEnd;
    EneCBKinProt[w][nEneCBKinProt[w]]     = EBufferEnd;
    EnePIDProt[w][nEnePIDProt[w]]         = EBufferEnd;
    EnePIDKinProt[w][nEnePIDKinProt[w]]   = EBufferEnd;
    EnePIDNormProt[w][nEnePIDNormProt[w]] = EBufferEnd;
    MissingEneProt[w][nMissingEneProt[w]] = EBufferEnd;
    CheckMEProt[w][nCheckMEProt[w]]       = EBufferEnd;
    CheckMEPion[w][nCheckMEPion[w]]       = EBufferEnd;
    EneCBPion[w][nEneCBPion[w]]           = EBufferEnd;
    EneCBKinPion[w][nEneCBKinPion[w]]     = EBufferEnd;
    EnePIDPion[w][nEnePIDPion[w]]         = EBufferEnd;
    EnePIDKinPion[w][nEnePIDKinPion[w]]   = EBufferEnd;
    EnePIDNormPion[w][nEnePIDNormPion[w]] = EBufferEnd;
    MissingEnePion[w][nMissingEnePion[w]] = EBufferEnd;
    EneCBStopPion[w][nEneCBStopPion[w]]   = EBufferEnd;
    EnePIDStopPion[w][nEnePIDStopPion[w]] = EBufferEnd;
    CheckEnePIDProt[w][nCheckEnePIDProt[w]] = EBufferEnd;
    CheckEnePIDPion[w][nCheckEnePIDPion[w]] = EBufferEnd;
    MissingEnePIDProt[w][nMissingEnePIDProt[w]] = EBufferEnd;
    MissingEnePIDPion[w][nMissingEnePIDPion[w]] = EBufferEnd;
    EnePIDKinStopPion[w][nEnePIDKinStopPion[w]] = EBufferEnd;
    EneCBKinStopPion[w][nEneCBKinStopPion[w]] = EBufferEnd;
    EneCBKinLowPion[w][nEneCBKinLowPion[w]]   = EBufferEnd;
    CheckEneCBStopPion[w][nCheckEneCBStopPion[w]] = EBufferEnd;
    OkEnePIDNormProt[w][nOkEnePIDNormProt[w]] = EBufferEnd;
    OkEneCBProt[w][nOkEneCBProt[w]] = EBufferEnd;
    OkEnePIDNormPion[w][nOkEnePIDNormPion[w]] = EBufferEnd;
    OkEneCBPion[w][nOkEneCBPion[w]] = EBufferEnd;
  }
}

// ================================================================
void TA2ExampleRangeFit::ParseMisc(char* line) {
  FILE* WinFile;
  Char_t sWord[256];

  // Get keyword
  if(sscanf(line, "%s", sWord)!=1) return;

  if(!strcmp("TimeWindows", sWord)) {
    sscanf(line, "%*s %s", WindowFilename);
    printf("Time windows from:\n %s\n", WindowFilename);
    WinFile = fopen(WindowFilename, "r");
    for(Int_t t=1; t<OPT_NWIND+1; t++) {
      fscanf(WinFile, "%lf %lf", &Window[t][0], &Window[t][1]);
      printf("Time window %d from %f to %f\n", t, Window[t][0], Window[t][1]);
    }
    fclose(WinFile);
    return;
  }
}


// ================================================================
void TA2ExampleRangeFit::Speed() 
{
  static Int_t iStart = time(NULL);
  static Int_t nEvents = 0;
  
  nEvents++;
  if(time(NULL) - iStart >= 60) {
    iStart = time(NULL);
    printf("Processed %d Events/sec\n", nEvents/60);
    nEvents = 0;
  }
}


// =============================================================
void TA2ExampleRangeFit::Read_Cal() {
  ifstream infile("/home/susanna/acquroot/a2/acqu/data/Tagger/TagCal883.dat");
  Double_t t2, t3, t1;

  for(int i=1; i<353; i++){
    infile>>t1>>t2>>t3;
    EgammaCal[i] = t2*0.001;
//     cout << EgammaCal[i] << endl;
  }
  infile.close();

  return;
}


// =============================================================
void TA2ExampleRangeFit::ComputeKineticEne(Double_t Egamma, Double_t th_prot, 
					Double_t &p_prot, Double_t &E_prot, Double_t &T_prot,
					Double_t &th_pi, Double_t &p_pi, Double_t &E_pi, Double_t &T_pion, 
					Int_t &ierr) {

  Double_t A = 0, B = 0, C = 0, delta = 0, dm = 0;  
  ierr = 0;

  //   cout << "========== theta prot " << th_prot << endl;
  th_prot = th_prot*DegToRad;
  A = Egamma * TMath::Cos(th_prot);
  B = Egamma * mproton + (mproton*mproton + mproton*mproton - mpi*mpi)/2.;
  C = Egamma + mproton;

  delta = (A*B)*(A*B) + (C*C-A*A)*(B*B-(C*mproton)*(C*mproton));
  if (delta<0) {
    ierr = 1;
    return;
  }
  delta = TMath::Sqrt(delta);

  Double_t p_ch1 = (A*B+delta)/(C*C-A*A);
  Double_t p_ch2 = (A*B-delta)/(C*C-A*A);
  p_prot = p_ch1;
  if (p_ch2 >= 0 )
    ierr = -1;
  if (p_prot < 0.) {
    ierr = 1;   
    return;
  }

  E_prot = TMath::Sqrt(p_prot*p_prot + mproton*mproton);
  E_pi = Egamma + mproton - E_prot;
  p_pi = TMath::Sqrt(E_pi*E_pi - mpi*mpi);
  dm = (mproton*mproton + mproton*mproton - mpi*mpi)/2.;
  th_pi = (dm - mproton*E_prot + Egamma*E_pi)/(Egamma*p_pi);
  T_prot = E_prot - mproton;
  T_pion = E_pi - mpi;
  th_pi = TMath::ACos(th_pi)*RadToDeg;
  //   cout << th_prot*RadToDeg << "\t" <<dm << "\t" << mproton << "\t" << E_prot << "\t" << Egamma << "\t" << E_pi << "\t" << p_pi << "\t" << th_pi << endl;

  return;
}


// =============================================================
void TA2ExampleRangeFit::ComputeKineticEneNPi(Double_t Egamma, Double_t th_pi, 
					  Double_t &p_pi, Double_t &E_pi, Double_t &T_pi,
					  Double_t &th_n, Double_t &p_n, Double_t &E_n, Double_t &T_n, 
					  Int_t &ierr2) {

  Double_t A = 0, B = 0, C = 0, delta = 0, dm = 0;  
  ierr2 = 0;

  th_pi = th_pi*DegToRad;
  A = Egamma * TMath::Cos(th_pi);
  B = Egamma * mproton + (mpi*mpi + mproton*mproton - mneutron*mneutron)/2.;
  C = Egamma + mproton;

  delta = (A*B)*(A*B) + (C*C-A*A)*(B*B-(C*mpi)*(C*mpi));
  if (delta<0) {
    ierr2 = 1;
    return;
  }
  delta = TMath::Sqrt(delta);

  Double_t p_ch1 = (A*B+delta)/(C*C-A*A);
  Double_t p_ch2 = (A*B-delta)/(C*C-A*A);
  p_pi = p_ch1;
  if (p_ch2 >= 0 )
    ierr2 = -1;
  if (p_pi < 0.) {
    ierr2 = 1;   
    return;
  }

  E_pi = TMath::Sqrt(p_pi*p_pi + mpi*mpi);
  E_n = Egamma + mproton - E_pi;
  p_n = TMath::Sqrt(E_n*E_n - mneutron*mneutron);
  dm = (mproton*mproton + mpi*mpi - mneutron*mneutron)/2.;
  th_n = (dm - mproton*E_pi + Egamma*E_n)/(Egamma*p_n);
  T_pi = E_pi - mpi;
  T_n = E_n - mneutron;
  th_n = TMath::ACos(th_n)*RadToDeg;
//   cout << th_pi*RadToDeg << "\t" << p_pi << "\t" <<  E_pi << "\t" << Egamma << "\t" << E_n << "\t" << p_n << "\t" << th_n << endl;

  return;
}



// =============================================================
Double_t TA2ExampleRangeFit::Read_CBProton(Double_t enecb_prot) {
//   Double_t par[2] = {28.6692, 0.754808};
//   Double_t par[2] = {2.151, 0.9773}; // 12/02
  Double_t par[2] = {2.874, 0.9623};

  Double_t newene;
  Double_t t = enecb_prot;
  newene = (t-par[0])/par[1];

  return newene;
}


// =============================================================
Double_t TA2ExampleRangeFit::Read_CBPion(Double_t enecb_pion) {
//   Double_t par[2] = {28.6692, 0.754808};
//   Double_t par[2] = {17.52, 0.833}; // 12/02
  Double_t par[2] = {17.25, 0.84};

  Double_t newene;
  Double_t t = enecb_pion;
  newene = (t-par[0])/par[1];
//   if (iWindow==1)
//     cout << "=================== Pion " << t << "\t" << newene << "\t true: " << enecb_kin << endl;
  return newene;
}


// =============================================================
Double_t TA2ExampleRangeFit::Read_PIDPion(Double_t enepid_pion) {

//   Double_t par[2] = {0.12, 0.1403};
//  Double_t par[2] = {0., 0.2857}; // 1/3.5
//   Double_t par[2] = {0.3989, 0.5009}; // 12/02
  Double_t par[2] = {0., 0.5};

  Double_t newene;
  Double_t t = enepid_pion;
  newene = (t-par[0])/par[1];
//   if (iWindow==1)
//     cout << "=================== Pion " << t << "\t" << newene << "\t true: " << enecb_kin << endl;
  return newene;
}



// =============================================================
Double_t TA2ExampleRangeFit::Read_PIDProton(Double_t enepid_prot) {

//   Double_t par[2] = {0.3989, 0.5009}; // 12/02
  Double_t par[2] = {0.2559, 0.5289};

  Double_t newene;
  Double_t t = enepid_prot;
  newene = (t-par[0])/par[1];
//   if (iWindow==1)
//     cout << "=================== Pion " << t << "\t" << newene << "\t true: " << enecb_kin << endl;
  return newene;
}
