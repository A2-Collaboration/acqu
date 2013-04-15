//--Author	JRM Annand         4th Feb 2003
//--Rev 	JRM Annand...     26th Feb 2003  1st production tagger
//--Rev 	K   Livingston... 21th May 2004  Major KL mods
//--Rev 	JRM Annand....... 10th Feb 2005  gcc 3.4 compatible
//--Update	JRM Annand....... 21st Apr 2005  IsPrompt(i) etc.,fMaxParticle
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2KensTagger
//
// Mainz photon-tagging spectrometer apparatus class
// Tagger can consist of 1 or 2 ladders (FP, Microscope) plus PbGlass
// (not yet implemened).
// FP and Micro are instances of the TA2KensLadder class
//

// A typical setup file is like this:
//#       NMR field       Beam Energy
//Misc:   1.01455         853.87
//#
//#         Name          Class Name      Setup file
//Detector: Ladder        TA2KensLadder   FPladderScal.cal
//Detector: Micro         TA2KensLadder   Microscope_kl.dat
//#Condition: xxx
//Reconstruct:
//Initialise:
//Display ...etc
//
// If 2 Ladders are used, then the main FP (ladder) must be before the Micro
// in the setup file

#include "TA2KensTagger.h"

enum { ETaggLadder, ETaggBeamMon, ETaggMicroscope, ETaggGenDet, ETaggPhotonEnergy };

ClassImp(TA2KensTagger)

//-----------------------------------------------------------------------------

TA2KensTagger::TA2KensTagger(const char* name, TA2System* analysis)
              :TA2Tagger(name, analysis)
{
  IsSim = false;

  CircularRuns = 0;
  CircularDeg[0] = 1.0;
  UseCircular = false;
  CurrentRun[0] = '\0';
}

//-----------------------------------------------------------------------------

TA2KensTagger::~TA2KensTagger()
{

}

//-----------------------------------------------------------------------------

void TA2KensTagger::LoadVariable()
{
  TA2Tagger::LoadVariable();

  TA2DataManager::LoadVariable("DegreeMAMI", &fBeamPol, EDSingleX);
  TA2DataManager::LoadVariable("SupWt",      &SupWt,    EDSingleX);
}

//-----------------------------------------------------------------------------

TA2Detector* TA2KensTagger::CreateChild(const char* name, Int_t dclass)
{
  // Create a TA2Detector class for use with the TA2KensTagger
  // Valid detector types are...
  // 1. TA2KensLadder

  if(!name) name = Map2String(dclass);
  switch (dclass)
  {
   case ETaggLadder:
    return new TA2KensLadder(name, this);
    break;
   default:
    return TA2Tagger::CreateChild(name, dclass);
    break;
  }
}

//-----------------------------------------------------------------------------

void TA2KensTagger::Reconstruct( )
{
  // From "hits" provided by the TA2Ladder detector classes
  // reconstruct valid Gamma-Ray energies for the event

  if(UseCircular)
  {
    //Get name of file being analysed
    gUAN->ReadRunName(RunName);
    //Check if file name has changed since last event (i.e. are we analysing a new file)
    if(strcmp(RunName, CurrentRun))
    {
      strcpy(CurrentRun, RunName);
      //Find current run in energy scale table
      for(iRun=0; iRun<CircularRuns; iRun++)
        if(!strcmp(RunName, CircularRun[iRun])) break;
    }
    //Get average MAMI polarisation
    fBeamPol = CircularDeg[iRun];
    SupWt = 1e-6;
  }

  Double_t *FPEnergy, *FPTime, *MicroEnergy, *MicroTime;
  Int_t nFPNhits, nMicroNhits;
  Int_t nPhoton = 0;
  Int_t i;
  UInt_t j;
  Bool_t ignoreMicro = kTRUE;

  Int_t nFPNmultihit;
  Int_t nNFPMhits;
  const Double_t* FPEnergyM;
  //Double_t** LaddTimeM;
  //UInt_t* LaddNhitsM;
  Int_t* LaddHitsM;

  nFPNmultihit = fFP->GetNMultihit();
  UInt_t LaddNhitsM[nFPNmultihit];
  Double_t* LaddTimeM[nFPNmultihit];

  FPEnergy = fFP->GetEelecOR();
  FPTime = fFP->GetTimeOR();
  nFPNhits = fFP->GetNhits();
  if(fIsMicro)
  {
    MicroEnergy = fMicro->GetEelecOR();
    MicroTime = fMicro->GetTimeOR();
    nMicroNhits = fMicro->GetNhits();
  }
  //if no microscope loop round the focal plane
  if((!fIsMicro) || ignoreMicro)
  {
    if(!nFPNmultihit)
    {
      for(i=0; i<nFPNhits; i++)
      {
        fPhotonEnergy[i] = fBeamEnergy - FPEnergy[i];
        fPhotonTime[i] = FPTime[i];
        fP4[i].SetPxPyPzE(0.0, 0.0, fPhotonEnergy[i], fPhotonEnergy[i]);
        SetParticleInfo(i, fFP->GetHits()[i]);
      }
      fPhotonEnergy[nFPNhits] = (Double_t)EBufferEnd;
      fPhotonTime[nFPNhits] = (Double_t)EBufferEnd;
      fNparticle = nFPNhits;
    }
    else
    {
      nNFPMhits = 0;
      //LaddNhitsM = new UInt_t[nFPNmultihit];
      //LaddTimeM = new Double_t*[nFPNmultihit];
      for(i=0; i<nFPNmultihit; i++)
      {
        LaddNhitsM[i] = fFP->GetNhitsM(i);
        LaddTimeM[i] = fFP->GetTimeORM(i);
        LaddHitsM = fFP->GetHitsM(i);
        FPEnergyM = fFP->GetECalibration();
        for(j=0; j<LaddNhitsM[i]; j++)
        {
          fPhotonEnergy[nNFPMhits] = fBeamEnergy - FPEnergyM[LaddHitsM[j]];
          fPhotonTime[nNFPMhits] = LaddTimeM[i][j];

          fP4[nNFPMhits].SetPxPyPzE(0.0, 0.0, fPhotonEnergy[nNFPMhits], fPhotonEnergy[nNFPMhits]);
          SetParticleInfo(nNFPMhits, LaddHitsM[j]);
          nNFPMhits++;
        }
        fPhotonEnergy[nNFPMhits] = (Double_t)EBufferEnd;
        fPhotonTime[nNFPMhits] = (Double_t)EBufferEnd;
        fNparticle = nNFPMhits;
      }
      //delete[] LaddNhitsM;
      //delete[] LaddTimeM;
    }
  }
  // if Micro, use micro values in the energy range it covers
  // ignore FP hits just above and below microscope range 
  // to prevent double counting.
  else
  {
    for(i=0; i<nFPNhits; i++)
    {
      if(FPEnergy[i] > 0.9*fMicroEcalib[0]) //if into Micro range
      {
        for(int j=0; j<nMicroNhits; j++) //do micro hits
        {
          fPhotonEnergy[nPhoton] = fBeamEnergy - MicroEnergy[j];
          fPhotonTime[nPhoton] = MicroTime[j];
          fP4[nPhoton].SetPxPyPzE(0.0, 0.0, fPhotonEnergy[nPhoton], fPhotonEnergy[nPhoton]);
          nPhoton++;
        }
        break;
      }
      else
      {
        fPhotonEnergy[nPhoton] = fBeamEnergy - FPEnergy[i];
        fPhotonTime[nPhoton] = FPTime[i];
        fP4[nPhoton].SetPxPyPzE(0.0, 0.0, fPhotonEnergy[nPhoton], fPhotonEnergy[nPhoton]);
        nPhoton++;
      }
    }
    for(int k=i;k<nFPNhits;k++)
    {
      if(FPEnergy[k] > 1.1*fMicroEcalib[fNMicroElem-1])
      {
        fPhotonEnergy[nPhoton] = fBeamEnergy - FPEnergy[k];
        fPhotonTime[nPhoton] = FPTime[k];
        fP4[nPhoton].SetPxPyPzE(0.0, 0.0, fPhotonEnergy[nPhoton], fPhotonEnergy[nPhoton]);
        SetParticleInfo(nPhoton, fFP->GetHits()[nPhoton]);
        fParticles[nPhoton].SetDetector(EDetMicro);
        nPhoton++;
      }
    }
    fNparticle = nPhoton;
  }
  fPhotonEnergy[fNparticle] = (Double_t)EBufferEnd;
  fPhotonTime[fNparticle] = (Double_t)EBufferEnd;

  // Do any linear polarisation stuff if there's a scaler read
  // and linear polarisation is defined
  if( gAR->IsScalerRead() && fIsLinPol ) ReconstructLinPol(); 
}

//-----------------------------------------------------------------------------

void TA2KensTagger::ParseMisc(char* line)
{
  FILE* CircFile;
  Char_t sWord[256];

  //Get keyword
  if(sscanf(line, "%s", sWord)!=1) return;

  if(!strcmp("CircularFile", sWord))
  {
    sscanf(line, "%*s %s", CircularFile);
    printf("Circular beam polarisation degrees from:\n %s\n", CircularFile);
    CircFile = fopen(CircularFile, "r");
    CircularRuns = 0;
    while(!feof(CircFile))
      if(fscanf(CircFile, "%s %lf", CircularRun[CircularRuns], &CircularDeg[CircularRuns])==2)
      {
        CircularRuns++;
        if(CircularRuns > MAXRUNS) break;
      }
    fclose(CircFile);
    UseCircular = true;
    return;
  }

  if(!strcmp("IsSimulation", sWord))
  {
    sscanf(line, "%*s %d", &IsSim);
    printf("Analysed tagger data is simulated: %d\n", IsSim);
    return;
  }
}

//-----------------------------------------------------------------------------
