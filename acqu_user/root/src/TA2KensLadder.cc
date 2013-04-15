//--Author      JRM Annand        12th Feb 2003
//--Rev         JRM Annand...     26th Feb 2003...1st production version
//--Rev         JRM Annand...      5th Feb 2004...3v8 compatibility
//--Rev         K   Livingston... 20th Apr 2004...changed fLadder to fHits
//--Rev         K   Livingston... 20th Apr 2004...PostInit() init fLaddDoubles
//--Rev         K   Livingston... 16th Jun 2004...Prompt and rand timing wind
//--Rev         K   Livingston... 25th Apr 2005...Multi random windows 10 max
//--Rev         JRM Annand...     16th May 2005...ReadDecoded for Monte Carlo
//--Update      JRM Annand...     25th Jul 2005...SetConfig hierarchy
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2KensLadder
//
// "Linear" array of simple detector elements which provide hit
// timing and rates info. The hit-position is related to an energy scale
// e.g. as in focal-plane plastic scintillator array
// of the Glasgow-Mainz photon tagger
// Expect this to be a high-rate device...so multiple hits.
//
// In general, a particle may go through a single element, or 2 overlapping
// elements. Elements are defined as for TA2Detector, but have extra setup
// parameters relating to energy calibration and scalers.
// An element can have any of the following:
// QDC, TDC, Position, Energy, OverlapEnergy, Scaler
// For details, see the source code and the setupfile.

#include "TA2KensLadder.h"
#include "TA2KensTagger.h"         // need this one for Monte Carlo decode

// constants for command-line maps below
enum { 
  ELaddMisc=100,
  ELaddTrig=101,
  ELaddInit=102,
  ELaddDisplay=103,
  ELaddWindow=104
};

ClassImp(TA2KensLadder)

//---------------------------------------------------------------------------

TA2KensLadder::TA2KensLadder(const char* name, TA2System* apparatus)
              :TA2Ladder(name, apparatus)
{

}

//---------------------------------------------------------------------------

TA2KensLadder::~TA2KensLadder()
{

}

//-----------------------------------------------------------------------------

void TA2KensLadder::LoadVariable()
{
  TA2Ladder::LoadVariable();

  TA2DataManager::LoadVariable("Livetime",      &fLivetime,     EDSingleX);
  TA2DataManager::LoadVariable("Livetagg",      &fLivetagg,     EDSingleX);
  TA2DataManager::LoadVariable("Livecorr",      &fLivecorr,     EDSingleX);
  TA2DataManager::LoadVariable("Interrupts",    &fInterrupts,   EDSingleX);
  TA2DataManager::LoadVariable("LadderP2",      &fLadderP2,     EDSingleX);
  TA2DataManager::LoadVariable("Faraday",       &fFaraday,      EDSingleX);
}

//---------------------------------------------------------------------------

void TA2KensLadder::SetConfig(char* line, int key)
{
  switch(key)
  {
   case ELaddMisc:
    ParseMisc(line); 
    break;
  }

  TA2Ladder::SetConfig(line, key);
}

//---------------------------------------------------------------------------

void TA2KensLadder::ParseMisc(char* line)
{
  Char_t keyword[256];

  //Get keyword
  if(sscanf(line, "%s", keyword)!=1) return;

  if(!strcmp("Dead", keyword))
  {
    sscanf(line, "%*s %d", &iDead[nDead]);
    printf("Assuming ladder channel %d as dead in simulation\n", iDead[nDead]);
    nDead++;
  }
}

//---------------------------------------------------------------------------

void TA2KensLadder::PostInit()
{
  TA2Ladder::PostInit();

  //Overall counter for scaler read events, used for...
  iEvents = 0;
  //...history spectra
  LadderP2Hist   = new TH1F("History_LadderP2",   "History_LadderP2",   1000000, 0, 1000000);
  InterruptsHist = new TH1F("History_Interrupts", "History_Interrupts", 1000000, 0, 1000000);
  LivetimeHist   = new TH1F("History_Livetime",   "History_Livetime",   1000000, 0, 1000000);
  ScalerSumHist  = new TH1F("History_ScalerSum",  "History_ScalerSum",  1000000, 0, 1000000);
  FaradayHist    = new TH1F("History_Faraday",    "History_Faraday",    1000000, 0, 1000000);

  //Total flux monitor from P2
  Char_t Buff[256];
  sprintf(Buff, "%s_P2Flux", fName.Data());
  P2Flux = new TH1F(Buff, Buff, 3, 0, 3);

  //RNG for simulation
  pRandoms = new TRandom3;
  pRandoms->SetSeed(0);

  nDead = 0;
  iDead = new UInt_t[fNelem];
}

//---------------------------------------------------------------------------

void TA2KensLadder::Decode()
{
  TA2Ladder::Decode();

  //Reset, so that not filled in normal events
  fLadderP2 = EBufferEnd;
  fLivetime = EBufferEnd;
  fLivetagg = EBufferEnd;
  fLivecorr = EBufferEnd;
  fInterrupts = EBufferEnd;
  fFaraday = EBufferEnd;

  // If scaler event, and fIsScaler, fill scaler buffs and hists
  if(fIsScaler && gAR->IsScalerRead())
  {
    Double_t LadderSum = 0.0;
    for(UInt_t i=0; i<fNelem; i++)
      LadderSum+=fScaler[fScalerIndex[i]];
    LadderSum/=fNelem;
    fLivetime = (1e0*fScaler[1])/((1e0*fScaler[0]));
    fLivetagg = (1e0*fScaler[14])/((1e0*fScaler[0]));
    fLivecorr = fLivetagg/fLivetime;
    fLadderP2 = (1e0*fScaler[139])/(LadderSum);
    fFaraday = (1e3*fScaler[138])/((1e0*fScaler[0]));
    fInterrupts = (1e6*fScaler[2])/((1e0*fScaler[0]));

    LadderP2Hist->Fill(iEvents+1, fLadderP2);
    LivetimeHist->Fill(iEvents+1, fLivetime);
    ScalerSumHist->Fill(iEvents+1, LadderSum);
    InterruptsHist->Fill(iEvents+1, fInterrupts);
    FaradayHist->Fill(iEvents+1, fFaraday);
    LadderP2Hist->GetXaxis()->SetRange(iEvents-1000, iEvents+1000);
    LivetimeHist->GetXaxis()->SetRange(iEvents-1000, iEvents+1000);
    ScalerSumHist->GetXaxis()->SetRange(iEvents-1000, iEvents+1000);
    InterruptsHist->GetXaxis()->SetRange(iEvents-1000, iEvents+1000);
    FaradayHist->GetXaxis()->SetRange(iEvents-1000, iEvents+1000);
    iEvents++; if(iEvents > 1000000) iEvents = 0;
    P2Flux->Fill(1.0, 1e-6*fScaler[139]);
  }
}

//---------------------------------------------------------------------------

void TA2KensLadder::ReadDecoded()
{
  // Read photon beam parameters from GEANT-3 simulation output
  // Convert to notional residual electron energy to be compatible with
  // what is expected by  TA2KensTagger apparatus
  // See MCBranchID.h for GEANT-3 output details

  // G-3 outputs photon energy in GeV...convert electron energy, MeV
  // Beam array...direction cosines followed by energy.
  Float_t* energy = (Float_t*)(fEvent[EI_beam]);
  Double_t beamEnergy = ((TA2KensTagger*)fParent)->GetBeamEnergy();

  //Search corresponding (imaginary) ladder channel
  UInt_t Channel = EBufferEnd;
  Double_t SimEnergy = beamEnergy - (energy[3] * 1000.0);
  Double_t DEnergy = 1e38;
  for(UInt_t t=0; t<fNelem; t++)
    if(fabs(SimEnergy - fECalibration[t]) < DEnergy)
    {
      Channel = t;
      DEnergy = fabs(SimEnergy - fECalibration[t]);
    }

  //Look, if channel should be dead
  Bool_t bKill = false;
  for(Int_t t=0; t<nDead; t++)
    if(Channel==iDead[t]) bKill = true;

  if(!bKill)
  {
    //Ensure some arrays properly terminated
    //Single prompt hit
    fNhits = 1;

    fHits[0] = Channel;
    fHitsAll[0] = Channel;
    fHitsPrompt[0] = Channel;
    fEelecOR[0] = fECalibration[Channel];
    fTimeOR[0] = pRandoms->Gaus(0.0, 0.75);

    fHits[1] = EBufferEnd;
    fHitsAll[1] = EBufferEnd;
    fHitsPrompt[1] = EBufferEnd;
    fEelecOR[1] = EBufferEnd;
    fTimeOR[1] = EBufferEnd;
  }
  else
  {
    fNhits = 0;

    fHits[0] = EBufferEnd;
    fHitsAll[0] = EBufferEnd;
    fHitsPrompt[0] = EBufferEnd;
    fEelecOR[0] = EBufferEnd;
    fTimeOR[0] = EBufferEnd;
  }

  //Terminate unused stuff
  fNDoubles = 0;
  if(fIsOverlap) fDoubles[0] = EBufferEnd;
  fHitsRand[0] = EBufferEnd;
  fWindows[0] = ELaddPrompt;
}

//---------------------------------------------------------------------------
