//--Author      JRM Annand        12th Feb 2003
//--Rev         JRM Annand...     26th Feb 2003...1st production version
//--Rev         JRM Annand...      5th Feb 2004...3v8 compatibility
//--Rev         K   Livingston... 20th Apr 2004...changed fLadder to fHits
//--Rev         K   Livingston... 20th Apr 2004...PostInit() init fLaddDoubles
//--Rev         K   Livingston... 16th Jun 2004...Prompt and rand timing wind
//--Rev         K   Livingston... 25th Apr 2005...Multi random windows 10 max
//--Rev         JRM Annand...     16th May 2005...ReadDecoded for Monter Carlo
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

#ifndef __TA2KensLadder_h__
#define __TA2KensLadder_h__

#include "TRandom3.h"

#include "TA2Ladder.h"        // basic detector
#include "MCBranchID.h"         // GEANT-3 output format defs

//---------------------------------------------------------------------------

class TA2KensLadder : public TA2Ladder
{
 protected:
  TRandom3* pRandoms;
  Double_t fLivetime;
  Double_t fLivetagg;
  Double_t fLivecorr;
  Double_t fInterrupts;
  Double_t fLadderP2;
  Double_t fFaraday;
  Double_t fP2Flux;
  UInt_t* iDead;
  Int_t nDead;
  Int_t iEvents;
  TH1F* LadderP2Hist;
  TH1F* InterruptsHist;
  TH1F* ScalerSumHist;
  TH1F* LivetimeHist;
  TH1F* FaradayHist;
  TH1F* P2Flux;
public:
  TA2KensLadder(const char*, TA2System*);// Normal use
  virtual ~TA2KensLadder();
  virtual void SetConfig(char*, int);   // decode and load setup info
  virtual void LoadVariable();          // display setup
  virtual void PostInit();              // initialise using setup info
  virtual void Decode();                // hits -> energy procedure
  virtual void ReadDecoded();           // read back previous analysis
  virtual void ParseMisc(char*);
  ClassDef(TA2KensLadder,1)
};

//---------------------------------------------------------------------------

#endif
