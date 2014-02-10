//--Author	JRM Annand   21st Nov 2004
//--Rev 	
//--Update	
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2GenericApp
//
// For general purpose detector use

#ifndef __TA2GenericApp_h__
#define __TA2GenericApp_h__

#include "TA2Apparatus.h"              // base apparatus class
#include "TA2GenericDetector.h"        // common detector stuff
#include "TA2LongScint.h"              // array long scintillator bars
#include "TA2PlasticPID.h"             // auxiliary trigger detectors
#include "TA2FPMicro.h"                // tagger microscope tests
#include "TA2PairSpec.h"               // Pair Spectrometer

class TA2GenericApp : public TA2Apparatus {
protected:
  TA2GenericDetector* fGenDet;             // "standard" detector
  TA2LongScint* fLScint;                   // Long scintillator bars
  TA2LongScint* fLTrig;                    // Long scntillator trigger
  TA2PlasticPID* fTrig;                    // Small trigger scint
  TA2FPMicro* fFPMic;                      // microscope tests
  TA2PairSpec* fPairSpec;                      // Pair spectrometer

public:
  TA2GenericApp( const char*, TA2System* ); // pass ptr to analyser
  virtual ~TA2GenericApp();
  virtual void PostInit();                 // some setup after parms read in
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void LoadVariable();             // display setup
  virtual void Reconstruct();
  virtual void Cleanup();                 // reset at end of event

  TA2GenericDetector* GetGenDet(){ return fGenDet; }
  TA2LongScint* GetLScint(){ return fLScint; }
  TA2LongScint* GetLTrig(){ return fLTrig; }
  TA2PlasticPID* GetTrig(){ return fTrig; }

  // Root needs this line for incorporation in dictionary
  ClassDef(TA2GenericApp,1)
};

//-----------------------------------------------------------------------------
inline void TA2GenericApp::Cleanup( )
{
  // Clear any arrays holding variables
 TA2Apparatus::Cleanup();
}

//-----------------------------------------------------------------------------
inline void TA2GenericApp::Reconstruct( )
{
  // Reconstruct track of cosmic-ray muon

}

#endif
