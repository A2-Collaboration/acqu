//--Author	JRM Annand   21st Nov 2004
//--Rev 	
//--Update	
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2CosmicCal
//
// For general description of calorimeter-type systems

#ifndef __TA2CosmicCal_h__
#define __TA2CosmicCal_h__

#include "TA2Apparatus.h"              // base apparatus class
#include "TA2LongScint.h"              // array long scintillator bars
#include "TA2PlasticPID.h"             // auxiliary trigger detectors

class TA2CosmicCal : public TA2Apparatus {
protected:
  TA2LongScint* fLScint;                   // Long scintillator bars
  TA2LongScint* fLTrig;                    // Long scntillator trigger
  TA2PlasticPID* fTrig;                    // Small trigger scint

public:
  TA2CosmicCal( const char*, TA2System* ); // pass ptr to analyser
  virtual ~TA2CosmicCal();
  virtual void PostInit();                 // some setup after parms read in
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void LoadVariable();             // display setup
  virtual void Reconstruct();
  virtual void Cleanup();                 // reset at end of event

  TA2LongScint* GetLScint(){ return fLScint; }
  TA2LongScint* GetLTrig(){ return fLTrig; }
  TA2PlasticPID* GetTrig(){ return fTrig; }

  // Root needs this line for incorporation in dictionary
  ClassDef(TA2CosmicCal,1)
};

//-----------------------------------------------------------------------------
inline void TA2CosmicCal::Cleanup( )
{
  // Clear any arrays holding variables
 TA2Apparatus::Cleanup();
}

//-----------------------------------------------------------------------------
inline void TA2CosmicCal::Reconstruct( )
{
  // Reconstruct track of cosmic-ray muon

}

#endif
