#ifndef __TA2BeamPolMon_h__
#define __TA2BeamPolMon_h__

#include "TA2Apparatus.h"
#include "TA2PairSpec.h"               
#include "TA2Moeller.h"               

class TA2BeamPolMon : public TA2Apparatus {
protected:
  TA2PairSpec* fPairSpec;                      // Pair spectrometer
  TA2Moeller*  fMoeller;                      // Moeller electron polarisation
  

public:
  TA2BeamPolMon( const char*, TA2System* ); // pass ptr to analyser
  virtual ~TA2BeamPolMon();
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void Reconstruct() {}
  // Root needs this line for incorporation in dictionary
  ClassDef(TA2BeamPolMon,1)
};



#endif
