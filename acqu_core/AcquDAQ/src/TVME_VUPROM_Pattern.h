//--Author	A Neiser      XXth Nov 2013 Adapt from TVME_VUPROM
//--Rev 	...
//--Update 	JRM Annand    30th Jan 2014 Fix register setup 
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_VUPROM_Pattern
// Hit pattern reading module, configurable via dat file

#ifndef __TVME_VUPROM_Pattern_h__
#define __TVME_VUPROM_Pattern_h__

#include <vector>
#include "TVMEmodule.h"


class TVME_VUPROM_Pattern : public TVMEmodule {
private:
  std::vector<VMEreg_t> fVUPROMregs; // instead of this enum/pointer hell

  std::vector< size_t > fPatternOffsets; // indices in fVUPROMregs 
  UInt_t fPatternOffset;
  
protected:

public:
  TVME_VUPROM_Pattern( Char_t*, Char_t*, FILE*, Char_t* );
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual Bool_t CheckHardID();
  virtual void ReadIRQ(void**); 
 
  ClassDef(TVME_VUPROM_Pattern,1)

};

#endif
