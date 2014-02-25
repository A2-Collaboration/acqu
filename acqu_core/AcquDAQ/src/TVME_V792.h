//--Author	JRM Annand   1st Feb 2007
//--Rev 	JRM Annand
//--Rev 	JRM Annand  11th Jan 2008  1st working version
//--Rev 	JRM Annand..28th Apr 2009..remove TDAQmemmap.h
//--Rev 	JRM Annand.. 2nd May 2011..no pedestal suppress option
//--Update	JRM Annand..13th Apr 2013..tidy up coding...fix pedestals
//--Update	D. Werthmueller..27th Jan 2014..added support for high threshold range
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_V792
// CAEN VMEbus 32-channel QDC
// Also forms basis for classes to handle similar CAEN ADC hardware
//

#ifndef __TVME_V792_h__
#define __TVME_V792_h__

#include "TVMEmodule.h"

enum { 
  EV7XX_NMem = 512,
  EV7XX_NThresh = 32,
};
enum {
  EV7XX_Outbuff,
  EV7XX_Firmware = EV7XX_Outbuff + EV7XX_NMem,
  EV7XX_BitSet1,
  EV7XX_BitClr1,
  EV7XX_Status1,
  EV7XX_Ctrl1,
  EV7XX_Status2,
  EV7XX_EvCntL,
  EV7XX_EvCntH,
  EV7XX_FCLRWind,
  EV7XX_BitSet2,
  EV7XX_BitClr2,
  EV7XX_EvCntRes,
  EV7XX_IPed,
  EV7XX_Thresholds,
  EV7XX_ID0 = EV7XX_Thresholds + EV7XX_NThresh,
  EV7XX_ID1,
  EV7XX_ID2,
};
//
class TVME_V792 : public TVMEmodule {
 protected:
  Bool_t fIsThreshold;
  Bool_t fIsHighThrRange;
  Int_t fThresh[EV7XX_NThresh];
 public:
  TVME_V792( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TVME_V792();
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual void ReadIRQ( void** );
  virtual Bool_t CheckHardID();
  virtual Bool_t IsThreshold(){ return fIsThreshold; }
  //
  virtual void ResetData(){ Write(EV7XX_BitSet2); Write(EV7XX_BitClr2); }
  virtual void Reset(){ Write(EV7XX_BitSet1); Write(EV7XX_BitClr1); }

  ClassDef(TVME_V792,1)   

    };

#endif

