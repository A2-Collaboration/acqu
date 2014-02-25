//--Author	A Neiser      XXth Nov 2013 Adapt from TVME_VUPROM
//--Rev 	...
//--Update 	JRM Annand    31st Jan 2014 Fix register setup 
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_VUPROM Moeller
// Readout of the Moeller data

#ifndef __TVME_VUPROM_Moeller_h__
#define __TVME_VUPROM_Moeller_h__

#include <vector>
#include "TVMEmodule.h"

enum { EVUM_NMoellerScaler = 96 };
enum {
  EVUM_Firmware,
  EVUM_GlobalEnable,
  EVUM_NimOutMode,
  EVUM_DAQ_status,
  EVUM_DAQ_reset,
  EVUM_DAQ_enable,
  EVUM_DAQ_counter_fix,
  EVUM_DAQ_counter_livetime,
  EVUM_DAQ_counter_livetime_gated,
  EVUM_VN2andVN1,
  EVUM_MamiEnergy,
  EVUM_InputChannelsDebugMode,
  EVUM_InputChannelsDebugLeftStart,
  EVUM_InputChannelsDebugRightStart,
  // scaler
  EVUM_ScalerClr,
  EVUM_ScalerLoad,
  EVUM_Scaler,
  // RAM interface
  EVUM_RAM_WEB = EVUM_Scaler + EVUM_NMoellerScaler,
  EVUM_RAM_AddrB,
  EVUM_RAM_DinB,
  EVUM_RAM_DOutB0
};

class TVME_VUPROM_Moeller : public TVMEmodule {
private:
  // fields to be configured
 UShort_t fNLeftChannels;
 UShort_t fNPairsPerCh; 
 UShort_t fNBins; // in read-out Moeller histogram (per channel pair)
 UInt_t fNReadsPerIRQ; // VME reads per IRQ are limited!
 UInt_t fScalerModulo; // scaler IRQs needed to trigger Moeller readout 
 
 Bool_t kCfgMode;
 UInt_t fOffsetLeft, fOffsetRight; // in config mode, holds the two offsets
 
 // internal fields
 UInt_t fNScalerIRQsSeen;
 Bool_t kReadoutStarted;
 UShort_t fReadoutsDone;
 UShort_t iBin, iHelicity, iLeftChannel, iPair;
 Int_t  iRamAddrPrev; // -1 is needed as invalid flag
 Bool_t kChainIsLast;
 
 
 void ClearIndices();
 void StartMoellerDAQ();
 void StopMoellerDAQ();
 
 
protected:

public:
  TVME_VUPROM_Moeller( Char_t*, Char_t*, FILE*, Char_t* );
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual Bool_t CheckHardID();
  virtual void ReadIRQScaler(void**); 
  virtual void ReadIRQ(void**); 
   
  ClassDef(TVME_VUPROM_Moeller,1)

};

#endif
