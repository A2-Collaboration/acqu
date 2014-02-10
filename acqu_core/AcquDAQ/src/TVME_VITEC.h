//--Author	JRM Annand   3rd Sep 2013
//--Rev 	JRM Annand
//--Update	JRM Annand  17th Sep 2013  Mods for new firmware
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_VITEC
// VITEC piggy back card to the KPh i686 VMEbus SBC
// Used for "interrupt" and event-ID processing
//

#ifndef __TVME_VITEC_h__
#define __TVME_VITEC_h__

#include "TVMEmodule.h"

enum{
  EVIT_ECLin,
  EVIT_ECLout,
  EVIT_NIMin,
  EVIT_NIMout,
  EVIT_EvIDlsb,
  EVIT_EvIDmsb,
  EVIT_Status,
  EVIT_FwID,
  EVIT_BitPatternLow,
  EVIT_BitPatternHigh
};
enum{ EVIT_EvIDTimeout = 4000 };

class TVME_VITEC : public TVMEmodule {
 protected:
  Bool_t    fIsIRQEnabled;
  UInt_t fEvID;                 // event ID
  UInt_t fReadoutPatternStatus; // status register (also for readout pattern)
  UInt_t fReadoutPatternOffset;     // ADC Index Offset for readout pattern
  UInt_t fReadoutPattern;           // ReadoutPattern (for TAPS fast readout)
 public:
  TVME_VITEC( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TVME_VITEC();
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual Bool_t CheckHardID();
  virtual void WaitIRQ();
  virtual void SetIRQ(){ Write(EVIT_NIMout,(UInt_t)0x1); }
  virtual void ResetIRQ(){ Write(EVIT_NIMout,(UInt_t)0x0); }
  virtual void EnableIRQ(){ fIsIRQEnabled = kTRUE; }
  virtual void DisableIRQ(){ fIsIRQEnabled = kFALSE; }
  virtual UInt_t GetEventID();
  virtual UInt_t GetReadoutPatternOffset();
  virtual UInt_t GetReadoutPatternStatus();
  virtual UInt_t GetReadoutPattern(); 
  Bool_t IsIRQEnabled(){ return fIsIRQEnabled; }

  ClassDef(TVME_VITEC,1)   

    };

//-----------------------------------------------------------------------------
inline void TVME_VITEC::WaitIRQ( )
{
  // Poll the state of the status register
  // When bit 15 set exit loop, 
  // read bitpattern,
  // and set NIM output 1: ACK
  UInt_t datum;
  for(;;){
    if( fIsIRQEnabled ){                 // "interrupt" software enabled?
      datum = Read(EVIT_Status);         // poll status register
      // this happens sometimes for whatever reason
      // then just try again!
      if( datum == 0xffff)
        continue;
      if( datum & 0x8000 ) break;         // exit when bit 15 set
    }
  }
  fReadoutPatternStatus = datum & 0xffff;
  // check bitpattern received
  if((datum & 0x1f00) == 0x1a00) {
    // yes, read it
   fReadoutPattern  =  Read(EVIT_BitPatternLow)  & 0xffff;        // lsb part
   fReadoutPattern |= (Read(EVIT_BitPatternHigh) & 0xffff) << 16; // msb part
  }
  else {
    // failed? nothing received?
    fReadoutPattern = 0xffffffff;
  }
  // set ACK (NIM output 1)
  SetIRQ();
}

//-----------------------------------------------------------------------------
inline UInt_t TVME_VITEC::GetEventID( )
{
  // Read the event ID
  //
  UInt_t datum,evid;
  for(;;){
    datum = Read(EVIT_Status);
    // this happens sometimes for whatever reason
    // then just try again!
    if( datum == 0xffff)
      continue;
    // check if bit4 is high
    if( (datum & 0x10) == 0x10 ) 
      break;
  }

  // something received, check if it's okay
  if((datum & 0x1f) != 0x1a){
    evid = 0xffffffff;
    return evid;
  }

  // event id ok
  evid = Read(EVIT_EvIDlsb);           // get the event ID
  datum = Read(EVIT_EvIDmsb);
  evid = evid | (datum << 16);
  return evid;
}

//-----------------------------------------------------------------------------
inline UInt_t TVME_VITEC::GetReadoutPatternOffset( )
{
  return fReadoutPatternOffset;
}

//-----------------------------------------------------------------------------
inline UInt_t TVME_VITEC::GetReadoutPatternStatus( )
{
  return fReadoutPatternStatus;
}

//-----------------------------------------------------------------------------
inline UInt_t TVME_VITEC::GetReadoutPattern( )
{
  return fReadoutPattern;
}

#endif
