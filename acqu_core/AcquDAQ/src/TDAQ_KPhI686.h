//--Author	JRM Annand     21st Jan 2006  
//--Rev 	JRM Annand...
//--Rev         B.Oussena       1st Jan Start updates for KpHBoard
//--Rev 	JRM Annand..   28th Apr 2009  memory mapping
//--Rev         B. Oussena     10th Jul 2009  Add  cmd EKPhI686Board
//--Rev         B. Oussena     15th Jul 2009  Add Init()
//--Rev         B. Oussena     20th Jul 2009  fix bugs in SetConfig(), MapMem()
//--Rev         B. Oussena     3rd Mar 2010   Set KpHboard  to perform 32/24bits
//--Rev         B. Oussena     17th Sep 2010  VME addr modifiers to ModIndex.h
//--Update      JRM Annand      8th Jul 2013  32-bit addr set range
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TDAQ_KPhI686
// Mainz Pentium M based VMEbus SBC
// "Direct bridge" to VMEbus, addresses set by firmware
// Default setup 32-bit addressing to be bits 28-31 = 0xe
// Other rangles, ie bits 28-31 = 0x1,2,3..... can be chosen by setting
// the fRange variable to 0x10000000,0x20000000....
// All 32-bit address slaves in the crate must follow that address scheme

#ifndef __TDAQ_KPhI686_h__
#define __TDAQ_KPhI686_h__

#include "TDAQmodule.h"

// VMEbus address mapping offsets
enum { EVMEbusA24 = 0x9e000000,
       EVMEbusA32 = 0x80000000,
       EVMEbusA16 = 0x9f000000 };

class TDAQ_KPhI686 : public TDAQmodule {
 protected:
  Int_t fMemFd;
  Int_t fRange;
 public:
  TDAQ_KPhI686( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TDAQ_KPhI686();
  void SetConfig( Char_t*, Int_t );         // configure v2718
  virtual void Init(); 
  virtual void Read(void*, void*, Int_t, Int_t);
  virtual void Write(void*, void*, Int_t, Int_t);
  virtual void Read(void*, void*);
  virtual void Write(void*, void*);
  virtual DAQMemMap_t* MapSlave(void*, Int_t, Int_t);
  Int_t GetMemFd(){ return fMemFd; }
  //
  ClassDef(TDAQ_KPhI686,1)   
    };

//-----------------------------------------------------------------------------
inline void TDAQ_KPhI686::Read( void* addr, void* data, Int_t am, Int_t dw )
{
#ifdef VME_HOST
  // General purpose read
  switch(dw){
  case 2:
  default:
    *((UShort_t*)data) = *((UShort_t*)addr);
    break;
  case 4:
    *((UInt_t*)data) = *((UInt_t*)addr);
    break;
  }    
#endif
}

//-----------------------------------------------------------------------------
inline void TDAQ_KPhI686::Write( void* addr, void* data, Int_t am, Int_t dw )
{
#ifdef VME_HOST
  // General purpose write
  switch(dw){
  case 2:
  default:
    *((UShort_t*)addr) = *((UShort_t*)data); 
    break;
  case 4:
    *((UInt_t*)addr) = *((UInt_t*)data);
    break;
  }    
#endif
}

//-----------------------------------------------------------------------------
inline void TDAQ_KPhI686::Read( void* addr, void* data )
{
  //
}

//-----------------------------------------------------------------------------
inline void TDAQ_KPhI686::Write( void* addr, void* data )
{
  //
}



#endif
