//--Author	JRM Annand   22nd Sep 2013 Adapt TVME_VUPROM for TAPS trigger
//--Rev 	...
//--Update 	JRM Annand  
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_VUPROMT
// GSI VUPROM FPGA based VMEbus module
// Has a multitude of functions, depending on the programing of the FPGA
// TAPS trigger condition master processor
// DAQ control master, of up to 14 VME systems
// Operate as stand-alone master or slave of overall experiment master
// 96 channel scaler
//

#ifndef __TVME_VUPROMT_h__
#define __TVME_VUPROMT_h__

#include "TVME_VUPROM.h"
//
enum{ EVUT_MaxCPU=14 };     // maximum number of CPUs VUPROM can handle

// enumerator specifies internal register index
// For register address offsets see VMEreg_t VUPROMTreg
enum {
  EVUT_Firmware,
  EVUT_TAPSMode,
  EVUT_TAPSStatus,
  EVUT_EvID,
  EVUT_TrigEvID,
  EVUT_StatEvID,
  EVUT_ReadDebug,
  EVUT_SetDebug0,
  EVUT_SetDebug1,
  EVUT_SetDebug2,
  EVUT_SetDebug3,
  EVUT_IntDelay0,
  EVUT_CPUctrl = EVUT_IntDelay0 + EVUT_MaxCPU,
  EVUT_TrigCtrl,
  EVUT_CPUreset,
  EVUT_ACK,
  EVUT_IRQ,
  //
  EVUT_CFDMask,
  EVUT_LED1Mask,
  EVUT_LED2Mask,
  EVUT_VetoMask,
  EVUT_PbWOMask,
  EVUT_PbWOVetoMask,
  EVUT_CFDPreSc,
  EVUT_LED1PreSc,
  EVUT_LED2PreSc, 
  EVUT_VetoPreSc,
  EVUT_CFDM2PreSc,
  EVUT_LED1M2PreSc,
  EVUT_LED2M2PreSc,
  EVUT_VetoM2PreSc,
  EVUT_PulserPreSc,
  EVUT_PbWOPreSc,
  EVUT_PbWOM2PreSc,
  EVUT_PbWOVetoPreSc,
  EVUT_PbWOVetoM2PreSc,
  EVUT_ReadoutBitpattern,  
  //
  EVUT_ScalerClr,
  EVUT_ScalerLoad,
  EVUT_Scaler
};

// For mask register offsets
enum {EM_CFD,EM_LED1,EM_LED2,EM_Veto,EM_PbWO4,EM_PbWO4V};
// For front-end mask registers
enum {EMF_BaF2=6, EMF_Veto=4};
// For Prescale register offsets
enum {EP_CFD,EP_LED1,EP_LED2,EP_Veto,
      EP_CFDM2,EP_LED1M2,EP_LED2M2,EP_VetoM2,
      EP_Pulser,
      EP_PbWO4,EP_PbWO4M2,EP_PbWO4V,EP_PbWO4VM2};

class TVME_VUPROMT : public TVME_VUPROM {
 protected:
  UInt_t fMask[EM_PbWO4V+1];   // mask values master trigger vuprom
  UInt_t fMaskF[EMF_BaF2];     // BaF2 mask values front-end trigger vuprom
  UInt_t fMaskFV[EMF_Veto];    // Veto mask values front-end trigger vuprom
  UInt_t fMaskFPbWO4;          // PbWO4 mask values front-end trigger vuprom
  UInt_t fPreSc[EP_PbWO4VM2+1];// prescale values master trigger vuprom
  Int_t fTAPSMode;             // master or slave
  Bool_t fIsIRQEnabled;        // Readout register activated?
  Bool_t fIsFrontEnd;          // Front end trigger VUPROM?
  Bool_t fIsBaF;               // Front-end BaF2 trigger module?
  Bool_t fIsVeto;              // Front-end Veto trigger module?
  Bool_t fIsPbWO4;             // Front-end PbWO4 trigger module?
 public:
  TVME_VUPROMT( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TVME_VUPROMT();
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual void WaitIRQ();
  virtual void ResetIRQ();
  virtual void EnableIRQ();
  virtual void DisableIRQ(){ fIsIRQEnabled = kFALSE; }
  virtual Bool_t CheckHardID();
  virtual void ReadIRQScaler(void**); 
  virtual void ReadIRQ(void**); 
  // All TDAQmodules have the following procedures, by default dummies
  // Here they provide control of trigger-related hardware
  virtual void RunTrigCtrl();     // Issue start to TCS and enable triggers
  virtual void EndTrigCtrl();     // Disable triggers and Issue stop to TCS
  virtual void ResetTrigCtrl();   // Stop/start the TCS, disable/enable triggers
  virtual void EnableTrigCtrl();  // Enable triggers
  virtual void DisableTrigCtrl(); // Disable triggers
  virtual void SetDebugOut(Int_t, Int_t);    // switch debug output
  virtual void SetIntDelay(Int_t, Int_t);    // interrupt delay
  //virtual void CmdExe(Char_t*);            // execute command from DAQ ctrl
  //virtual Int_t* GetParms(Int_t);          // provide trigger parameters
  virtual void SendEventID(UInt_t);          // send event ID to remote
  virtual UInt_t GetEventID();               // receive event ID from remote
  Bool_t IsFrontEnd(){ return fIsFrontEnd; }
  //
  ClassDef(TVME_VUPROMT,1)   

    };


//-----------------------------------------------------------------------------
inline void TVME_VUPROMT::WaitIRQ( )
{
  // Poll the state of the IRQ-4 flag, bit 1 of the CSR register
  // Set IRQ ACK-2 NIM output if IRQ-4 received (acquisition busy)
  Int_t datum;
  for(;;){
    if( fIsIRQEnabled ){       // "interrupt" enabled?
      datum = Read(EVUT_IRQ);   // poll IRQ register
      if( datum )break;        // exit when IRQ sensed
    }
  }
}

//-----------------------------------------------------------------------------
inline void TVME_VUPROMT::ResetIRQ( )
{
  // Reset IRQ ACK-2 (busy) flag;
  //
  Write(EVUT_ACK, 1);
}

//-----------------------------------------------------------------------------
inline void TVME_VUPROMT::EnableIRQ( )
{
  // Enable triggers via INT-4 input and ACK-2 output
  //
  fIsIRQEnabled = kTRUE;
}

//-----------------------------------------------------------------------------
inline void TVME_VUPROMT::SendEventID(UInt_t ID)
{
  // Send Event ID to remote system
  Write(EVUT_EvID, ID);        // load the register
  Write(EVUT_TrigEvID, 1);     // trigger the sending
  UInt_t id = Read(EVUT_EvID);
  if(id != ID)
    printf("VUPROM event id reg: wrote %d, read back %d\n",ID,id);
}

//-----------------------------------------------------------------------------
inline UInt_t TVME_VUPROMT::GetEventID()
{
  // Read Event ID from remote system
  // Time out after 100 tries
  UInt_t id;
  Int_t i;
  for(i=0;i<100;i++){
    id = Read(EVUT_EvID);           // read the ID register
    if( id != 0xf0000000 ) break;  // 0xffffffff = nothing there yet
    usleep(1);
  }
  Write(EVUT_TrigEvID, 1);          // reset the ID register
  return id;
}

#endif
