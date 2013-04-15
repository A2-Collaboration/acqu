#ifndef __TA2TriggerPhysics_h__
#define __TA2TriggerPhysics_h__

#include "TAcquRoot.h"
#include "TA2Analysis.h"
#include "TA2CB.h"
#include "TA2Taps.h"

#define APP_CB   1
#define APP_TAPS 2

class TA2Apparatus;

//-----------------------------------------------------------------------------

class TA2TriggerPhysics : public TA2Physics
{
  protected:
    TA2CB* CB;           //Crystal Ball apparatus
    TA2Taps* TAPS;       //TAPS apparatus
    TA2CalArray* NaI;    //NaI array
    TA2TAPS_BaF2* BaF2;  //BaF2 array
    Double_t* ThresNaI;  //Energy spectra for M-trigger-relevant NaI hits
    Double_t* ThresBaF2; //Energy spectra for M-trigger-relevant BaF2 hits

    void Speed();                    //Do speed measurement (Events/sec)
    void Swap(Double_t*, Double_t*); //Swap array fields (Double_t)
    void Swap(Int_t*, Int_t*);       //Swap array fields (Int_t)
  public:
    TA2TriggerPhysics(const char*, TA2Analysis*);
    virtual ~TA2TriggerPhysics();
    virtual void LoadVariable();
    virtual void PostInit();
    virtual void Reconstruct();
    virtual TA2DataManager* CreateChild(const char*, Int_t){return NULL;}

  ClassDef(TA2TriggerPhysics,1)
};

//-----------------------------------------------------------------------------

#endif
