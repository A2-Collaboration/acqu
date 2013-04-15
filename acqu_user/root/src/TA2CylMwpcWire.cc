// SVN info: $Id: TA2CylMwpcWire.cc 68 2011-10-07 19:16:03Z mushkar $
#include "TA2CylMwpcWire.h"

ClassImp(TA2CylMwpcWire)

//---------------------------------------------------------------------------
TA2CylMwpcWire::TA2CylMwpcWire(const char* name, Int_t nelem, Int_t maxclust, Int_t maxclsize, void* det, const Double_t* layerparm)
               :TA2WCLayerSven(name, nelem, maxclust, maxclsize, det)
{
  // Store dimensions, alignment and correction factors
  fRadius   = layerparm[0];
  fLength   = layerparm[1];
  fPhiSpace = 2.*kPi/fNElement;
  
  //
  fNClust = 0;
}
