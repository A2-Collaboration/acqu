//--Author	V Lisin      28th Jun 2004  original DAPHNE fortran -> C
//--Update	JRM Annand... 8th Jul 2004  AcquRoot C++ class
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2WCLayer
//
// General methods for single layer wire chamber
//

#include "TA2WCLayerSven.h"
#include "TA2WireChamber.h"

ClassImp(TA2WCLayerSven)

//---------------------------------------------------------------------------
  TA2WCLayerSven::TA2WCLayerSven(const char* name, Int_t nelem, Int_t maxclust, Int_t maxclsize, void* det)
                 :TA2WCLayer(name, nelem, maxclust, maxclsize, det)
{

}
