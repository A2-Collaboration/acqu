//--Author      V Lisin      28th Jun 2004  original DAPHNE fortran -> C
//--Update      JRM Annand... 8th Jul 2004  AcquRoot C++ class
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2CylWire
//
// Helical cathode strips on a cylindrical surface
//

#ifndef __TA2CylWireSven_h__
#define __TA2CylWireSven_h__

#include "TA2WCLayerSven.h"

class TA2CylWireSven : public TA2WCLayerSven
{
 protected:
  Double_t fRadius;                // radius mm
  Double_t fLength;                // length mm
  Double_t fPhiCor[3];             // z correction coeffs
  Double_t fPhiSpace;              // phi spacing of wires
 public:
  TA2CylWireSven(const char*, Int_t, Int_t, Int_t, void*, Double_t*);
  //  virtual void CGCluster( Int_t );
  virtual void CGCluster( Int_t );
  virtual Double_t GetRadius(){ return fRadius; }
  virtual Double_t GetLength(){ return fLength; }
  virtual Double_t* GetPhiCor(){ return fPhiCor; }
  virtual Double_t GetPhiCor(Int_t i){ return fPhiCor[i]; }

  ClassDef(TA2CylWireSven,1)
};

//---------------------------------------------------------------------------
inline void TA2CylWireSven::CGCluster( Int_t ic )
{
  // Center of gravity of cluster (index ic) of adjacent hit wires
  // No pulse height to weight the distribution.
  // Convert wire # to phi

  Double_t phiCG = fPhiSpace * ( fClust[ic] + 0.5*(fLenClust[ic] - 1) );
  Double_t twoPi = TMath::TwoPi();
  phiCG += fPhiCor[0];
  if( phiCG >= twoPi ) phiCG -= twoPi;
  else if( phiCG < 0 ) phiCG += twoPi;
  //fCGClust[ic] = phiCG + fPhiCor[0];
}

#endif
