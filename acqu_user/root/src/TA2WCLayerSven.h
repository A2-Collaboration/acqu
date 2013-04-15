//--Author      V Lisin      28th Jun 2004  original DAPHNE fortran -> C
//--Update      JRM Annand... 8th Jul 2004  AcquRoot C++ class
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2WCLayer
//
// General methods for single layer wire chamber
//

#ifndef __TA2WCLayerSven_h__
#define __TA2WCLayerSven_h__

#include "TA2WCLayer.h"

#define MARKED 32767

class TA2WCLayerSven : public TA2WCLayer {
 protected:

 public:
  TA2WCLayerSven(const char*, Int_t, Int_t, Int_t, void*); // Normal use
  virtual Int_t GetNElement() {return fNElement;}
  virtual Int_t GetNHits() {return fNHits;}
  virtual Int_t* GetHits() {return fHits;}
  virtual Int_t GetMaxClust() {return fMaxClust;}
  virtual void DecodeCluster( Int_t&, Int_t** );
  ClassDef(TA2WCLayerSven,1)
};

//---------------------------------------------------------------------------
inline void TA2WCLayerSven::DecodeCluster( Int_t& nHit, Int_t** phit )
{
  fIdist = 2;

  Int_t nh;
  Int_t* hit = *phit;

  fNClust = 0;
  //for(Int_t t=0; t<fNElement+1; t++)
    //fLayerHits[t] = EBufferEnd; //Seems not to be necessary...

  // Get number hits nh in particular layer
  for(nh=0; nh<nHit; nh++)
  {
    if(hit[nh] >= fIend) break;;
    fLayerHits[nh] = hit[nh] - fIstart;
  }
  fLayerHits[nh] = EBufferEnd;  // mark end of local hits buffer
  fNHits = nh;                  // store local # hits

  if(nh==0) return;
  nHit -= nh;                // subtract hits in current layer
  *phit += nh;               // point past hits in current layer

  if(nh==1) // trivial case...1 hit only in layer
  {
    fNClust = 1;
    fLenClust[0] = 1;
    fClust[0] = fLayerHits[0];
    return;
  }

  //----------------------------------------------------------------------------
  //This is Sven's rewritten cluster algorithm for strips/wires
  Int_t Temp[nh+2];
  Int_t Cnt;
  Int_t Clr;
  for(Int_t t=0; t<nh; t++) //Copy hits to a 2nd field for
    Temp[t] = fLayerHits[t];  //marking connected elements
  Temp[nh] = Temp[0];   //This "closes the ring"
  Temp[nh+1] = Temp[1]; //This "closes the ring"

  for(Int_t t=0; t<nh-1; t++)                 //In 2nd field, all elements connected to a
    if(fLayerHits[t+1] - fLayerHits[t]<=fIdist) //previous element are marked. The starting
      Temp[t+1]=MARKED;                       //element of the cluster stays unmarked.

  Cnt = 0;
  for(Int_t t=0; t<nh; t++) //Count all unmarked (=starting elements of a
    if(Temp[t]!=MARKED)     //cluster) elements in the hit pattern.
    {
      fClust[Cnt] = Temp[t];
      fLenClust[Cnt] = 1;
      Clr = 1;
      while(Temp[t+Clr]==MARKED) //The number of marked elements
      {                          //after an unmarked one is the
        fLenClust[Cnt]++;        //length of the detected cluster
        Clr++;
      }                     //Count all unmarked (=starting elements of a
      Cnt++;                //cluster) elements in the hit pattern.
      if(Cnt>=fMaxClust) //Exit if memory is full
      {
        fNClust = Cnt;
        return;
      }
    }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //This is Sven's Debug Code. Don't use it!
  /*for(Int_t t=0; t<nh; t++)
    printf("%d\t",fLayerHits[t]);
  printf("\n");
  for(Int_t t=0; t<nh; t++)
    printf("%d\t",Temp[t]);
  printf("%d", Cnt);
  for(Int_t t=0; t<Cnt; t++)
    printf("\n -> %d: %d",t, fLenClust[t]);
  printf("\n-------------------------------\n");*/
  //----------------------------------------------------------------------------

  fNClust = Cnt;
}

#endif
