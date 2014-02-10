#ifndef __TA2PairSpec_h__
#define __TA2PairSpec_h__

#include "TA2Detector.h"            // Acqu-Root histogrammer

class TA2PairSpec : public TA2Detector
{
private:
  UInt_t fNchannels;
  UInt_t fNHistograms;
  UInt_t fScalerOffset;
  UInt_t fScalerBlockSize;
  
  UInt_t* fScalerOpen;
  UInt_t* fScalerGated;    // Gated with PairSpec
  UInt_t* fScalerGatedDly; // Gated, but out of coincidence with PairSpec
  
  Double_t* fScalerSumOpen;
  Double_t* fScalerSumGated;    // Gated with PairSpec
  Double_t* fScalerSumGatedDly; // Gated, but out of coincidence with PairSpec
  
  UInt_t fScalerBlockIdx; // used in SetConfig, counts ScalerBlocks within Vuprom
  UInt_t fVupromSize;    // number of ScalerBlocks in Vuprom
  UInt_t fVupromSumSize; // used in SetConfig, sum of VupromSize
public:
  TA2PairSpec( const char*, TA2System* );  // pass ptr to analyser
  virtual ~TA2PairSpec();                  // destructor
  virtual void PostInit();                    // some setup after parms read in
  virtual void LoadVariable();              // display setup
  virtual void Decode();                // event by event analysis
  virtual void Cleanup() {}                     // reset at end of event
  virtual void SetConfig( Char_t*, Int_t );   // setup decode in implement
  virtual void SaveDecoded() {}
  virtual void ReadDecoded() {}
 
  // Root needs this line for incorporation in dictionary
  ClassDef(TA2PairSpec,1)
};

#endif
