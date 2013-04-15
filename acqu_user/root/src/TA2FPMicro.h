//--Author	JRM Annand   17th Mar 2006
//--Rev 	
//--Update	
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2FPMicro
//
// Tests of Gla/Mz tagger with high-resolution FP detectors
//

#ifndef __TA2FPMicro_h__
#define __TA2FPMicro_h__

#include "MCBranchID.h"
#include "TA2Detector.h"

class TA2FPMicro : public TA2Detector {
 protected:
  Int_t fNFPelem;
  Int_t fNFPhits;
  Int_t* fFPhits; //[fNFPelem]
 public:
  TA2FPMicro( const char*, TA2System* );// Normal use
  virtual ~TA2FPMicro();
  virtual void PostInit( );            // initialise using setup info
  virtual void ParseDisplay(char*);    // display setup
  virtual void Decode( );              // hits -> energy procedure
  virtual void Cleanup( );             // end-of-event cleanup
  virtual void SaveDecoded( );         // save local analysis
  virtual void ReadDecoded( ){}        // read back previous analysis
  virtual void LoadVariable();         // register variables for plot
  //
  Int_t* GetFPhits(){ return fFPhits; }
  Int_t GetNFPhits(){ return fNFPhits; }
  Int_t GetNFPelem(){ return fNFPelem; }

  ClassDef(TA2FPMicro,1)
};

//---------------------------------------------------------------------------
inline void TA2FPMicro::Decode( )
{
  // Decode the ADC and TDC information in DecodeBasic()
  // Rearrange hits according to adjacent double hits or not
  //

  TA2Detector::DecodeBasic();          // ADC/TDC -> Energy/Time
  UInt_t i,j,k,l,m;
  fNFPhits = 0;
  // loop over "standard" hits
  for( i=0; i<fNhits; i++ ){
    j = fHits[i];
    k = fHits[i+1];
    if( k == (j+1) ) l = 2*j +1;      // adjacent double, assumed overlap
    else l = 2*j;                     // isolated single hit
    if( i ){
      m = fHits[i-1];
      if( m == (j-1) ) continue;
    }
    fFPhits[fNFPhits] = l;
    fNFPhits++;
  }
  fFPhits[fNFPhits] = EBufferEnd;    // end markers in buffers
}

//---------------------------------------------------------------------------
inline void TA2FPMicro::Cleanup( )
{
  // end-of-event cleanup
  TA2Detector::Cleanup();
  // "zero" hits buffer
  fFPhits[0] = EBufferEnd;
}
#endif
