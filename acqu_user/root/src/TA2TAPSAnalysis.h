// SVN Info: $Id: TA2MyRateEstimation.h 455 2009-06-21 19:27:15Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2013
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TA2TAPSAnalysis                                                      //
//                                                                      //
// Class used for the analysis of TAPS stand-alone data.                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef __TA2TAPSAnalysis_h__
#define __TA2TAPSAnalysis_h__


#include "TA2Physics.h"
#include "TA2GenericApp.h"
#include "TA2GenericDetector.h"


class TA2TAPSAnalysis : public TA2Physics
{

private:
    TA2GenericApp* fTAPS;                                   // TAPS apparatus
    TA2GenericDetector* fBaF2_LG;                           // BaF2 LG detector
    TA2GenericDetector* fBaF2_LGS;                          // BaF2 LGS detector
    TA2GenericDetector* fBaF2_SG;                           // BaF2 SG detector
    TA2GenericDetector* fBaF2_SGS;                          // BaF2 SGS detector
    TA2GenericDetector* fVeto;                              // Veto detector
    TA2GenericDetector* fPWO;                               // PWO detector
    TA2GenericDetector* fPWO_S;                             // PWO sensitive detector
    TA2GenericDetector* fPWO_Veto;                          // PWO Veto detector
    TA2GenericDetector* fPWO_Veto_S;                        // PWO Veto sensitive detector
    
    TA2BitPattern* fTAPSPattern;                            // pointer to the TAPS bit pattern
    Int_t fBaF2CFDPattern;                                  // index of the BaF2 CFD pattern in the TAPS bit pattern
    Int_t fBaF2LED1Pattern;                                 // index of the BaF2 LED1 pattern in the TAPS bit pattern
    Int_t fBaF2LED2Pattern;                                 // index of the BaF2 LED2 pattern in the TAPS bit pattern
    Int_t fVetoLEDPattern;                                  // index of the Veto LED pattern in the TAPS bit pattern
    UInt_t fBaF2CFDNhits;                                   // number of CFD hits in the BaF2
    Int_t* fBaF2CFDHits;                                    // pointer to the BaF2 CFD hits
    UInt_t fBaF2LED1Nhits;                                  // number of LED1 hits in the BaF2
    Int_t* fBaF2LED1Hits;                                   // pointer to the BaF2 LED1 hits
    UInt_t fBaF2LED2Nhits;                                  // number of LED2 hits in the BaF2
    Int_t* fBaF2LED2Hits;                                   // pointer to the BaF2 LED2 hits
    UInt_t fVetoLEDNhits;                                   // number of LED hits in the Veto
    Int_t* fVetoLEDHits;                                    // pointer to the Veto LED hits
    
    TH1** fH_BaF2_LG;                                       // BaF2 LG ADC 
    TH1** fH_BaF2_LG_CFD;                                   // BaF2 LG ADC with CFD cut
    TH1** fH_BaF2_LG_LED1;                                  // BaF2 LG ADC with LED1 cut
    TH1** fH_BaF2_LG_LED2;                                  // BaF2 LG ADC with LED2 cut
    TH1** fH_BaF2_Time;                                     // BaF2 LG TDC
    
    TH1** fH_BaF2_LGS;                                      // BaF2 LGS ADC 
    TH1** fH_BaF2_LGS_CFD;                                  // BaF2 LGS ADC with CFD cut
    TH1** fH_BaF2_LGS_LED1;                                 // BaF2 LGS ADC with LED1 cut
    TH1** fH_BaF2_LGS_LED2;                                 // BaF2 LGS ADC with LED2 cut
    
    TH1** fH_BaF2_SG;                                       // BaF2 SG ADC 
    TH1** fH_BaF2_SGS;                                      // BaF2 SG ADC 
    
    TH1** fH_Veto;                                          // Veto ADC 
    TH1** fH_Veto_Time;                                     // Veto TDC 
    
    TH1** fH_PWO;                                           // PWO ADC
    TH1** fH_PWO_CFD;                                       // PWO ADC with CFD (TDC) cut
    TH1** fH_PWO_Time;                                      // PWO TDC
    
    TH1** fH_PWO_S;                                         // PWO ADC sensitive
    TH1** fH_PWO_S_CFD;                                     // PWO ADC sensitive with CFD (TDC) cut
    TH1** fH_PWO_S_Time;                                    // PWO TDC
    
    TH1** fH_PWO_Veto;                                      // PWO Veto ADC
    TH1** fH_PWO_Veto_CFD;                                  // PWO Veto ADC with CFD (TDC) cut
    TH1** fH_PWO_Veto_Time;                                 // PWO Veto TDC
    
    TH1** fH_PWO_Veto_S;                                    // PWO ADC sensitive
    TH1** fH_PWO_Veto_S_CFD;                                // PWO ADC sensitive with CFD (TDC) cut
    TH1** fH_PWO_Veto_S_Time;                               // PWO Veto TDC
    
    Int_t GetPatternIndex(const Char_t* patternName, TA2BitPattern* pattern);
    void FillADCHits(TA2Detector* det, TH1** h);
    void FillTDCHits(TA2Detector* det, TH1** h);
    void FillTDCHitsM(TA2Detector* det, TH1** h);
    
public:
    TA2TAPSAnalysis(const char* name, TA2Analysis* analysis);
    virtual ~TA2TAPSAnalysis();

    void SetConfig(Char_t* line, Int_t key);
    void PostInit();
    void Reconstruct();
    
    ClassDef(TA2TAPSAnalysis, 1)
};

#endif

