//--Author	JRM Annand    4th Apr 2003
//--Rev		JRM Annand...
//--Update	JRM Annand... 1st Dec 2004  Remove MC analysis class
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data.
//
// TA2UserControl
//
// Provide a ROOT command-line interface for Acqu-Root analyser.
// Create TAcquRoot Acqu-Root interface and TA2Analysis data analyser.
// Attempt to configure and run the system as a separate thread.
// Allow user defined analyses....enter analysis class in Map_t kKnownAnalysis
// and in CreateAnalyser() method

#ifndef __TA2UserControl_h__
#define __TA2UserControl_h__


#include "TROOT.h"
#include "TInterpreter.h"
#include "TAcquRoot.h"
#include "TA2Control.h"
#include "TA2Analysis.h"
#include "TA2UserAnalysis.h"

R__EXTERN class TA2UserAnalysis* gUAN;

class TA2UserControl : public TA2Control
{
 protected:
 public:
  TA2UserControl(const char*, int*, char**, void* = NULL, int = -1, Bool_t = 0);
  virtual void CreateAnalyser();

  ClassDef(TA2UserControl,1)
};

#endif
