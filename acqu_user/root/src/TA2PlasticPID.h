//--Author	JRM Annand   30th Sep 2003
//--Rev 	JRM Annand...15th Oct 2003 ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004 3v8 compatible
//--Rev 	JRM Annand...19th Feb 2004 User code
//--Rev 	JRM Annand...16th May 2005 ReadDecoded (bug G3 output)
//--Rev 	JRM Annand... 2nd Jun 2005 ReadRecoded bug fix (D.Glazier)
//--Update	D.Glazier....24th Aug 2007 Add time for ReadDecoded
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2PlasticPID
//
// Internal Particle Identification Detector for the Crystal Ball
// Cylindrical array of plastic scintillators
//

#ifndef __TA2PlasticPID_h__
#define __TA2PlasticPID_h__

#include "MCBranchID.h"
#include "TA2Detector.h"

class TA2PlasticPID : public TA2Detector
{
 private:
 protected:
  UInt_t fUseSigmaEnergy; //Use energy resoution smearing for MC?
  UInt_t fUseSigmaTime;   //Use time resolution smearing in MC?
  Double_t fSigmaEnergy;  //Sigma for energy resolution
  Double_t fPowerEnergy;  //Sigma for energy resolution
  Double_t fSigmaTime;    //Sigma for time resolution
  Double_t fThrMean;
  Double_t fThrSigma;
 public:
  TA2PlasticPID(const char*, TA2System*); //Normal use
  virtual ~TA2PlasticPID();
  virtual void LoadVariable();            //display/cut setup
  virtual void SetConfig(char*, int);     //read in specific parameters
  virtual void Decode();                  //hits -> energy procedure
  virtual void SaveDecoded();             //save local analysis
  virtual void ReadDecoded();             //read back previous analysis

  ClassDef(TA2PlasticPID,1)
};

//---------------------------------------------------------------------------


#endif

