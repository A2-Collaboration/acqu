//--Author	JRM Annand   18th Feb 2004   Use def physics
//--Rev
//--Update
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2UserPhysics
//
// User-defined physics class
//

#ifndef __TA2UserPhysics_h__
#define __TA2UserPhysics_h__

#include "TA2Physics.h"

class TA2Apparatus;

class TA2UserPhysics : public TA2Physics {
 protected:
  TA2Apparatus* fTAGG;
  TA2Apparatus* fCB;
  TA2Apparatus* fTAPS;
  TLorentzVector* fP4beam;
  TLorentzVector fP4input;
  Int_t fNbeam;
  Double_t fDpnAngle;
  Double_t* fDpnAngleCM;
  Double_t fPhiDiff;
  Double_t fThetaDiff;
  Double_t fThetaP;
  Double_t fThetaN;
  Double_t* fThetaPCM;
  Double_t* fThetaNCM;
  Double_t* fThetaDiffCM;
  Double_t fDpnMinv;
  Double_t* fDpnMmiss;
  Double_t fCBTAPS2PhotonMinv;
 public:
  TA2UserPhysics( const char*, TA2Analysis* );
  virtual ~TA2UserPhysics();
  virtual void LoadVariable();            // variables for display/cuts
  virtual void PostInit( );               // init after parameter input
  virtual void Reconstruct();             // reconstruct detector info
  virtual TA2DataManager* CreateChild( const char*, Int_t ){ return NULL;}

  ClassDef(TA2UserPhysics,1)
};

#endif
