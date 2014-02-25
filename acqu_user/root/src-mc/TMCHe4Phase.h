
//--Author	JRM Annand   19th Oct 2005  Skeleton user version	
//--Rev      
//--Update      
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TMCHe4Phase
//
// Skeleton User Particle specification and kinematics generation for 
// the Monte Carlo event generator TMCGenerator.
// Users can redefine the behaviour of TMCParticle here

#ifndef __TMCHe4Phase_h__
#define __TMCHe4Phase_h__

#include "TA2System.h"
#include "TRandom.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"

static const Double_t Mp = 0.938272029;
static const Double_t Mn = 0.9395535;
static const Double_t Md = 1.87561282;
static const Double_t Mamu = 0.931494043;
static const Double_t Me = 0.000510998;


class TMCHe4Phase : public TA2System {
 protected:
  TRandom* fRand;
  TGenPhaseSpace fEvent2;                        // 2-body final state
  TGenPhaseSpace fEvent3;                        // 3-body final state
  TLorentzVector fTarget;
  TLorentzVector fBeam;
  TLorentzVector fW;
  TLorentzVector* fPn;
  Double_t fQ2;
  Double_t fQ3;
  Double_t fTh_n;
  Double_t fT_n;
  Double_t fT_nr;
  Double_t fEbeam;
  Double_t fEb;
  Double_t fdEb;
  Double_t fM2[2];
  Double_t fM3[3];
  Double_t* fThDet;
  Double_t fdThDet;
  Double_t* fArcMi;
  Double_t* fArcPl;
  Double_t* fXpl;
  Double_t* fXmi;
  Double_t* fYpl;
  Double_t* fYmi;
  Double_t* fGrad;
  Double_t fRfp;
  Double_t fRb;
  Double_t fRt;
  Double_t fRdet;
  Double_t fDfp;
  Double_t fTfp;
  Double_t fdTfp;
  Int_t fNdet;
  Int_t fIdet;
  Int_t fNevent;
  Bool_t fIsHist;
  Bool_t fIsHit;
  //
  // Histograms
  //
  TH1F* hEbeam;
  TH2F* hDalitz;
  //   TH2F h2_12("T-Theta_2","h12", 200,0.0,50., 180,0.0,180.);
  //   TH2F h2_13("Time-Theta_2","h13", 200,20.0,120., 180,0.0,180.);
  //   TH2F h2_14("Tr-Theta_2","h12", 200,0.0,50., 180,0.0,180.);
  //   TH2F h3_12("T-Theta_3","h12", 200,0.0,50., 180,0.0,180.);
  //   TH2F h3_13("Time-Theta_3","h13", 200,20.0,120., 180,0.0,180.);
  //   TH2F h3_14("Tr-Theta_3","h12", 200,0.0,50., 180,0.0,180.);
  //   TH2F h23_12("T-Theta_23","h12", 200,0.0,50., 180,0.0,180.);
  //   TH2F h23_13("Time-Theta_23","h13", 200,20.0,120., 180,0.0,180.);
  //   TH2F h23_14("Tr-Theta_23","h12", 200,0.0,50., 180,0.0,180.);
  TH1F** hTm2;  // time
  TH1F** hTm3;
  TH1F** hDfp;
  TH1F** hTm23;
  TH1F** hT2B;   // kinetic energy
  TH1F** hT2;   // kinetic energy
  TH1F** hT3;
  TH1F** hT23;
 public:
  TMCHe4Phase( const Char_t*, Int_t=999999 );
  virtual ~TMCHe4Phase();
  virtual void SetConfig(Char_t*, Int_t);
  virtual void PostInit();
  virtual void FlightPath(TLorentzVector*);
  virtual void Generate();
  virtual void GenTotal();
  virtual void Gen2Body();
  virtual void Gen3Body();
  // Getters
  //  virtual Double_t Mass( ){ return fMass; }
  //
  ClassDef(TMCHe4Phase,1)

};

//-----------------------------------------------------------------------------
inline void TMCHe4Phase::GenTotal()
{
  // Generate total 4 momentum
  fEbeam = fRand->Uniform( fEb-fdEb, fEb+fdEb );
  fBeam.SetPxPyPzE( 0.0, 0.0, fEbeam, fEbeam );
  fW = fBeam + fTarget;
  if( fIsHist ){
    hEbeam->Fill( fEbeam*1000.0 );
  }
}

//-----------------------------------------------------------------------------
inline void TMCHe4Phase::Gen2Body()
{
  // Generate total 2-body phase space
  // Save neutron 4 momentum, polar angle and kinetic energy

  Double_t rwt,weight;               // to iron out weighting
  fEvent2.SetDecay(fW, 2, fM2);
  do{
    weight = fEvent2.Generate();
    rwt = fRand->Uniform(0.0, fEvent2.GetWtMax());
  }while( rwt > weight );
  fPn    = fEvent2.GetDecay(1);
  fTh_n = fPn->Theta() * TMath::RadToDeg();
  fT_n = 1000.0 * (fPn->E() - Mn);
  if( fIsHist ){
    //    h2_12.Fill( T_n, Theta_n );
    //    h23_12.Fill( T_n, Theta_n );
  }
}
//-----------------------------------------------------------------------------
inline void TMCHe4Phase::Gen3Body()
{
  // Generate 3-body phase space
  // Save neutron 4 momentum, polar angle and kinetic energy

  TLorentzVector P12;                // for Dalitz plot check
  TLorentzVector P13;                // for Dalitz plot check
  Double_t rwt,weight;               // to iron out weighting
  fEvent3.SetDecay(fW, 3, fM3);
  do{
    weight = fEvent3.Generate();
    rwt = fRand->Uniform(0.0, fEvent3.GetWtMax());
  }while( rwt > weight );
  //
  fPn    = fEvent3.GetDecay(2);
  fTh_n = fPn->Theta() * TMath::RadToDeg();
  fT_n = 1000.0 * (fPn->E() - Mn);
  P12 = *(fEvent3.GetDecay(0)) + *(fEvent3.GetDecay(1));
  P13 = *(fEvent3.GetDecay(0)) + *(fEvent3.GetDecay(2));
  if( fIsHist ){
     hDalitz->Fill(P12.M2(), P13.M2());
     //     h3_12.Fill( T_n, Theta_n );
     //     h23_12.Fill( T_n, Theta_n );
  }
}

#endif

