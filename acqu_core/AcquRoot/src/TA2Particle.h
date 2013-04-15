//--Author	S Schumann   29th Oct 2008
//--Rev         E J Downie   16th mar 2009 combine SS & JRMA updates
//--Update      JRM Annand   23rd Apr 2009 call TObj constructor
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2Particle
//
// Interface class for particle information transfer between <TA2Apparatus>
// and <TA2Physics>

#ifndef __TA2Particle_h__
#define __TA2Particle_h__

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TA2System.h"
#include "TA2Analysis.h"
#include "TA2Apparatus.h"
#include "TA2ClusterDetector.h"
#include "HitCluster_t.h"

//-----------------------------------------------------------------------------

enum
{
  ECircNeg = -1,
  ECircNone = 0,
  ECircPos = 1
};

enum
{
  ELinNone = 0,
  ELinPara = 1,
  ELinPerp = 2
};

enum
{
  EWindowNone,
  EWindowPrompt,
  EWindowRandom
};

enum
{
 EDetUndef = 0,
 EDetNaI = 1,
 EDetPID = 2,
 EDetMWPC = 4,
 EDetBaF2 = 8,
 EDetPbWO4 = 16,
 EDetVeto = 32,
 EDetLadder = 64,
 EDetMicro = 128,
 EDetEndpoint = 256,
 EDetTOFBar = 512,
 EDetCerenkov = 1024
};

enum
{
 EAppUndef = 0,
 EAppCB = 1,
 EAppTAPS = 2,
 EAppTagger = 4,
 EAppTOF = 8
};

//-----------------------------------------------------------------------------

class TA2Particle : public TObject
{
 private:
 protected:
  TDatabasePDG* fPDG;
  TLorentzVector fP4;
  Double_t fTime;
  Double_t fVetoTime;
  Double_t fVetoEnergy;
  Double_t fCircDegree;
  Double_t fLinDegree;
  Double_t fPSAShort;
  Double_t fSigmaE;
  Double_t fSigmaPhi;
  Double_t fSigmaTheta;
  Double_t fSigmaTime;
  Int_t fClusterSize;
  Int_t fParticleID;
  Int_t fSecondID;
  Int_t fVetoIndex;
  Int_t fCentralIndex;
  Int_t fCircState;
  Int_t fLinState;
  Int_t fWindow;
  Int_t fApparati;
  Int_t fDetectors;
  Bool_t fUnclear;
  Int_t    fIintersMwpc[2]; //Indexes of intersections in the internal (0) and external (1) MWPC
  Double_t fTrackEnergy;    //Sum of MWPCs intersections pulse amplitudes
  Double_t fTrackTime;      //To-do 
  //
 public:
  TA2Particle();
  TA2Particle(TDatabasePDG* pPDG); // Use fParticleID->GetPDG() of TA2Apparatus class
  virtual ~TA2Particle();
  TLorentzVector GetP4(){ return fP4; }
  TLorentzVector* GetP4A(){ return &fP4; }
  TVector3 GetP3(){ return fP4.Vect(); }
  TVector3 GetVect(){ return fP4.Vect(); }
  TVector3 GetDirection(){ return fP4.Vect().Unit(); }
  Double_t GetEnergy(){ return fP4.E(); }
  Double_t GetTotal(){ return fP4.E(); }
  Double_t GetKinetic(){ return (fP4.E() - fP4.M()); }
  Double_t GetMass(){ return fP4.M(); }
  Double_t GetMomentum(){ return fP4.P(); }
  Double_t GetE(){ return fP4.E(); }
  Double_t GetT(){ return (fP4.E() - fP4.M()); }
  Double_t GetM(){ return fP4.M(); }
  Double_t GetP(){ return fP4.P(); }
  Double_t GetPx(){ return fP4.Px(); }
  Double_t GetPy(){ return fP4.Py(); }
  Double_t GetPz(){ return fP4.Pz(); }
  Double_t GetPhi(){ return fP4.Phi(); }
  Double_t GetTheta(){ return fP4.Theta(); }
  Double_t GetPhiDg(){ return (fP4.Phi() * TMath::RadToDeg()); }
  Double_t GetThetaDg(){ return (fP4.Theta() * TMath::RadToDeg()); }
  Double_t GetTime(){ return fTime; }
  Double_t GetLinDegree(){ if((fApparati & EAppTagger))
      return fLinDegree; else return 0.0; }
  Double_t GetCircDegree(){ if((fApparati & EAppTagger))
      return fCircDegree; else return 0.0; }
  Double_t GetPSAShort(){ if((fDetectors & EDetBaF2))
      return fPSAShort; else return (Double_t)EBufferEnd; }
  Double_t GetPSALong(){ if((fDetectors & EDetBaF2)) return GetKinetic();
    else return (Double_t)EBufferEnd; }
  Double_t GetVetoEnergy(){ if((fDetectors & EDetVeto)||(fDetectors & EDetPID)) return fVetoEnergy;
    else return (Double_t)EBufferEnd; }
  Double_t GetVetoTime(){ if((fDetectors & EDetVeto)||(fDetectors & EDetPID)) return fVetoTime;
    else return (Double_t)EBufferEnd; }
  Double_t GetSigmaE(){ return fSigmaE; }
  Double_t GetSigmaPhi(){ return fSigmaPhi; }
  Double_t GetSigmaTheta(){ return fSigmaTheta; }
  Double_t GetSigmaPhiDg(){ return (fSigmaPhi * TMath::RadToDeg()); }
  Double_t GetSigmaThetaDg(){ return (fSigmaTheta * TMath::RadToDeg()); }
  Double_t GetSigmaTime(){ return fSigmaTime; }
  Double_t GetAngle(TA2Particle pParticle){
    return fP4.Vect().Angle(pParticle.GetVect()); }
  Double_t GetAngle(TA2Particle* pParticle){
    return fP4.Vect().Angle(pParticle->GetVect()); }
  Double_t GetAngle(TLorentzVector pP4){ return fP4.Vect().Angle(pP4.Vect()); }
  Double_t GetAngle(TLorentzVector* pP4){
    return fP4.Vect().Angle(pP4->Vect()); }
  Double_t GetAngle(TVector3 pVect){ return fP4.Vect().Angle(pVect); }
  Double_t GetAngle(TVector3* pVect){ return fP4.Vect().Angle(*pVect); }
  Double_t GetAngleDg(TA2Particle pParticle){
    return (fP4.Vect().Angle(pParticle.GetVect()) * TMath::RadToDeg()); }
  Double_t GetAngleDg(TA2Particle* pParticle){
    return (fP4.Vect().Angle(pParticle->GetVect()) * TMath::RadToDeg()); }
  Double_t GetAngleDg(TLorentzVector pP4){
    return (fP4.Vect().Angle(pP4.Vect()) * TMath::RadToDeg()); }
  Double_t GetAngleDg(TLorentzVector* pP4){
    return (fP4.Vect().Angle(pP4->Vect()) * TMath::RadToDeg()); }
  Double_t GetAngleDg(TVector3 pVect){
    return (fP4.Vect().Angle(pVect) * TMath::RadToDeg()); }
  Double_t GetAngleDg(TVector3* pVect){
    return (fP4.Vect().Angle(*pVect) * TMath::RadToDeg()); }
  Int_t GetVetoIndex(){ if((fDetectors & EDetVeto)||(fDetectors & EDetPID)) return fVetoIndex; else return EBufferEnd; }
  Int_t GetCentralIndex(){ if((fDetectors & EDetBaF2)||(fDetectors & EDetNaI)||(fDetectors & EDetPbWO4)||(fDetectors & EDetLadder)
                            ||(fDetectors & EDetMicro)||(fDetectors & EDetEndpoint)||(fDetectors & EDetTOFBar))return fCentralIndex;
    else return EBufferEnd; }
  Int_t GetLinState(){ if((fApparati & EAppTagger)) return fLinState;
    else return ELinNone; }
  Int_t GetCircState(){ if((fApparati & EAppTagger)) return fCircState;
    else return ECircNone; }
  Int_t GetClusterSize(){ return fClusterSize; }
  Int_t GetParticleID(){ return fParticleID; }
  Int_t GetSecondID(){ return fSecondID; }
  Int_t GetWindow(){ if((fDetectors & EDetLadder)||(fDetectors & EDetMicro)) return fWindow; else return EWindowNone; }
  Int_t GetApparati(){ return fApparati; }
  Int_t GetDetectors(){ return fDetectors; }
  Bool_t HasApparatus(Int_t pApparatus);
  Bool_t HasApparatusA(Int_t pApparatus){ return (fApparati & pApparatus); }
  Bool_t HasDetector(Int_t pDetector);
  Bool_t HasDetectorA(Int_t pDetector){ return (fDetectors & pDetector); }
  Bool_t IsUnclear(){ return fUnclear; }
  const Int_t  &GetTrackIntersect(const Int_t iCh) const { return fIintersMwpc[iCh]; }
  const Int_t* GetTrackIntersects() const { return fIintersMwpc; }
  const Double_t &GetTrackEnergy() const { return fTrackEnergy; }
  const Double_t &GetTrackTime()   const { return fTrackTime; }
  void SetTime(Double_t pTime){ fTime = pTime; }
  void SetVetoIndex(Int_t pVetoIndex){ fVetoIndex = pVetoIndex; }
  void SetCentralIndex(Int_t pCentralIndex){ fCentralIndex = pCentralIndex; }
  void SetClusterSize(Int_t pClusterSize){ fClusterSize = pClusterSize; }
  void SetLinState(Int_t pLinState);
  void SetCircState(Int_t pCircState);
  void SetLinDegree(Double_t pLinDegree){ fLinDegree = pLinDegree; }
  void SetCircDegree(Double_t pCircDegree){ fCircDegree = pCircDegree; }
  void SetPSAShort(Double_t pPSAShort){ fPSAShort = pPSAShort; }
  void SetVetoEnergy(Double_t pVetoEnergy){ fVetoEnergy = pVetoEnergy; }
  void SetVetoTime(Double_t pVetoTime){ fVetoTime = pVetoTime; }
  void SetSigmaE(Double_t pSigmaE){ fSigmaE = pSigmaE; }
  void SetSigmaPhi(Double_t pSigmaPhi){ fSigmaPhi = pSigmaPhi; }
  void SetSigmaTheta(Double_t pSigmaTheta){ fSigmaTheta = pSigmaTheta; }
  void SetSigmaPhiDg(Double_t pSigmaPhi){ fSigmaPhi = (pSigmaPhi * TMath::DegToRad()); }
  void SetSigmaThetaDg(Double_t pSigmaTheta){ fSigmaTheta = (pSigmaTheta * TMath::DegToRad()); }
  void SetSigmaTime(Double_t pSigmaTime){ fSigmaTime = pSigmaTime; }
  void SetWindow(Int_t pWindow);
  void SetApparatus(Int_t pApparatus);
  void SetApparatusA(Int_t pApparatus){ fApparati+=pApparatus; };
  void SetDetector(Int_t pDetector);
  void SetDetectorA(Int_t pDetector){ fDetectors+=pDetector; }
  void SetUnclear(Bool_t pUnclear){ fUnclear = pUnclear; }
  void SetP4(TLorentzVector pP4){ fP4 = pP4; }
  void SetP4(TLorentzVector* pP4){ if(pP4) fP4 = *pP4; }
  void SetP4(Int_t pParticleID, Double_t pT, TVector3 pP_);
  void SetP4(Int_t pParticleID, Double_t pT, TVector3* pP_){
    if(pP_) SetP4(pParticleID, pT, *pP_); }
  void SetParticleID(Int_t pParticleID);
  void SetParticleIDA(Int_t pParticleID){ fParticleID = pParticleID; }
  void SetSecondID(Int_t pSecondID){ fSecondID = pSecondID; }
  void SetMass(Double_t pM);
  void SetKinetic(Double_t pT);
  void SetTotal(Double_t pE);
  void SetEnergy(Double_t pE){ SetTotal(pE); }
  void SetMomentum(Double_t pP);
  void SetDirection(TVector3 pP_){ fP4.SetVect(fP4.P() * pP_.Unit()); }
  void SetDirection(TVector3* pP_){
    if(pP_) fP4.SetVect(fP4.P() * pP_->Unit()); }
  void SetP3(TVector3 pP3);
  void SetP3(TVector3* pP3) { if(pP3) SetP3(*pP3); }
  void SetM(Double_t pM) { SetMass(pM); }
  void SetT(Double_t pT){ SetKinetic(pT); }
  void SetE(Double_t pE){ SetTotal(pE); }
  void SetP(Double_t pP){ SetMomentum(pP); }
  void SetPx(Double_t pPx);
  void SetPy(Double_t pPy);
  void SetPz(Double_t pPz);
  void SetPhi(Double_t pPhi){ fP4.SetPhi(pPhi); }
  void SetTheta(Double_t pTheta){ fP4.SetTheta(pTheta); };
  void SetPhiDg(Double_t pPhi){ fP4.SetPhi(pPhi * TMath::DegToRad()); }
  void SetThetaDg(Double_t pTheta){ fP4.SetTheta(pTheta * TMath::DegToRad()); }
  void SetPDG(TDatabasePDG* pPDG){ fPDG = pPDG; }
  void SetTrackIntersects(const Int_t i0, const Int_t i1) { fIintersMwpc[0] = i0; fIintersMwpc[1] = i1; }
  void SetTrackIntersect(const Int_t iCh, const Int_t i) { fIintersMwpc[iCh] = i; }
  void SetTrackEnergy(const Double_t &e) { fTrackEnergy = e; }
  void SetTrackTime(const Double_t &t) { fTrackTime = t; }
  void Reset();
  TA2Particle Add(TA2Particle pParticle);
  TA2Particle Add(TA2Particle* pParticle){ if(pParticle) return Add(*pParticle); else return *this; }
  TA2Particle Add(TA2Particle pParticle1, TA2Particle pParticle2);
  TA2Particle Add(TA2Particle* pParticle1, TA2Particle* pParticle2){ if(pParticle1||pParticle2) return Add(*pParticle1, *pParticle2); else return *this; }
  TA2Particle Add(TA2Particle pParticle1, TA2Particle* pParticle2){ if(pParticle2) return Add(pParticle1, *pParticle2); else return *this; }
  TA2Particle Add(TA2Particle* pParticle1, TA2Particle pParticle2){ if(pParticle1) return Add(*pParticle1, pParticle2); else return *this; }
  TA2Particle Add(TA2Particle* pParticles, Int_t pNumber);
  HitCluster_t* FindCluster(const Char_t*, const Char_t*);
  TDatabasePDG* GetPDG(){ return fPDG; }

  ClassDef(TA2Particle , 1)
};

//-----------------------------------------------------------------------------

#endif
