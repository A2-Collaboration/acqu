//--Author      R Gregor   3th Jul 2005
//
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TAPS
//
// User coded version of TAPS apparatus for identification
// of charged particles.

#ifndef __TA2TAPS_h__
#define __TA2TAPS_h__

#include "TA2Particle.h"

#include "TA2Apparatus.h"            // Acqu-Root histogrammer
#include "TA2TAPS_BaF2.h"
#include <iostream>

class TA2KensTagger;
class TA2KensLadder;

class TA2TAPS : public TA2Apparatus
{
 protected:
  TA2KensTagger* fTAGG;                     // Glasgow photon tagger
  TA2KensLadder* fLadder;
  TA2TAPS_BaF2* fBaF2;    // Pointer to BaF2 array
  Int_t* RealVetoHits;
  Bool_t* IsProton;
  Int_t* VetoMap;
  UInt_t fNWVetosAmount;
  Double_t fTapsSigma1;
  Double_t fTapsSigma2;
  UInt_t fNWVetos;
  UInt_t fVetoStatus;
  Double_t fTapsFudge;
  Int_t fTapsTof;
  Int_t fTapsPsa;
  Double_t fTapsVetodet;
  Int_t fNCharged;
  Bool_t* isVCharged;
  Int_t* psaInfo;
  Double_t* fTestLG;
  Double_t* fTestSG;
  Double_t* psarad;
  Double_t* psaang;
  Int_t* minRad;
  Int_t* maxRad;
  Double_t* angle1;
  Double_t* angle2;
  Double_t* angle3;
  Double_t* angle4;
  Double_t* angle5;
  Double_t* angle6;
  Double_t* angle7;
  Double_t* angle8;
  Double_t* angle9;
  Double_t* angle10;
  Double_t* angle11;
  Double_t* sigma1;
  Double_t* sigma2;
  Double_t* sigma3;
  Double_t* sigma4;
  Double_t* sigma5;
  Double_t* sigma6;
  Double_t* sigma7;
  Double_t* sigma8;
  Double_t* sigma9;
  Double_t* sigma10;
  Double_t* sigma11;
  Int_t crystal;
  Double_t* fClusterTime;      // copied from ED 31/07/05 ;-)
  TFile* fCutFile;                          // File used for TAPS t-o-f cut
  TCutG* fProtonCut;                        // This cut
  Int_t ntaps;
  Int_t particle;
  Int_t fNtapsneutron;
  Int_t fNtapsproton;
  Int_t fNtapselectron;
  Int_t fNtapsgamma;
  Int_t fNtapsrootino;
  Double_t PosShift[4];
  Double_t TimeShift;
  Int_t* fSize;
  Double_t STesting[4];
  UInt_t* fVetoIndex;
  TA2Particle* fParticles;               // use this to store particle-information

public:
  TA2TAPS(const char*, TA2System*);  // pass ptr to analyser
  virtual ~TA2TAPS();                  // destructor
  virtual void PostInit();                    // some setup after parms read in
  virtual TA2DataManager* CreateChild(const char*, Int_t);
  virtual void LoadVariable();              // display setup
  virtual void Reconstruct();                // event by event analysis
  virtual void Cleanup();                     // reset at end of event
  virtual void SetConfig(Char_t*, Int_t);   // setup decode in implement
  virtual void ReadDecoded();            // read back previous analysis

  TA2Particle* GetParticles(){ return fParticles; }
  TA2Particle GetParticles(Int_t index){ return fParticles[index]; }
  UInt_t GetNparticles(){ return fNparticle; }
  UInt_t GetNParticles(){ return fNparticle; }
  void SetParticleInfo(Int_t pNum);

  // in your analyis you have to cast the "Taps-Apparatus" to get this functions
  TA2TAPS_BaF2* GetBaF2() {return fBaF2; }
  Int_t* GetVetos() { return RealVetoHits; }                           // Transform Vetoinformation in correct order
  Double_t* GetClusterTime(){ return fClusterTime; }

  Bool_t* GetVetoInfo() {return isVCharged;}
  Bool_t* GetTofInfo() {return IsProton;}
  Int_t* GetPsaInfo() {return psaInfo;}

  ClassDef(TA2TAPS,1)
};

//---------------------------------------------------------------------------

inline void TA2TAPS::ReadDecoded()
{

}

/*inline void TA2TAPS::ReadDecoded()
{
  // Read back "decoded" data for the Veto array
  // In this case produced by the GEANT3 CB simulation
  Int_t fNvetohits = *(Int_t*)(fEvent[EI_ntaps]);
  Float_t* energy = (Float_t*)(fEvent[EI_eveto]);
  Int_t* index = (Int_t*)(fEvent[EI_iveto]);
  Int_t vi;
  for(vi=0; vi<fNvetohits; vi++ )
  {
      printf("Veto %i fired!\n", *index);
      if (*energy++ > 0) RealVetoHits[vi] = *index++;
  }
  RealVetoHits[vi] = EBufferEnd;
}*/

//-----------------------------------------------------------------------------

inline void TA2TAPS::SetParticleInfo(Int_t pNum)
{
  UInt_t* id_clBaF2 = fBaF2->GetClustHit();                //central detector array
  HitCluster_t** clBaF2 = fBaF2->GetCluster();             //cluster array

  fParticles[pNum].Reset();
  fParticles[pNum].SetP4(fP4[pNum]);
  fParticles[pNum].SetTime(fBaF2->GetTime(id_clBaF2[pNum]));
  fParticles[pNum].SetClusterSize(fSize[pNum]);
  fParticles[pNum].SetParticleIDA(fPDG_ID[pNum]);
  fParticles[pNum].SetCentralIndex(id_clBaF2[pNum]);
  fParticles[pNum].SetApparatus(EAppTAPS);
  fParticles[pNum].SetDetector(EDetBaF2);

  if(fVetoIndex[pNum]!=EBufferEnd)
  {
    fParticles[pNum].SetVetoIndex(fVetoIndex[pNum]);
    fParticles[pNum].SetDetector(EDetVeto);
  }
  if(fBaF2->GetMaxSGElements())
    fParticles[pNum].SetPSAShort(fBaF2->GetSGEnergy(id_clBaF2[pNum]));

  fParticles[pNum].SetSigmaE(fBaF2->GetEnergyResolution(clBaF2[id_clBaF2[pNum]]->GetEnergy()));
  fParticles[pNum].SetSigmaPhi(fBaF2->GetPhiResolution());
  fParticles[pNum].SetSigmaTheta(fBaF2->GetThetaResolution());
  fParticles[pNum].SetSigmaTime(fBaF2->GetTimeResolution());
}

//-----------------------------------------------------------------------------

#endif
