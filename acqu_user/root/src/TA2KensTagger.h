//--Author      JRM Annand         4th Feb 2003
//--Rev         JRM Annand...     26th Feb 2003  1st production tagger
//--Rev         K   Livingston... 21th May 2004  Major KL mods
//--Rev         JRM Annand....... 10th Feb 2005  gcc 3.4 compatible
//--Update      JRM Annand....... 21st Apr 2005  IsPrompt(i) etc.,fMaxParticle
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TTA2KensTagger
//
// Mainz photon-tagging spectrometer apparatus class
// Tagger can consist of 1 or 2 ladders (FP, Microscope) plus PbGlass
// (not yet implemened).
// FP and Micro are instances of the TA2KensLadder class
//

// A typical setup file is like this:
//#       NMR field       Beam Energy
//Misc:   1.01455         853.87
//#
//#         Name          Class Name      Setup file
//Detector: Ladder        TA2KensLadder   FPladderScal.cal
//Detector: Micro         TA2KensLadder   Microscope_kl.dat
//#Condition: xxx
//Reconstruct:
//Initialise:
//Display ...etc
//
// If 2 Ladders are used, then the main FP (ladder) must be before the Micro
// in the setup file

#ifndef __TA2KensTagger_h__
#define __TA2KensTagger_h__

#define MAXRUNS 2047

#include "TA2KensLadder.h"
#include "TA2UserControl.h"

//FUCK OFF!!!
#ifndef private
#define private protected
#endif
#include "TA2Tagger.h"
#ifdef private
#undef private
#endif

class TA2KensTagger : public TA2Tagger
{
 protected:
  Bool_t IsSim;
  Bool_t UseCircular;
  Int_t iRun;
  Int_t CircularRuns;
  Char_t RunName[1024];
  Char_t CurrentRun[1024];
  Char_t CircularFile[1024];
  Char_t CircularRun[MAXRUNS+1][256];
  Double_t CircularDeg[MAXRUNS+1];
  Double_t SupWt;

 public:
  TA2KensTagger(const char*, TA2System*); // construct called by TAcquRoot
  virtual ~TA2KensTagger();
  virtual TA2Detector* CreateChild(const char*, Int_t);
  virtual void ParseMisc(char* line);
  virtual void LoadVariable();
  virtual void Reconstruct();
  virtual void SetParticleInfo(Int_t, Int_t);
  Int_t GetNparticles(){ return fNparticle; }
  Int_t GetNParticles(){ return fNparticle; }

  ClassDef(TA2KensTagger,1)
};

//-----------------------------------------------------------------------------

inline void TA2KensTagger::SetParticleInfo(Int_t pNum, Int_t iHit)
{
  fParticles[pNum].Reset();
  fParticles[pNum].SetP4(fP4[pNum]);
  fParticles[pNum].SetTime(fPhotonTime[pNum]);
  fParticles[pNum].SetParticleIDA(kGamma);
  fParticles[pNum].SetSecondID(kGamma);
  fParticles[pNum].SetUnclear(false);
  fParticles[pNum].SetApparatus(EAppTagger);
  fParticles[pNum].SetDetector(EDetLadder);
  fParticles[pNum].SetCentralIndex(iHit);
  fParticles[pNum].SetClusterSize(1);

  fParticles[pNum].SetWindow(EWindowNone);
  for(UInt_t s=0; s<fFP->GetNhitsRand(); s++)
    if(iHit==fFP->GetHitsRand()[s]) fParticles[pNum].SetWindow(EWindowRandom);
  for(UInt_t s=0; s<fFP->GetNhitsPrompt(); s++)
    if(iHit==fFP->GetHitsPrompt()[s]) fParticles[pNum].SetWindow(EWindowPrompt);

  Int_t nValue;
  Int_t nCircPol[2];
  Bool_t bResult;
  Double_t E_e, E_gamma, Pol_e, Pol_gamma;
  fParticles[pNum].SetCircState(ECircNone);
  fParticles[pNum].SetCircDegree(0.0);
  if(fADC)
  {
    //ADC 6 carries helicity information:
    //Bit 0: MAMI source is switching (1) or stable (0)
    //Bit 1: "Random flipper has generated positive helicity" (according to Tagger Manual)
    //Bit 2: "Random flipper has generated negative helicity" (according to Tagger Manual)
    //Bit 3: Positive (1) or negative (0) helicity as reported from MAMI source (maybe other way round...)
    //Bit 4: Positive (0) or negative (1) helicity as reported from MAMI source (maybe other way round...)
    nValue = fADC[6];
    bResult = !(nValue & 0x0001); //Polarisation is good when MAMI source indicates stable state (bit 0 unset)
    if(bResult)
    {
      //Use polarisation information as reported from MAMI source
      nCircPol[0] = ((nValue & 0x0008) >> 3) + 1; //Transforms helicity bit (bit 3) into values of 1 and 2 for neg/pos polarisation
      nCircPol[0] = nCircPol[0] ^ 0x0000003; //Do an XOR with the 2 LSBs: Inversion
      //Use polarisation information as reported from random flipper
      nCircPol[1] = ((nValue & 0x0006) >> 1); //Transforms helicity bits (bits 1 and 2) into values of 1 and 2 for pos/neg polarisation

      //If in doubt, use random flipper information
      if(nCircPol[0]!=nCircPol[1]) nCircPol[0] = nCircPol[1];

      if(nCircPol[0]==1)
        fParticles[pNum].SetCircState(ECircPos);
      else if(nCircPol[0]==2)
        fParticles[pNum].SetCircState(ECircNeg);

      E_e = fBeamEnergy;
      E_gamma = fPhotonEnergy[pNum];
      Pol_e = fBeamPol;
      Pol_gamma = Pol_e * E_gamma*(E_e + (1./3.)*(E_e - E_gamma))
                        / (E_e*E_e + (E_e - E_gamma)*(E_e - E_gamma)
                        - (2./3.)*E_e*(E_e - E_gamma));
      fParticles[pNum].SetCircDegree(Pol_gamma);
    }
  }
  if(IsSim)
  {
    E_e = fBeamEnergy;
    E_gamma = fPhotonEnergy[pNum];
    Pol_e = fBeamPol;
    Pol_gamma = Pol_e * E_gamma*(E_e + (1./3.)*(E_e - E_gamma))
                      / (E_e*E_e + (E_e - E_gamma)*(E_e - E_gamma)
                      - (2./3.)*E_e*(E_e - E_gamma));
    fParticles[pNum].SetCircDegree(Pol_gamma);
  }
  
  fParticles[pNum].SetLinState(ELinNone);
  fParticles[pNum].SetLinDegree(1.0); //IMPROVE ME!
  if(fADC)
    if((fADC[30]==1)||(fADC[30]==2)) fParticles[pNum].SetLinState(fADC[30]);

  fParticles[pNum].SetSigmaTime(0.75);
  fParticles[pNum].SetSigmaPhi(TMath::TwoPi());
  fParticles[pNum].SetSigmaTheta(0.511 / (fBeamEnergy - 0.511));
  fParticles[pNum].SetSigmaE(fFP->GetECalibration()[pNum] - fFP->GetECalibration()[pNum-1]);
}

//-----------------------------------------------------------------------------


#endif
