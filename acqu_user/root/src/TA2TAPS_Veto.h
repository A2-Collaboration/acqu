//--Author	S Schumann   10th Jun 2009
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TAPS_Veto
//
// TAPS Veto detectors
// Plastic scintillator tiles in front of TAPS crystals
//

#ifndef __TA2TAPS_Veto_h__
#define __TA2TAPS_Veto_h__

#include "MCBranchID.h"
#include "TA2Detector.h"

class TA2TAPS_Veto : public TA2Detector
{
 private:
 protected:
  UInt_t fUseSigmaEnergy; //Use energy resoution smearing for MC?
  UInt_t fUseSigmaTime;   //Use time resolution smearing in MC?
  Double_t fSigmaEnergy;  //Sigma for energy resolution
  Double_t fSigmaTime;    //Sigma for time resolution
  void Do0Rings();        //read back previous analysis
  void Do1Rings();        //read back previous analysis
  void Do2Rings();        //read back previous analysis
 public:
  TA2TAPS_Veto(const char*, TA2System*); //Normal use
  virtual ~TA2TAPS_Veto();
  virtual void LoadVariable();           //display/cut setup
  virtual void SetConfig(char*, int);    //read in specific parameters
  virtual void Decode();                 //hits -> energy procedure
  virtual void SaveDecoded();            //save local analysis
  virtual void ReadDecoded();            //read back previous analysis

  ClassDef(TA2TAPS_Veto,1)
};

//---------------------------------------------------------------------------

inline void TA2TAPS_Veto::ReadDecoded()
{
  switch(fNelem)
  {
   case 402: //1 PbWO4 ring: 378 BaF2 + 24 PbWO4
    Do1Rings();
    break;
   case 438: //2 PbWO4 rings: 366 BaF2 + 72 PbWO4
    Do2Rings();
    break;
   default:  //No PbWO4 rings: 384 or 510 BaF2, depending on TAPS setup
    Do0Rings();
    break;
  }
}

//---------------------------------------------------------------------------

inline void TA2TAPS_Veto::Do0Rings()
{
  //Read back MC simulation results
  Double_t E = 0.0;
  Double_t T = 0.0;
  Double_t Lo, Hi;
  Double_t total = 0.0;
  UInt_t i,k;
  Int_t j;
  Int_t* index;
  Float_t* energy;

  fNhits = *(Int_t*)(fEvent[EI_nvtaps]);
  index = (Int_t*)(fEvent[EI_ivtaps]);
  energy = (Float_t*)(fEvent[EI_evtaps]);

  for(i=0,k=0; i<fNhits; i++)
  {
    j = *index++;
    if(!j) //Is it a real hit
    {
      i--;
      energy++;
      continue;
    }
    j--;

    //Get energy and time information
    if(fIsEnergy)
    {
      E = (*energy++) * 1000.0 * fEnergyScale; //Convert GeV to MeV and correct bad simulation energies
      if(fUseSigmaEnergy) E+=gRandom->Gaus(0.0, E * fSigmaEnergy);
      Lo = fElement[j]->GetEnergyLowThr();
      Hi = fElement[j]->GetEnergyHighThr();
      if((E < Lo) || (E > Hi)) continue; //Energy inside thresh?
      fEnergy[j] = E;
      fEnergyOR[k] = E;
      total+=E;
    }
    if(fIsTime)
    {
      T = 0.0; //NOT IMPLEMENTED IN G4-SIMULATION!
      if(fUseSigmaTime) T+=gRandom->Gaus(0.0, fSigmaTime);
      Lo = fElement[j]->GetTimeLowThr();
      Hi = fElement[j]->GetTimeHighThr();
      if((T < Lo) || (T > Hi)) continue; //Time inside thresh?
      fTime[j] = T;
      fTimeOR[k] = T;
    }
    fHits[k] = j;
    k++;
  }

  fHits[k] = EBufferEnd;
  fEnergyOR[k] = EBufferEnd;
  fTimeOR[k] = EBufferEnd;
  fTotalEnergy = total;
  fNhits = k;
}

//---------------------------------------------------------------------------

inline void TA2TAPS_Veto::Do1Rings()
{
  //Read back MC simulation results
  Double_t E = 0.0;
  Double_t T = 0.0;
  Double_t Lo, Hi;
  Double_t total = 0.0;
  UInt_t i,k;
  Int_t j;
  Int_t l;
  Int_t* index;
  Float_t* energy;

  fNhits = *(Int_t*)(fEvent[EI_nvtaps]);
  index = (Int_t*)(fEvent[EI_ivtaps]);
  energy = (Float_t*)(fEvent[EI_evtaps]);

  for(i=0,k=0; i<fNhits; i++)
  {
    j = *index++;
    if(!j) //Is it a real hit
    {
      i--;
      energy++;
      continue;
    }
    j--;

    //Get energy and time information
    if(fIsEnergy)
    {
      E = (*energy++) * 1000.0 * fEnergyScale; //Convert GeV to MeV and correct bad simulation energies
      if(fUseSigmaEnergy) E+=gRandom->Gaus(0.0, E * fSigmaEnergy);
    }
    if(fIsTime)
    {
      T = 0.0; //NOT IMPLEMENTED IN G4-SIMULATION!
      if(fUseSigmaTime) T+=gRandom->Gaus(0.0, fSigmaTime);
    }

    //For vetos in front of PbWO4 modules, the veto information is replicated to a hit veto in front of every PbWO4 crystal
    if((j % 64)==0) //Module 0 within a 64-crystal block of TAPS is replaced by PbWO4
    {
      for(Int_t t=0; t<4; t++) //Replicate veto information for all crystals within a PbWO4 module
      {
        //Recalculate element number for vetos in front of PbWO4 crystals:
        // a) (j / 64) gives you the TAPS block number (A...F) based on a geometry with 384 (BaF2) crystals
        // b) 67 * (j / 64) gives you the first element of TAPS block number (A...F) based on a geometry with 402 (PbWO4 & BaF2) crystals
        // c) 4 PbWO crystal instead of 1 BaF2 crystal
        // d) (j % 64) gives you the index within one of the 64-crystal blocks of TAPS
        // e) Information is used for all 4 crystals
        l = (67 * (j / 64)) + (4 * (j % 64)) + t;

        if(fIsEnergy)
        {
          Lo = fElement[l]->GetEnergyLowThr();
          Hi = fElement[l]->GetEnergyHighThr();
          if((E < Lo) || (E > Hi)) continue; //Energy inside thresh?
          fEnergy[l] = E;
          fEnergyOR[k] = E;
          total+=E;
        }
        if(fIsTime)
        {
          Lo = fElement[l]->GetTimeLowThr();
          Hi = fElement[l]->GetTimeHighThr();
          if((T < Lo) || (T > Hi)) continue; //Time inside thresh?
          fTime[l] = T;
          fTimeOR[k] = T;
        }
        fHits[k] = l;
        k++;
      }
    }
    //For vetos in front of BaF2 modules, the element numbers need bo be recalculated in order to use agreements between BaF2 and veto
    else
    {
      //Recalculate element number for vetos in front of BaF2 crystals:
      // a) (j / 64) gives you the TAPS block number (A...F) based on a geometry with 384 (BaF2) crystals
      // b) 67 * (j / 64) gives you the first element of TAPS block number (A...F) based on a geometry with 438 (PbWO4 & BaF2) crystals
      // c) 1 BaF2 crystal has been replaced by 4 PbWO4 crystals, so there are 4-1=3 crystals more at the beginning of a block
      // d) (j % 64) gives you the index within one of the 64-crystal blocks of TAPS
      l = (67 * (j / 64)) + 3 + (j % 64);

      if(fIsEnergy)
      {
        Lo = fElement[l]->GetEnergyLowThr();
        Hi = fElement[l]->GetEnergyHighThr();
        if((E < Lo) || (E > Hi)) continue; //Energy inside thresh?
        fEnergy[l] = E;
        fEnergyOR[k] = E;
        total+=E;
      }
      if(fIsTime)
      {
        Lo = fElement[l]->GetTimeLowThr();
        Hi = fElement[l]->GetTimeHighThr();
        if((T < Lo) || (T > Hi)) continue; //Time inside thresh?
        fTime[l] = T;
        fTimeOR[k] = T;
      }
      fHits[k] = l;
      k++;
    }
  }

  fHits[k] = EBufferEnd;
  fEnergyOR[k] = EBufferEnd;
  fTimeOR[k] = EBufferEnd;
  fTotalEnergy = total;
  fNhits = k;
}

//---------------------------------------------------------------------------

inline void TA2TAPS_Veto::Do2Rings()
{
  //Read back MC simulation results
  Double_t E = 0.0;
  Double_t T = 0.0;
  Double_t Lo, Hi;
  Double_t total = 0.0;
  UInt_t i,k;
  Int_t j;
  Int_t l;
  Int_t* index;
  Float_t* energy;

  fNhits = *(Int_t*)(fEvent[EI_nvtaps]);
  index = (Int_t*)(fEvent[EI_ivtaps]);
  energy = (Float_t*)(fEvent[EI_evtaps]);

  for(i=0,k=0; i<fNhits; i++)
  {
    j = *index++;
    if(!j) //Is it a real hit
    {
      i--;
      energy++;
      continue;
    }
    j--;

    //Get energy and time information
    if(fIsEnergy)
    {
      E = (*energy++) * 1000.0 * fEnergyScale; //Convert GeV to MeV and correct bad simulation energies
      if(fUseSigmaEnergy) E+=gRandom->Gaus(0.0, E * fSigmaEnergy);
    }
    if(fIsTime)
    {
      T = 0.0; //NOT IMPLEMENTED IN G4-SIMULATION!
      if(fUseSigmaTime) T+=gRandom->Gaus(0.0, fSigmaTime);
    }

    //For vetos in front of PbWO4 modules, the veto information is replicated to a hit veto in front of every PbWO4 crystal
    if(((j % 64)==0) || ((j % 64)==1) || ((j % 64)==2)) //Modules 0, 1, 2 within a 64-crystal block of TAPS are replaced by PbWO4s
    {
      for(Int_t t=0; t<4; t++) //Replicate veto information for all crystals within a PbWO4 module
      {
        //Recalculate element number for vetos in front of PbWO4 crystals:
        // a) (j / 64) gives you the TAPS block number (A...F) based on a geometry with 384 (BaF2) crystals
        // b) 73 * (j / 64) gives you the first element of TAPS block number (A...F) based on a geometry with 438 (PbWO4 & BaF2) crystals
        // c) 4 PbWO crystal instead of 1 BaF2 crystal
        // d) (j % 64) gives you the index within one of the 64-crystal blocks of TAPS
        // e) Information is used for all 4 crystals
        l = (73 * (j / 64)) + (4 * (j % 64)) + t;

        if(fIsEnergy)
        {
          Lo = fElement[l]->GetEnergyLowThr();
          Hi = fElement[l]->GetEnergyHighThr();
          if((E < Lo) || (E > Hi)) continue; //Energy inside thresh?
          fEnergy[l] = E;
          fEnergyOR[k] = E;
          total+=E;
        }
        if(fIsTime)
        {
          Lo = fElement[l]->GetTimeLowThr();
          Hi = fElement[l]->GetTimeHighThr();
          if((T < Lo) || (T > Hi)) continue; //Time inside thresh?
          fTime[l] = T;
          fTimeOR[k] = T;
        }
        fHits[k] = l;
        k++;
      }
    }
    //For vetos in front of BaF2 modules, the element numbers need bo be recalculated in order to use agreements between BaF2 and veto
    else
    {
      //Recalculate element number for vetos in front of BaF2 crystals:
      // a) (j / 64) gives you the TAPS block number (A...F) based on a geometry with 384 (BaF2) crystals
      // b) 73 * (j / 64) gives you the first element of TAPS block number (A...F) based on a geometry with 438 (PbWO4 & BaF2) crystals
      // c) 3 BaF2 crystals have been replaced by 3*4=12 PbWO4 crystals, so there are 12-3=9 crystals more at the beginning of a block
      // d) (j % 64) gives you the index within one of the 64-crystal blocks of TAPS
      l = (73 * (j / 64)) + 9 + (j % 64);

      if(fIsEnergy)
      {
        Lo = fElement[l]->GetEnergyLowThr();
        Hi = fElement[l]->GetEnergyHighThr();
        if((E < Lo) || (E > Hi)) continue; //Energy inside thresh?
        fEnergy[l] = E;
        fEnergyOR[k] = E;
        total+=E;
      }
      if(fIsTime)
      {
        Lo = fElement[l]->GetTimeLowThr();
        Hi = fElement[l]->GetTimeHighThr();
        if((T < Lo) || (T > Hi)) continue; //Time inside thresh?
        fTime[l] = T;
        fTimeOR[k] = T;
      }
      fHits[k] = l;
      k++;
    }
  }

  fHits[k] = EBufferEnd;
  fEnergyOR[k] = EBufferEnd;
  fTimeOR[k] = EBufferEnd;
  fTotalEnergy = total;
  fNhits = k;
}

//---------------------------------------------------------------------------

#endif
