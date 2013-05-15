//--Author	JRM Annand         4th Feb 2003
//--Rev 	JRM Annand...     26th Feb 2003  1st production tagger
//--Rev 	K   Livingston... 21th May 2004  Major KL mods
//--Rev 	JRM Annand....... 10th Feb 2005  gcc 3.4 compatible
//--Rev 	JRM Annand....... 21st Apr 2005  IsPrompt(i) etc.,fMaxParticle
//--Rev 	JRM Annand....... 10th Nov 2006  Eelect from ladder
//--Rev 	JRM Annand....... 13th Dec 2008  Add Ken's TA2LinearPol
//--Rev 	JRM Annand....... 14th Apr 2009  Add some Getters
//--Rev 	JRM Annand....... 21st May 2009  Add TA2GenericDetector to list
//--Rev 	JRM Annand....... 12th Nov 2009  TA2Particle support
//                                               fBadScalerChan init
//--Update	JRM Annand.......  1st Apr 2011  Class variables protected
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2Tagger
//
// Mainz photon-tagging spectrometer apparatus class
// Tagger can consist of 1 or 2 ladders (FP, Microscope) plus PbGlass
// (not yet implemened).
// FP and Micro are instances of the TA2Ladder class
// TA2LinearPol is now incorporated into TA2Tagger. The methods used are as
// coded originally by K.Livingston
//

#ifndef __TA2Tagger_h__
#define __TA2Tagger_h__

#include "TA2Apparatus.h"
#include "TA2Particle.h"
#include "TA2Ladder.h"
#include "TF1.h"

// defs of various goni adcs
// copied from TA2LinearPol
enum { EGoniModeADC=0,	//mode see enum below 
       EGoniStepADC=1,	//step no in scan
       EGoniAOffADC=2,	//Azi offset
       EGoniVOffADC=3,	//Voffset
       EGoniHOffADC=4,	//Hoffset
       EGoniAPosADC=5,	//Azi pos
       EGoniVPosADC=6,	//Vpos
       EGoniHPosADC=7,	//Vpos
       EGoniRadiusADC=8,	//Hpos
       EGoniPlaneADC=9,	//Angle of pol plane
};
enum { EPara	=1,	//defs of goni modes
       EPerp	=2,
       EScan	=100
};
 // allow multiple calls to PostInit()
enum { EInitLevel0,EInitLevel1,EInitLevel2 };
enum { EHelicityADC = 6 };

class TA2Tagger : public TA2Apparatus {

protected:
  Double_t fNMR;                        // NMR measurement of field strength
  Double_t fBeamEnergy;                 // incident electron-beam energy
  Double_t fBeamPol;                    // incident electron-beam polarisation
  Double_t* fPhotonEnergy;              // photon energies
  Double_t* fPhotonTime;                // photon times
  TA2Ladder* fFP;                       // main focal-plane array
  TA2Ladder* fMicro;                    // focal-plane microscope
  TA2Detector* fPbGlass;                // photon monitor
  const Double_t *fFPEcalib;		// FP energy calibration table	
  const Double_t *fMicroEcalib;		// Microscope energy calibration table
  Bool_t fIsFP;
  Bool_t fIsMicro;
  Bool_t fIsPbGlass;
  Bool_t fIsLinPol;
  Bool_t fIsHelicityAmbiguous;
  Int_t fNFPElem;
  Int_t fNMicroElem;
  Int_t fHelicityADC;
  TA2Particle* fParticles;               // to store particle-information
  //
  // Linear polarisation stuff
  //
  Int_t    fInitLevel;		// To allow multiple passes at Init
  Int_t    fBaseADC;		// Index of the 1st ADC with goni / scan data
  UShort_t *fGoniADC;		// Array of goni ADC vals
  Int_t    fScanFlagADC;	// Is scanning on
  Int_t    fStepADC;		// Step no in scan
  //  Int_t    fTaggerChannels;	// no of tagger scalers
  Int_t    fNormChannel;	// Channel to use for normalisation
  Double_t fNormEnergy;		// Photon energy for norm of  "enhancement"
  Double_t *fEnergyCalib;	// Photon energy map
  UInt_t *fScalerCurr;	        // ptr to ladder current scaler buffer
  Char_t* fRefFileName;	        // File name for reference tagger scaler dump
  Double_t *fIncSpectrum;	// Array to hold Inc ref. spectrum
  Double_t *fCohSpectrum;	// Array to hold Coh spectrum
  Double_t *fEnhSpectrum;	// Array to hold Enhancement spectrum
  Bool_t   *fBadScalerChan;	// Hold map of bad channels
  Bool_t   fHaveIncScaler;	// Flag to say we have ref from scalers.
  Bool_t   fDoingScalers;	// Flag ot say if we're handling scalers
  Bool_t   fDoingTDCs;		// Flag to say if we're handling TDCs
  Char_t* fTDCRefFileName;      // File name for reference tagger TDC based.
  Bool_t   *fBadTDCChan;	// Hold map of bad channels
  Int_t    fFillFreq;	   // # scaler reads between fills of TDC based spectra
  Int_t    fFillCounter;	// Count no fills of TDC based spectra
  Double_t *fPrompt;		// for spectrum from prompt peak
  Double_t *fRand;		// for background
  Double_t *fTDCIncSpectrum;	// Array to hold Inc ref. spectrum
  Double_t *fTDCCohSpectrum;	// Array to hold Coh spectrum
  Double_t *fTDCEnhSpectrum;	// Array to hold Enhancement spectrum
  Bool_t   fHaveIncTDC;		// Flag to say we have ref from TDC.
  Double_t fPromptRandRatio;// Ratio of the widths of the Prompt/rand time bins
  Double_t fEdge;		// for position of coherent edge chans
  Double_t fEdgeEnergy;		// for position of coherent edge MeV
  UShort_t fGoniAOff;	    // for the varoius bits of info from the goniometer
  UShort_t fGoniVOff;
  UShort_t fGoniHOff;
  UShort_t fGoniPlane;
  UShort_t fGoniRadius;
  UShort_t fGoniMode;

  Double_t fDInc;	   // dummy variables to set up 1D and 2D scaler hists
  Double_t fDCoh;
  Double_t fDEnh;
  Double_t fDCohPara;
  Double_t fDEnhPara;
  Double_t fDCohPerp;
  Double_t fDEnhPerp;

  Double_t fDTDCInc;	   //dummy variables to set up 1D and 2D tdc hists
  Double_t fDTDCCoh;
  Double_t fDTDCEnh;
  Double_t fDTDCCohPara;
  Double_t fDTDCEnhPara;
  Double_t fDTDCCohPerp;
  Double_t fDTDCEnhPerp;
  Double_t fDEdge;
  Double_t fDEdgeEnergy;		// for position of coherent edge MeV

  Double_t fDScanStep;
  UInt_t   fScalerCount;		//running counter of scaler events

  Double_t *fLadderTime;		//from the ladder
  Int_t *fLadderHits;			//from the ladder
  UInt_t *fLadderWindows;
  Char_t* fEdgeFitName;
  Int_t fEdgeChMin;
  Int_t fEdgeChMax;

  TH1F *fHInc;			//ptrs to the hists if they're defined
  TH1F *fHCoh;
  TH1F *fHEnh;
  TH1F *fHCohPara;
  TH1F *fHEnhPara;
  TH1F *fHCohPerp;
  TH1F *fHEnhPerp;
  TH1F *fHTDCInc;			//ptrs to the hists if they're defined
  TH1F *fHTDCCoh;
  TH1F *fHTDCEnh;
  TH1F *fHTDCCohPara;
  TH1F *fHTDCEnhPara;
  TH1F *fHTDCCohPerp;
  TH1F *fHTDCEnhPerp;
  TH2F *fHScan;
  TH1F *fHEdge;
  TH1F *fHEdgeEnergy;
  TF1*  fEdgeFit;

public:
  TA2Tagger( const Char_t*, TA2System* ); // construct called by TAcquRoot
  virtual ~TA2Tagger();
  virtual TA2Detector* CreateChild( const Char_t*, Int_t );
  virtual void LoadVariable(  );          // display setup
  virtual void PostInitialise( );               // finish off initialisation
  virtual void InitLinPol();
  virtual void InitLinPolHist();
  virtual void SetConfig( Char_t*, Int_t );
  //  virtual void ParseDisplay( Char_t* );
  virtual void Reconstruct( );
  virtual void ReconstructLinPol( );
  virtual void ReconstructMicro( );
  TA2Ladder* GetFP(){ return fFP; }
  TA2Ladder* GetMicro(){ return fMicro; }
  Int_t* GetLadderHits(){ return fLadderHits; }
  Int_t GetLadderHits(Int_t i){ return fLadderHits[i]; }
  Double_t* GetLadderTime(){ return fLadderTime; }
  Double_t GetBeamEnergy(){ return fBeamEnergy; }
  virtual Int_t GetWindow(Int_t i){ return  *(fFP->GetWindows() + i); }
  virtual Bool_t IsAny(Int_t i){
    if( *(fFP->GetWindows() + i) == ELaddAny ) return ETrue;
    else return EFalse;
  }
  virtual Bool_t IsPrompt(Int_t i){
    if( *(fFP->GetWindows() + i) == ELaddPrompt ) return ETrue;
    else return EFalse;
  }
  virtual Bool_t IsRandom(Int_t i){
    if( *(fFP->GetWindows() + i) >= ELaddRand ) return ETrue;
    else return EFalse;
  }
  virtual Bool_t IsHelicityAmbiguous(){ return fIsHelicityAmbiguous; }
  TA2Particle* GetParticles(){ return fParticles; }
  TA2Particle GetParticles(Int_t index){ return fParticles[index]; }
  //	Int_t GetNParticles(){ return fNparticle; }
  void SetParticleInfo( Int_t );

  Bool_t IsLinPol(){ return fIsLinPol; }

  Int_t GetGoniAOff(){return fGoniAOff-10000;}
  Int_t GetGoniVOff(){return fGoniVOff-10000;}
  Int_t GetGoniHOff(){return fGoniHOff-10000;}
  Int_t GetGoniPlane(){return fGoniPlane;}
  Int_t GetGoniRadius(){return (int)fGoniRadius;}
  Int_t GetGoniMode(){return (int)fGoniMode;}
  Double_t GetCohEdge(){return fEdge;}
  Double_t GetCohEdgeEnergy(){return fEdgeEnergy;}
  ClassDef(TA2Tagger,1)
};

//-----------------------------------------------------------------------------
inline void TA2Tagger::SetParticleInfo(Int_t pNum)
{
  // Adapted from Sven Schumann's code in TA2KensTagger
  //
  fParticles[pNum].Reset();
  fParticles[pNum].SetP4(fP4[pNum]);
  fParticles[pNum].SetTime(fPhotonTime[pNum]);
  fParticles[pNum].SetParticleIDA(kGamma);
  fParticles[pNum].SetSecondID(kGamma);
  fParticles[pNum].SetUnclear(false);
  fParticles[pNum].SetApparatus(EAppTagger);
  fParticles[pNum].SetDetector(EDetLadder);
  fParticles[pNum].SetCentralIndex(fFP->GetHits()[pNum]);
  fParticles[pNum].SetClusterSize(1);

  Int_t iHit;
  fParticles[pNum].SetWindow(EWindowNone);
  for(UInt_t t=0; t<fFP->GetNhits(); t++)
  {
    iHit = fFP->GetHits()[t];
    for(UInt_t s=0; s<fFP->GetNhitsRand(); s++)
      if(iHit==fFP->GetHitsRand()[s])
	fParticles[pNum].SetWindow(EWindowRandom);
    for(UInt_t s=0; s<fFP->GetNhitsPrompt(); s++)
      if(iHit==fFP->GetHitsPrompt()[s])
	fParticles[pNum].SetWindow(EWindowPrompt);
  }

  fParticles[pNum].SetCircState(ECircNone);
  fParticles[pNum].SetCircDegree(0.0);
  if(fADC)
  {
    // fADC[fHelicityBit] carries helicity information:
    // fHelicityBit = 6 (2009 experiments)
    // Bit 0: MAMI source is switching (1) or stable (0). Don't use switching
    // Bit 1: Random flipper has generated positive helicity
    // Bit 2: Random flipper has generated negative helicity
    // Bit 3: +ve (1) or -ve (0) helicity return from MAMI (maybe other way round...)
    Int_t nValue = fADC[fHelicityADC] & 0xf;
    Bool_t bResult = !(nValue & 0x0001);
    fIsHelicityAmbiguous = kFALSE;
    if(bResult)
    {
      Double_t Pol_e = fBeamPol;
      switch( nValue ){
      case 4:
	// binary 0100: unambiguous -ve helicity
	fParticles[pNum].SetCircState(ECircNeg);
	break;
      case 10:
	// binary 1010: unambiguous +ve helicity
	fParticles[pNum].SetCircState(ECircPos);
	break;
      case 12:
	// binary 1100: ambiguous flipper/mami status...choose flipper = +ve
	fParticles[pNum].SetCircState(ECircPos);
	fIsHelicityAmbiguous = kTRUE;
	break;
      case 2:
	// binary 0010: ambiguous flipper/mami status...choose flipper = -ve
	fParticles[pNum].SetCircState(ECircNeg);
	fIsHelicityAmbiguous = kTRUE;
	break;
      default:
	// anything else doesn't make any sense
	Pol_e = 0.0;
	fIsHelicityAmbiguous = kTRUE;
	break;
      }
      Double_t E_e = fBeamEnergy;
      Double_t E_gamma = fPhotonEnergy[pNum];
      Double_t Pol_gamma = Pol_e * E_gamma*(E_e + (1.0/3.0)*(E_e - E_gamma)) /
	(E_e*E_e + (E_e - E_gamma)*(E_e - E_gamma) - (2.0/3.0)*E_e*(E_e - E_gamma));
      fParticles[pNum].SetCircDegree(Pol_gamma);
    }
  }

  fParticles[pNum].SetLinState(ELinNone);
  fParticles[pNum].SetLinDegree(1.0);
  if(fADC)
    if((fADC[30]==1)||(fADC[30]==2)) fParticles[pNum].SetLinState(fADC[30]);

  fParticles[pNum].SetSigmaTime(0.75);
  fParticles[pNum].SetSigmaPhi(TMath::TwoPi());
  fParticles[pNum].SetSigmaTheta(0.511 / (fBeamEnergy - 0.511));
  fParticles[pNum].SetSigmaE(fFP->GetECalibration()[pNum] - fFP->GetECalibration()[pNum-1]);
}

#endif
