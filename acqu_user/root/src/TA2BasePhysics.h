#ifndef __TA2BasePhysics_h__
#define __TA2BasePhysics_h__

#include "TAcquRoot.h"
#include "TA2Analysis.h"
#include "TA2AccessSQL.h"
#include "TA2Particle.h"
#include "TA2Tagger.h"
#include "TA2CrystalBall.h"
#include "TA2CB.h"
#include "TA2CentralApparatus.h"
#include "TA2Taps.h"
#include "TA2Event.h"
#include "TFile.h"

//Defines the number of separate time windows (for prompt/random hits)
#define NWINDOW 11

class TA2Apparatus;

//-----------------------------------------------------------------------------

class TA2BasePhysics : public TA2AccessSQL
{
  protected:
    //ROOT file output (carrying TA2Particles)
    TFile* EventFile;
    TTree* EventTree;
    TA2Event Event;
    Bool_t bDoROOT;
    Char_t OutputFilename[256];

    //Apparati and detectors
    //TA2Tagger* Tagger;        //Tagger apparatus													//defined in AccessSQL
    //TA2CentralApparatus* CB3; //CB apparatus (NaI, PID, MWPC; variant 2: TA2CentralApparatus)		//defined in AccessSQL
    //TA2Taps* TAPS;            //TAPS apparatus (BaF2, Veto)										//defined in AccessSQL
    //TA2CalArray* NaI;         //NaI array from CB													//defined in AccessSQL
    //TA2TAPS_BaF2* BaF2;       //BaF2 array from TAPS												//defined in AccessSQL

    //Trigger stuff
    Bool_t DoTrigger;         //Flag whether trigger information should be used
    Bool_t SimulateTrigger;
    Bool_t trigger;  // flag which indicates whether the event was triggered or not
    Bool_t UseSumModel;
    Int_t UseM[5];
    Double_t ESumThres;       //Threshold for energy sum in CB
    Double_t ESumSigma;       //Deviation for energy sum in CB
    Double_t ESum;            //L1 Energy sum in CB
    Double_t SumModel[2001];  //
    Double_t GainsNaI[720];   //Relative calibration of NaI crystals (used in energy sum calculation)
    Double_t ThresNaI[720];   //Discriminator thresholds for each single CB NaI channel
    Double_t SigmaNaI[720];   //Discriminator deviations for each single CB NaI channel
    Double_t ThresBaF2[438];  //Discriminator thresholds for each single TAPS BaF2 channel
    Double_t SigmaBaF2[438];  //Discriminator deviations for each single TAPS BaF2 channel
    Double_t ScaleNaI;
    Char_t SumFilename[256];  //Filename 
    Char_t NaIFilename[256];  //Filename to single CB NaI channel information
    Char_t BaF2Filename[256]; //Filename to single TAPS BaF2 channel information
    UInt_t Mult;              //L2 Multiplicity
    UInt_t PrescaleM[5];      //Prescaling for multiplicity trigger (M1-M4)
    UInt_t Rings;             //TAPS rings not contributing to trigger
    UInt_t L1Pattern;
    UInt_t L2Pattern;
    UInt_t NewTAPS;

    //Time windows (for prompt/random hits)
    Char_t WindowFilename[256];    //Filename to time window definition text file
    Double_t Window[NWINDOW+1][2]; //Lower [0] and upper [1] value for each time window

    //Incoming tagged photons
    Int_t nBeam;           //Maximum number of tagged beam photons
    Int_t nTagged;         //Number of detected tagged photons
    TA2Particle* Tagged;   //TA2Particle carrying all information about tagged beam photon

    //Outgoing particles (protons, photons, mesons)
    Int_t nBall;           //Maximum number of detected particles in (Crystal) Ball
    Int_t nWall;           //Maximum number of detected particles in (TAPS) Wall
    Int_t nPart;           //Maximum number of detected particles (CB + TAPS)
    Int_t nPhoton;         //Number of detected photons
    Int_t nProton;         //Number of detected protons
    Int_t nPiPlus;         //Number of detected pi+s
    TA2Particle* Photon;   //TA2Particle carrying all information about detected photons (CB/TAPS)
    TA2Particle* Proton;   //TA2Particle carrying all information about detected protons (CB/TAPS)
    TA2Particle* PiPlus;   //TA2Particle carrying all information about detected pi+s (CB/TAPS)

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 protected:
    //Trigger functions
    void TriggerProcessSW();    //Generates the trigger for the MC
    void TriggerProcessHW();    //Generates the trigger from hardware information
    Bool_t TriggerDecode();     //Decodes trigger L1/L2 patterns

  public:
    TA2BasePhysics(const char*, TA2Analysis*);
    virtual ~TA2BasePhysics();
    virtual void LoadVariable();            //Creates histograms
    virtual void SetConfig(Char_t*, Int_t); //Parses general information from configuration file
    virtual void ParseMisc(char* line);     //Parses additional information from configuration file
    virtual void Reconstruct();             //Event reconstruction
    virtual void PostInit();                //Initialisation etc.

  ClassDef(TA2BasePhysics,1)
};

//-----------------------------------------------------------------------------

#endif

