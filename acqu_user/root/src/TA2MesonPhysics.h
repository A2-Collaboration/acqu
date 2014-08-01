#ifndef __TA2MesonPhysics_h__
#define __TA2MesonPhysics_h__

#include "TA2BasePhysics.h"
#include "TA2KFParticle.h"
#include "TA2CBKinematicFitter.h"

class TA2Apparatus;

class TA2MesonPhysics : public TA2BasePhysics
{
  protected:
	static const double MASS_PIZERO = 134.9766;
	static const double MASS_ETA = 547.862;

    TA2CBKinematicFitter* KinFitter;
    TA2KFParticle KFPhoton[2];
    TA2KFParticle KFMeson;

    //Analysis parameters
    Double_t Threshold; //Energy threshold for single meson prduction reaction
    Double_t MinvLo;    //Lower value for 2gamma invariant mass cut
    Double_t MinvHi;    //Upper value for 2gamma invariant mass cut
    Double_t MmissLo;   //Lower value for missing mass cut
    Double_t MmissHi;   //Upper value for missing mass cut
    Int_t kMeson;       //Particle ID code for meson (pi0=111, eta=221)

    TA2Particle Meson;     //TA2Particle carrying all information about reconstructed mesons
    TA2Particle MesonCorr; //TA2Particle carrying all information about reconstructed mesons with corrected mass
    Int_t nMeson;          //Number of reconstructed mesons

    //Analysis spectra
    Double_t MesonMinv; //Invariant mass of 2gamma

    Double_t* MesonMmiss[NWINDOW+1]; //Missing mass calculated for each tagged beam photon
    Int_t nMesonMmiss[NWINDOW+1];    //Counter for corresponding array

    //Final state particle angular and energy distributions
    //These spectra can have multiple values for each defined time window
    Double_t* CMEProton[NWINDOW+1];      //CM energy of proton (if detected)
    Double_t* CMEMeson[NWINDOW+1];       //CM energy of meson
    Double_t* CMThetaProton[NWINDOW+1];  //CM polar angle of proton (if detected)
    Double_t* CMThetaMeson[NWINDOW+1];   //CM polar angle of meson
    Double_t* CMPhiProton[NWINDOW+1];    //CM azimzuthal angle of proton (if detected)
    Double_t* CMPhiMeson[NWINDOW+1];     //CM azimzuthal angle of meson
    Double_t* CMETotal[NWINDOW+1];       //Total CM energy W
    Int_t nCM[NWINDOW+1];                //Counter for corresponding arrays

    Double_t* LabEProton[NWINDOW+1];     //Lab energy of proton (if detected)
    Double_t* LabEMeson[NWINDOW+1];      //Lab energy of meson
    Double_t* LabThetaProton[NWINDOW+1]; //Lab polar angle of proton (if detected)
    Double_t* LabThetaMeson[NWINDOW+1];  //Lab polar angle of meson
    Double_t* LabPhiProton[NWINDOW+1];   //Lab azimzuthal angle of proton (if detected)
    Double_t* LabPhiMeson[NWINDOW+1];    //Lab azimzuthal angle of meson
    Double_t* LabEBeam[NWINDOW+1];       //Beam energy in Lab
    Int_t* LabChBeam[NWINDOW+1];         //Beam energy channel in Lab
    Int_t nLab[NWINDOW+1];               //Counter for corresponding arrays

    //Some properties of photons and protons as available from TA2Particle
    Double_t PhotonThetaCB[2];    //Polar angle of photons in CB
    Double_t ProtonThetaCB;       //Polar angle of protons in CB
    Double_t PhotonThetaTAPS[2];  //Polar angle of photons in TAPS
    Double_t ProtonThetaTAPS;     //Polar angle of protons in TAPS
    Double_t PhotonEnergyCB[2];   //Energy of photons in CB
    Double_t ProtonEnergyCB;      //Energy of protons in CB
    Double_t PhotonEnergyTAPS[2]; //Energy of photons in TAPS
    Double_t ProtonEnergyTAPS;    //Energy of protons in TAPS
    Int_t PhotonSizeCB[2];        //Cluster size of photons in CB
    Int_t ProtonSizeCB;           //Cluster size of protons in CB
    Int_t PhotonSizeTAPS[2];      //Cluster size of photons in TAPS
    Int_t ProtonSizeTAPS;         //Cluster size of protons in TAPS
    Int_t PhotonCentralCB[2];     //Central detector index of photons in CB
    Int_t ProtonCentralCB;        //Central detector index of protons in CB
    Int_t PhotonCentralTAPS[2];   //Central detector index of photons in TAPS
    Int_t ProtonCentralTAPS;      //Central detector index of protons in TAPS
    Int_t ProtonVetoCB;           //PID index of protons in CB
    Int_t ProtonVetoTAPS;         //Veto index of protons in TAPS
    Double_t ProtonDeltaECB;      //PID energy loss of protons in CB
    Double_t ProtonDeltaETAPS;    //Veto energy loss of protons in TAPS

    Double_t* EventTime; //Event time defined as difference between tagged beam photon time and meson time
    Int_t nEventTime;    //Counter for corresponding array

    //General functions
    void VarInit();                    //Clear all used variables
    void TermArrays();                 //Terminate arrays with EBufferEnd markers
    void ParticleProperties();         //Demonstrate some TA2Particle properties
    void FillLabSpectra(Int_t, Int_t); //Calculate Lab observables
    void FillCMSpectra(Int_t, Int_t);  //Calculate CM observables
    void CorrectMatrix();              //Perform matrix-formalism kinematic fit

  public:
    TA2MesonPhysics(const char*, TA2Analysis*);
    virtual ~TA2MesonPhysics();
    virtual void LoadVariable();            //Creates histograms
    virtual void SetConfig(Char_t*, Int_t); //Parses general information from configuration file
    virtual void ParseMisc(char* line);     //Parses additional information from configuration file
    virtual void Reconstruct();             //Event reconstruction
    virtual void PostInit();                //Initialisation etc.

  ClassDef(TA2MesonPhysics,1)
};

//-----------------------------------------------------------------------------

#endif

