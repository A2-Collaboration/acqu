#ifndef __TA2SaschaPhysics_h__
#define __TA2SaschaPhysics_h__

#include "TA2BasePhysics.h"
#include "TA2KFParticle.h"
#include "TA2CBKinematicFitter.h"

class TA2Apparatus;

class TA2SaschaPhysics : public TA2BasePhysics
{
	typedef std::map<int, const char*> IntCharMap;
	typedef std::pair<int, const char*> ICPair;
	typedef std::map<int, const char*>::iterator ICIter;

  protected:
	// some general constants
	static const int N_FINAL_STATE = 3;
	static const double MASS_ETAPRIME = 957.78;
	static const double MASS_PROTON = 938.272;
	static const double MASS_ETA = 547.862;
	static const double MASS_PIZERO = 134.9766;
	static const double R2D = 180./3.14159265359;
	// for examining spectra in degree-steps
	static const int MAX_STEPS = 20;
	// number of prompt and random windows
	static const int N_WINDOWS = 2;//3
	// at the moment I only use one random array, for this define PROMPT and RANDOM
	static const int PROMPT = 0;
	static const int RANDOM = 1;

	// enum for numbering of different cuts and cut combinations
	enum cuts {
		protE,
		copl,
		balance,
		dAlpha,
		missM,
		invM,
		copl_balance,
		balance_missM,
		protE_copl,
		copl_missM,
		missM_invM,
		balance_dAlpha,
		copl_dAlpha,
		protE_copl_balance,
		copl_balance_dAlpha,
		all,
		unknown
	};

	// Declare a map which is used to easily manage the used cuts and their combinations
	IntCharMap cut_names;

	// Apparati for the different used detectors
	TA2Tagger* fTAG;  // Tagger
	TA2CentralApparatus* fCB;  // Crystal Ball
	TA2Taps* fTAPS;  // TAPS
	TA2Ladder* fFPD;  // Ladder
	TA2Apparatus* fPID;
	TA2CalArray* fNaI;
	TA2TAPS_BaF2* fBaF2;
	TA2Apparatus* fVeto;

	// Flag for checking what data type will be processed
	Bool_t MC;

	// Print speed information if specified in config file
	bool speed_info;

	// Print some output for debugging if set to true via config
	bool dbg;

	// Kinematic fit
	TA2CBKinematicFitter* KinFitter;
	TA2KFParticle KFPhoton[2];
	TA2KFParticle KFMeson;

	// Prompt and random time window values
	Double_t promptLow, promptHigh;
	Double_t random1Low, random1High, random2Low, random2High;
	Double_t promptRandomRatio;
	// Counters to terminate arrays according to the entries correctly
	UInt_t nPrompt, nRandom;
	// ... and for the different cut combinations
	UInt_t* n_cuts[unknown];  // unknown is the last entry in the cut enum, thus the number of cuts

	// The polygon cut which will be applied to the energy and momentum balance
	TCutG* cutBalance;

	// bools to indicate if the current event passed the cuts
	bool passedProtonEnergy;
	bool passedCoplanarity;
	bool passedBalance;
	bool passedDAlphaProtTAPS;
	bool passedMissMass;
	bool passedInvMass;
	bool passedAllCuts;

	// Stuff used for a kinematic fit
	TA2CBKinematicFitter* KinFit;

	// ...
	TA2Particle* particles;
	UInt_t nParticles, maxParticles;
	UInt_t nParticlesCB, nParticlesTAPS, maxParticlesCB, maxParticlesTAPS;
	UInt_t nBeamPhotons, maxBeamPhotons;

	//just for testing
	TH1F *taggerTime;
	Double_t invMass2g, invMass2g1p, invMass6g, invMass6g1p, invMass2CB, invMass2CB1TAPS, invMass6CB, invMass6CB1TAPS;

	// Variables to store the true information from the simulation
	TLorentzVector trueP4Beam;
	TLorentzVector trueP4Target;
	Int_t trueNPart;
	Int_t* trueIDPart;
	TLorentzVector* trueP4Gamma;
	TLorentzVector* trueP4Elec;
	TLorentzVector* trueP4Posi;
	TLorentzVector* trueP4Prot;
	TLorentzVector* trueP4MuPls;
	TLorentzVector* trueP4MuMns;
	TLorentzVector* trueP4PiPls;
	TLorentzVector* trueP4PiMns;
	Int_t trueNGamma;
	Int_t trueNPosi;
	Int_t trueNElec;
	Int_t trueNProt;
	Int_t trueNMuPls;
	Int_t trueNMuMns;
	Int_t trueNPiPls;
	Int_t trueNPiMns;
	Double_t* trueProtMass;
	Double_t* trueProtEnergy;
	Double_t* trueProtTheta;
	Double_t* trueProtPhi;
	Double_t* trueGammaMass;
	Double_t* trueGammaEnergy;
	Double_t* trueGammaTheta;
	Double_t* trueGammaPhi;
	Double_t* truePosiMass;
	Double_t* truePosiEnergy;
	Double_t* truePosiTheta;
	Double_t* truePosiPhi;
	Double_t* trueElecMass;
	Double_t* trueElecEnergy;
	Double_t* trueElecTheta;
	Double_t* trueElecPhi;
	Double_t* trueMuPlsMass;
	Double_t* trueMuPlsEnergy;
	Double_t* trueMuPlsTheta;
	Double_t* trueMuPlsPhi;
	Double_t* trueMuMnsMass;
	Double_t* trueMuMnsEnergy;
	Double_t* trueMuMnsTheta;
	Double_t* trueMuMnsPhi;
	Double_t* truePiPlsMass;
	Double_t* truePiPlsEnergy;
	Double_t* truePiPlsTheta;
	Double_t* truePiPlsPhi;
	Double_t* truePiMnsMass;
	Double_t* truePiMnsEnergy;
	Double_t* truePiMnsTheta;
	Double_t* truePiMnsPhi;

	UInt_t nCharged;
	Double_t invM_allPart, invM_2charged;
	Double_t timeProton, timePhoton, timeFS;
	Double_t* timeLeptons;
	Double_t* timeTagger;
	Double_t* dTimeFS;
	Double_t coplanarity, protEnergyReconstr;

	Double_t* invM_cuts[unknown][N_WINDOWS];
	Double_t* missM_cuts[unknown][N_WINDOWS];
	Double_t* balancePx[N_WINDOWS];
	Double_t* balancePy[N_WINDOWS];
	Double_t* balancePz[N_WINDOWS];
	Double_t* balanceE[N_WINDOWS];
	Double_t* missMass[N_WINDOWS];
	Double_t* protonEnergyExpect[N_WINDOWS];
	Double_t* protDAlphaTAPSCl[N_WINDOWS];

	//General functions
	void VarInit();  //Clear all used variables
	void TermArrays();  //Terminate arrays with EBufferEnd markers
	// Create an array holding all possible permutations of particles
	int GetPermutations(Int_t perm[][N_FINAL_STATE]);
	// Returns the binomial of the two passed integers
	unsigned int binomial(int n, int k);
	// Get (and plot) particle properties from simulation
	void GetTrueBeam();
	void SetTrueParticles();
	void PlotTrueParticles();
	// Get particles for analysis independent of previous identification
	void GetParticles();
	void GetTrueParticles();
	// Methods concerning the analysis
	bool ApplyCuts();
	bool KinematicFit();
	void PlotData();

  public:
	TA2SaschaPhysics(const char*, TA2Analysis*);
	virtual ~TA2SaschaPhysics();
	virtual void LoadVariable();  //Creates histograms
	virtual void SetConfig(Char_t*, Int_t);  //Parses general information from configuration file
	virtual void ParseMisc(char* line);  //Parses additional information from configuration file
	virtual void Reconstruct();  //Event reconstruction
	virtual void PostInit();  //Initialisation etc.

  ClassDef(TA2SaschaPhysics, 1)
};

//-----------------------------------------------------------------------------

#endif

