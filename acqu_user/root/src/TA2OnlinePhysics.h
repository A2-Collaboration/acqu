#ifndef __TA2OnlinePhysics_h__
#define __TA2OnlinePhysics_h__

#include "TA2AccessSQL.h"

class TA2Apparatus;

class TA2OnlinePhysics : public TA2AccessSQL
{
  protected:
	// Apparati for the different used detectors
	TA2Tagger* fTAG;  // Tagger
	TA2CentralApparatus* fCB;  // Crystal Ball
	TA2Taps* fTAPS;  // TAPS
	TA2Ladder* fFPD;  // Ladder
	TA2Apparatus* fPID;
	TA2CalArray* fNaI;
	TA2TAPS_BaF2* fBaF2;
	TA2Apparatus* fVeto;

	TA2Particle* particles;
	UInt_t nParticles, maxParticles;
	UInt_t nParticlesCB, nParticlesTAPS, maxParticlesCB, maxParticlesTAPS;
	UInt_t nBeamPhotons, maxBeamPhotons;

	// quantities to be plotted
	Double_t invMass2g, invMass2g1p, invMass3g, invMass3g1p, invMass6g, invMass6g1p, 
			invMass2CB, invMass2CB1TAPS, invMass3CB, invMass3CB1TAPS, invMass6CB, invMass6CB1TAPS;

	//General functions
	void VarInit();  //Clear all used variables
	//void TermArrays();  //Terminate arrays with EBufferEnd markers
	// Get particles for analysis independent of previous identification
	void GetParticles();

  public:
	TA2OnlinePhysics(const char*, TA2Analysis*);
	virtual ~TA2OnlinePhysics();
	virtual void LoadVariable();  //Creates histograms
	virtual void SetConfig(Char_t*, Int_t);  //Parses general information from configuration file
	virtual void ParseMisc(char* line);  //Parses additional information from configuration file
	virtual void Reconstruct();  //Event reconstruction
	virtual void PostInit();  //Initialisation etc.

  ClassDef(TA2OnlinePhysics, 1)
};

//-----------------------------------------------------------------------------

#endif

