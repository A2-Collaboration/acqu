/**
 * TA2OnlinePhysics
 *
 * Physics class for simple online histograms used with CheckSystem macro
 */

#include "TA2OnlinePhysics.h"

ClassImp(TA2OnlinePhysics)

//-----------------------------------------------------------------------------

TA2OnlinePhysics::TA2OnlinePhysics(const char* Name, TA2Analysis* Analysis) : TA2AccessSQL(Name, Analysis)
{

}

//-----------------------------------------------------------------------------

TA2OnlinePhysics::~TA2OnlinePhysics()
{
	delete[] particles;
}

//---------------------------------------------------------------------------

void TA2OnlinePhysics::SetConfig(Char_t* line, Int_t key)
{
	// default SetConfig()
	TA2AccessSQL::SetConfig(line, key);
}

//---------------------------------------------------------------------------

void TA2OnlinePhysics::LoadVariable()
{
    //Call default LoadVariable()
    TA2AccessSQL::LoadVariable();

	TA2DataManager::LoadVariable("invMass2g", &invMass2g, EDSingleX);
	TA2DataManager::LoadVariable("invMass2g1p", &invMass2g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass3g", &invMass3g, EDSingleX);
	TA2DataManager::LoadVariable("invMass3g1p", &invMass3g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass6g", &invMass6g, EDSingleX);
	TA2DataManager::LoadVariable("invMass6g1p", &invMass6g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass2CB", &invMass2CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass2CB1TAPS", &invMass2CB1TAPS, EDSingleX);
	TA2DataManager::LoadVariable("invMass3CB", &invMass3CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass3CB1TAPS", &invMass3CB1TAPS, EDSingleX);
	TA2DataManager::LoadVariable("invMass6CB", &invMass6CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass6CB1TAPS", &invMass6CB1TAPS, EDSingleX);
}

//---------------------------------------------------------------------------

void TA2OnlinePhysics::PostInit()
{
	/* Get pointers to the different apparati */
	// Tagger
	fTAG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("EPT");
	if (!fTAG) {
		fTAG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
		if (fTAG)
			std::cout << "Found Tagger" << std::endl;
		else
			PrintError("", "No Tagger found!", EErrFatal);
	} else
		std::cout << "Found EPT" << std::endl;

	// Crystal Ball
	fCB = (TA2CentralApparatus*)((TA2Analysis*)fParent)->GetChild("CB");
	if (fCB)
		std::cout << "Found Crystal Ball" << std::endl;
	else
		PrintError("", "No Crystal Ball found!", EErrFatal);

	// TAPS
	fTAPS = (TA2Taps*)((TA2Analysis*)fParent)->GetChild("TAPS");
	if (fTAPS)
		std::cout << "Found TAPS" << std::endl;
	else
		std::cout << "TAPS will not be included in the analysis" << std::endl;

	// Get maximum particles
	maxBeamPhotons = fTAG->GetMaxParticle();
	maxParticlesCB = fCB->GetMaxParticle();
	maxParticlesTAPS = fTAPS ? fTAPS->GetMaxParticle() : 0;
	maxParticles = maxParticlesCB + maxParticlesTAPS;

	particles = new TA2Particle[maxParticles];

	// Call default PostInit()
	// Important!: Call this AFTER the previous arrays are defined, otherwise LoadVariable gets called and all arrays are still NULL
	TA2AccessSQL::PostInit();
}

//-----------------------------------------------------------------------------

/* Main method for applying the physics related stuff */
void TA2OnlinePhysics::Reconstruct()
{
	//Perform basic physics tasks
	TA2AccessSQL::Reconstruct();

	//Initialise array counters
	VarInit();

	/* Get the number of particles from the detectors */
	nBeamPhotons = fTAG->GetNparticle();
	nParticlesCB = fCB->GetNparticle();
	nParticlesTAPS = fTAPS ? fTAPS->GetNparticle() : 0;
	nParticles = nParticlesCB + nParticlesTAPS;

	GetParticles();

	// fast check for multiple gamma/CB events with/out proton/TAPS-Cluster
	int nPhoton = 0;
	bool hasProton = false;
	TA2Particle photons[maxParticles];
	for (unsigned int i = 0; i < nParticles; i++) {
		if (particles[i].GetParticleID() == kGamma)
			photons[nPhoton++] = particles[i];
		else if (particles[i].GetParticleID() == kProton && !hasProton)
			hasProton = true;
	}

	/* now sort the events */
	// use number of identified photons
	TLorentzVector tmp(0., 0., 0., 0.);
	switch (nPhoton) {
	case 2:
		invMass2g = (photons[0].GetP4()+photons[1].GetP4()).M();
		if (hasProton)
			invMass2g1p = invMass2g;
		break;
	case 3:
		for (int i = 0; i < 3; i++)
			tmp += photons[i].GetP4();
		invMass3g = tmp.M();
		if (hasProton)
			invMass3g1p = invMass3g;
		break;
	case 6:
		for (int i = 0; i < 6; i++)
			tmp += photons[i].GetP4();
		invMass6g = tmp.M();
		if (hasProton)
			invMass6g1p = invMass6g;
		break;
	default:
		break;
	}
	// use number of clusters in CB and TAPS
	tmp.SetXYZT(0., 0., 0., 0.);
	switch (nParticlesCB) {
	case 2:
		invMass2CB = (particles[0].GetP4()+particles[1].GetP4()).M();
		if (nParticlesTAPS == 1)
			invMass2CB1TAPS = invMass2CB;
		break;
	case 3:
		for (int i = 0; i < 3; i++)
			tmp += particles[i].GetP4();
		invMass3CB = tmp.M();
		if (nParticlesTAPS == 1)
			invMass3CB1TAPS = invMass3CB;
		break;
	case 6:
		for (int i = 0; i < 6; i++)
			tmp += particles[i].GetP4();
		invMass6CB = tmp.M();
		if (nParticlesTAPS == 1)
			invMass6CB1TAPS = invMass6CB;
		break;
	default:
		break;
	}

	//Clean up and set EBufferEnd end markers for arrays
	//TermArrays();
}

//-----------------------------------------------------------------------------

/* Collect all particles from data in one array */
void TA2OnlinePhysics::GetParticles()
{
	nParticles = 0;
	nParticlesCB = 0;
	nParticlesTAPS = 0;

	// Collect detected particles from CB
	for (Int_t nCB = 0; nCB < fCB->GetNparticle(); nCB++) {
		particles[nParticles++] = fCB->GetParticles(nCB);  // copy to TA2Particle array
		nParticlesCB++;
	}

	// Collect detected particles from TAPS, if available
	if (fTAPS)
		for (Int_t nTAPS = 0; nTAPS < fTAPS->GetNparticle(); nTAPS++) {
			particles[nParticles++] = fTAPS->GetParticles(nTAPS);  // copy to TA2Particle array
			nParticlesTAPS++;
		}
}

//-----------------------------------------------------------------------------

void TA2OnlinePhysics::VarInit()
{
	//Initally, set single value histograms to EBufferEnd (this value will not be plotted)
	invMass2g = EBufferEnd;
	invMass2g1p = EBufferEnd;
	invMass3g = EBufferEnd;
	invMass3g1p = EBufferEnd;
	invMass6g = EBufferEnd;
	invMass6g1p = EBufferEnd;
	invMass2CB = EBufferEnd;
	invMass2CB1TAPS = EBufferEnd;
	invMass3CB = EBufferEnd;
	invMass3CB1TAPS = EBufferEnd;
	invMass6CB = EBufferEnd;
	invMass6CB1TAPS = EBufferEnd;
}

//-----------------------------------------------------------------------------

/* Terminate all arrays with EBufferEnd markers at the end of the Reconstruct() method */
/*void TA2OnlinePhysics::TermArrays()
{

}*/

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void TA2OnlinePhysics::ParseMisc(char* line)
{
	TA2AccessSQL::ParseMisc(line);
}

