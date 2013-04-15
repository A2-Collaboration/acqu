//--Author	JRM Annand   20th Dec 2002
//--Rev 	JRM Annand...26th Feb 2003...1st production version
//--Rev 	JRM Annand... 1st Oct 2003...Incorp. TA2DataManager
//--Rev 	JRM Annand...12th Dec 2003...fMulti 2D plots
//--Rev 	JRM Annand...19th Feb 2004...User code
//--Rev 	JRM Annand...12th May 2004...General CB,TAPS etc. use
//--Rev 	JRM Annand...10th Feb 2005...gcc 3.4 compatible
//--Update	JRM Annand...12th Jul 2005...New TAPS detector class (R.Gregor)
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2Calorimeter
//
// For general description of calorimeter-type systems

#ifndef __TA2Calorimeter_h__
#define __TA2Calorimeter_h__

#include "TA2Apparatus.h"              // base apparatus class
#include "TA2ClusterDetector.h"        // calorimeter array

//class TA2ClusterDetector;

class TA2Calorimeter : public TA2Apparatus {
protected:
  TA2ClusterDetector* fCal;            // calorimeter
  TA2Detector* fPID;                   // charged-particle ID
  TA2Detector* fPTracker;              // charged-particle tracker
  TLorentzVector* fP4_Nphoton;         // Total invariant mass N photons
  Double_t* fM_Nphoton;                // Inv mass for N cluster totals

public:
  TA2Calorimeter( const char*, TA2System* ); // pass ptr to analyser
  virtual ~TA2Calorimeter();
  virtual void PostInit();                // some setup after parms read in
  virtual TA2DataManager* CreateChild( const char*, Int_t );
  virtual void LoadVariable(  );          // display setup
  virtual void PhotonReconstruct( );      // photon 4 momenta
  virtual void Reconstruct(){ PhotonReconstruct(); }
  virtual TLorentzVector SetLorentzVector( Double_t, HitCluster_t* );
  virtual void Cleanup();                 // reset at end of event

  TA2ClusterDetector* GetCal(){ return fCal; }
  TA2Detector* GetPID(){ return fPID; }
  TA2Detector* GetPTracker(){ return fPTracker; }
  TLorentzVector* GetP4_Nphoton(){ return fP4_Nphoton; }
  Double_t* GetM_Nphoton(){ return fM_Nphoton; }

  // Root needs this line for incorporation in dictionary
  ClassDef(TA2Calorimeter,1)
};

//-----------------------------------------------------------------------------
inline void TA2Calorimeter::Cleanup( )
{
  // Clear any arrays holding variables
 TA2DataManager::Cleanup();

 if( !fNparticle ) return;
 Int_t i;
 for( i=0; i<fNparticle; i++ ){
   fP4[i].SetXYZT(0.,0.,0.,0.);
   fMinv[i] = 0.0;
 }
 fP4_Nphoton[i-1].SetXYZT(0.,0.,0.,0.);
 fM_Nphoton[i-1] = 0.0;
 fNparticle = 0;
}

//-----------------------------------------------------------------------------
inline void TA2Calorimeter::PhotonReconstruct( )
{
  // Get energies and "center of charge" vectors from calorimeter arrays
  // epitomised by TA2ClusterDetector.
  // Create 4 vectors assuming cluster hits are photons

  fNparticle = fCal->GetNCluster();            // # cluster hits
  if( !fNparticle ) return;                    // nothing to do

  UInt_t* id_clCal = fCal->GetClustHit();      // indices of hit clusters
  HitCluster_t** clCal = fCal->GetCluster();   // -> cluster structs
  HitCluster_t* cl;                            // cluster struct

  fP4tot.SetXYZT(0,0,0,0);                     // init total 4 momentum
  // Loop round cluster hits
  // for each create 4 momentum from cluster energy and direction
  // assuming its a photon
  Int_t i;
  for (i=0; i<fNparticle; i++){
    cl = clCal[id_clCal[i]];                   // ith hit cluster struct
    fPDG_ID[i] = kGamma;
    // Set 4 vector on the basis of ID above
    fParticleID->SetP4( fP4+i,fPDG_ID[i],cl->GetEnergy(),
			cl->GetMeanPosition() );
    fP4tot += fP4[i];                          // accumulate total
    fMinv[i] = fP4[i].M();
  }
  fP4_Nphoton[i - 1]= fP4tot;                  // total 4-mom i-hit event
  fM_Nphoton[i - 1]= fP4tot.M();               // total inv mass i-hit event
  fMtot = fP4tot.M();
}

//-----------------------------------------------------------------------------
inline TLorentzVector 
TA2Calorimeter::SetLorentzVector(Double_t mass, HitCluster_t* cl)
{
  // Get cluster kinetic energy and particle direction
  // Construct 4 momentum using input assumed mass

  Double_t E = cl->GetEnergy() + mass;
  Double_t P = TMath::Sqrt(E*E - mass*mass);
  TLorentzVector p4;
  p4.SetE(E);
  p4.SetVect( (cl->GetMeanPosition())->Unit() * P );
  return p4;
}


#endif
