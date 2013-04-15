#ifndef _TA2CBKinematicFitter_h
#define _TA2CBKinematicFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TA2KFParticle.h"

class TA2CBKinematicFitter {

 public:

  TA2CBKinematicFitter(){};
  TA2CBKinematicFitter(Int_t npart, Int_t ncon, Int_t unk);
  virtual ~TA2CBKinematicFitter(){};

  Int_t Solve(); //do the least squares fit

  //Form the D and d matrixes for the fit
  void AddInvMassConstraint(Double_t Minv);   //based on Invariant mass of added particles 
  void AddSubInvMassConstraint(Int_t Np, Int_t pid[], Double_t Minv);//Add invariant mass constraint to subset of particles
  void AddTotEnergyConstraint(Double_t Etot);  //based on total energy of added particles
  void AddTotMomentumConstraint(TVector3 mom); //based on total 3 momentum of added particles
  void AddSubMissMassConstraint(TLorentzVector Mom, Int_t Np, Int_t pid[], Double_t MissMass); // Based on missing mass of particles in subset

  void AddPosKFParticle(TA2KFParticle kfp);//Add TA2KFParticles to be fitted contribute + to fPtot
  void AddNegKFParticle(TA2KFParticle kfp);//Add TA2KFParticles to be fitted contribute - to fPtot

  TA2KFParticle GetTotalFitParticle();//returns alpha=fPtot and sum of error matrices  from each particle
  TLorentzVector Get4Vector(){return fPtot;}
  TA2KFParticle GetParticle(Int_t ip);
  TA2KFParticle GetInitialParticle(Int_t ip);
  Double_t GetChi2(){return fchi2;};

  Double_t ConfidenceLevel(){ return TMath::Prob(fchi2,fNcon-fNunKnown);}//Note should be Ncon-Nunknowns
  Double_t Pull(Int_t i) {return (fmAlpha0[i][0]-fmAlpha[i][0])/sqrt(fmV_Alpha0[i][i]-fmV_Alpha[i][i]);}

  void ResetConstraints(){fNconi=0;}
  void ResetParticles(){fNpari=0;fNparti=0;fPtot.SetXYZT(0,0,0,0);}
  void ResetMatrices();
  void Reset(){ResetConstraints();ResetParticles();ResetMatrices();}
  void Debug();

 private:
  Int_t fNvar; //Number of variables (=4 for 4 vector)
  Int_t fNpar; //Number of parameters=Npart*fNvar
  Int_t fNpart; //Number of particles
  Int_t fNcon; //Number of constraints
  Int_t fNparti; //count Number of particles added
  Int_t fNpari; //count Number of particles added
  Int_t fNconi; // countNumber of constraints added
  Int_t fNiter; // Number of times Solve has been called
  Int_t fNunKnown; // Number of unknowns
  TMatrixD fmAlpha0; //original parameters
  TMatrixD fmAlpha;  //fitted parameters
  TMatrixD fmV_Alpha0;//Covariance matrix for original parameters
  TMatrixD fmV_Alpha; //Covariance matrix for fitted parameters
  TMatrixD fmD;      //Matrix of constraint derivitives
  TMatrixD fmd;      //Vector of evaluated constraints
  TMatrixD fmlamda;  //Vector of lagrangian multipliers
  TMatrixD fmV_D;    //Covariance matrix of constraints (TO BE INVERTED)
  Double_t fchi2;
  TMatrixD fT;     //Overall transforamtin matrix from Spherical->Cart
 
  TLorentzVector fPtot;

  ClassDef(TA2CBKinematicFitter, 1);

};

#endif
