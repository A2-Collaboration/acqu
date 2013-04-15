#ifndef _TA2KFParticle_h
#define _TA2KFParticle_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"

class TA2KFParticle {

 public:
  TA2KFParticle(){fNvar=4;fAlpha.ResizeTo(fNvar,1);fV_Alpha.ResizeTo(fNvar,fNvar);fT.ResizeTo(fNvar,fNvar-1);};
  TA2KFParticle(TLorentzVector lv,Double_t sig_th,Double_t sig_ph,Double_t sig_e);
  TA2KFParticle(TLorentzVector *lv);

  // TA2KFParticle(TA2KFParticle );
  virtual ~TA2KFParticle(){};

  void Set4Vector(TLorentzVector lv){flv=lv;SetAlpha();};
  void SetAlpha();//Always call this via Set4Vector
  void SetVAlpha(TMatrixD V){fV_Alpha=V;};
  void SetResolutions(Double_t sig_th,Double_t sig_ph,Double_t sig_E){ Polar2Cartesian(sig_th,sig_ph,sig_E);};
  TMatrixD GetVAlpha(){return fV_Alpha;};
  TMatrixD GetAlpha(){return fAlpha;};
  TMatrixD GetT();
  TLorentzVector Get4Vector();
  Int_t GetNVar(){return fNvar;};

  TA2KFParticle Add(TA2KFParticle p1); //returns this+p1
  TA2KFParticle Subtract(TA2KFParticle p1);//returns this-p1

  void Polar2Cartesian(Double_t sig_th,Double_t sig_ph,Double_t sig_e);//Calculate the error matrix in cart. from sig_th etc.

 
 private:
  TMatrixD fAlpha;
  TMatrixD fV_Alpha;
  TMatrixD fT;  //Transformation matrix S. polar ->Cartesian
  TLorentzVector flv;
  Int_t fNvar; //Number of variables =4 for 4 vector

  ClassDef(TA2KFParticle, 1);
};

#endif
