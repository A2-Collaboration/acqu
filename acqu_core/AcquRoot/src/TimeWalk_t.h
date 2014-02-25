//--Author	JRM Annand   27th May 2004
//--Rev
//--Rev         JRM Annand    8th Sep 2006  Add SetWalk()
//--Update      JRM Annand    3rd Dec 2008  TimeWalkSergey_t, TimeWalkSven_t
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TimeWalk_t
// Correct hit time obtained from TDC for leading-edge walk
// As of AcquRoot 4v3 3 versions provided
// Basic by JRM Annand
// TimeWalkSergey_t Extended by S Prakov
// TimeWalkSven_t Extended by S Schuman
//

#ifndef __TimeWalk_t_h__
#define __TimeWalk_t_h__

//-----------------------------------------------------------------------------
// Basic AcquRoot time-walk correction class
//-----------------------------------------------------------------------------
class TimeWalk_t {
 protected:
  Double_t fRiseTime;                       // pulse rise time (ns)
  Double_t fThreshold;                      // discriminator threshold (MeV)
 public:
  virtual ~TimeWalk_t() {}
  virtual void SetWalk( Double_t* par ){
    fRiseTime = par[0];
    fThreshold = par[1];
  }
  virtual Double_t Tcorr( Double_t energy, Double_t time ){
    // correct the input time on the basis of the input energy
    return ( time  - fRiseTime * TMath::Sqrt( fThreshold/energy ) );
  }
};

//-----------------------------------------------------------------------------
// Extended AcquRoot time-walk correction class by Sergey Prakov
//-----------------------------------------------------------------------------
class TimeWalkSergey_t : public TimeWalk_t {
 private:
  Double_t fShift;     
  Double_t fPower;     
 public:
  virtual ~TimeWalkSergey_t() {}
  virtual void SetWalk( Double_t* par ){
    fRiseTime = par[0];
    fThreshold = par[1];
    fShift = par[2];
    fPower = par[3];
  }
  Double_t Tcorr( Double_t energy, Double_t time ){
    // correct the input time on the basis of the input energy
    return ( time  - ( fRiseTime + fThreshold/TMath::Power(energy+fShift,fPower)));
  }
};

//-----------------------------------------------------------------------------
// Extended AcquRoot time-walk correction class by Sven Schuman
//-----------------------------------------------------------------------------
class TimeWalkSven_t : public TimeWalk_t
{
 private:
  Double_t fPar[6];
 public:
  virtual ~TimeWalkSven_t() {}
  virtual void SetWalk( Double_t* par ){
    for( Int_t i=0; i<6; i++ ) fPar[i] = par[i];
  }
  Double_t GetPar(Int_t t) { return fPar[t]; }

  Double_t Tcorr( Double_t energy, Double_t time ){
    //Shift time to zero, if energy is above cutoff...
    if(energy>fPar[3])
      time = time - fPar[4];
    //...or do walk correction, if energy is below cutoff
    else{
      time -= fPar[0] + fPar[1]*TMath::Power(energy, fPar[2]);
    }
    return time + fPar[5]; //Shift away from 0
  }
};

#endif
