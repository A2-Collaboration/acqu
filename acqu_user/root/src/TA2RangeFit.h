//--Author      Erik Heid
//--Rev
//--Update......
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2RangeFit
//
// User-defined physics class

#ifndef __TA2RangeFit_h__
#define __TA2RangeFit_h__

// Definitions

#define RADDEG       TMath::RadToDeg()
#define DEGRAD       TMath::DegToRad()

#define MASS_PROTON  938.27203
#define MASS_PION    139.57018

#define min_(a,b)    ((a) <= (b) ? (a) : (b))
#define max_(a,b)    ((a) >= (b) ? (a) : (b))
#define dmin_(a,b)   (double)min_(a,b)

#include "CERNfunctions.h"
#include "TA2Physics.h"

class TA2RangeFit:public TA2Physics {
protected:

// RangeFit structures

  Double_t ener_ran, sigm_ran, enran_min, enran_max, chi2_ran, cum_ran;
  Int_t nst_ran, ier_ran;

  Double_t rad_targ;
  Int_t    med_targ;
  Double_t thick_mwpc;
  Double_t thick_tube;
  Double_t thick_dres;
  Int_t    maxiter;
  Double_t chi_cut[6];

  Double_t totpath[7];
  Int_t    medium[7];
  Double_t edep[7];
  Int_t    ispos[3];

  Double_t epar[33];
  Int_t    np_ran;
  Double_t pr_range[429],                                  /* was [33][13] */
    pi_range[429],                                         /* was [33][13] */
    geka, gekb, a2cf[806],                                 /* was [31][13][2] */
    b2cf[806],                                             /* was [31][13][2] */
    c2cf[806];                                             /* was [31][13][2] */

  Int_t itms;
  Int_t itimes;

  Double_t ekmin, ekmax;
  Int_t    nebin;

  Double_t r_bar_min_2[7], r_bar_max_2[7];
  Double_t r_bar_max[7];

// RangeFit routines

  Int_t range_fit_core( Int_t, Int_t *, Int_t *, Double_t *, Double_t *,
                        Double_t *, Int_t *, Double_t *, Int_t *, Double_t * );
  Int_t    noquench_los( Int_t *, Double_t *, Double_t *, Int_t *, Double_t * );
  Int_t    path_calc( Int_t *, Double_t *, Int_t *, Double_t * );
  Double_t frange( Double_t *, Int_t *, Int_t *, Double_t * );
  Double_t frnoque( Double_t *, Int_t *, Int_t *, Double_t * );
  Double_t quench( Double_t * );
  Double_t deloss( Double_t *, Int_t *, Double_t * );
  Double_t givran( Double_t *, Int_t * );
  Double_t givene( Double_t *, Int_t * );
  Int_t    range_dat( void );
  Int_t    quad2( Double_t *, Double_t *, Double_t *, Double_t *, Double_t *,
                  Double_t *, Double_t *, Double_t *, Double_t * );
  Int_t    range_geo( Int_t *, Double_t *, Double_t *, Double_t * );
  Int_t    pb_geom( Double_t *, Double_t *, Double_t * );
  Int_t    pfit( Double_t *, Double_t *, Double_t *, Int_t *, Int_t *, Double_t *,
                 Double_t *, Double_t *, Int_t *, Double_t * );
  Int_t    gridls( Double_t *, Double_t *, Double_t *, Int_t *, Int_t *, Int_t *,
                   Double_t *, Double_t *, Double_t *, Double_t *, Double_t *,
                   Int_t * );
  Double_t fchisq( Double_t *, Double_t *, Int_t *, Int_t *, Int_t *,
                   Double_t * );
  Double_t fbin( Double_t *, Double_t *, Int_t * );
  Double_t ekprot_max( Double_t * );

  TRandom *RAND;

public:
    TA2RangeFit( const char *, TA2Analysis * );

    virtual ~ TA2RangeFit();
    virtual void LoadVariable();                           // variables for display/cuts
    virtual void PostInit();                               // init after parameter input

  virtual Int_t range_fit( Int_t *, Int_t *, Double_t *, Double_t *, Double_t *,
                           Int_t * );
  virtual Int_t range_los( Int_t *, Double_t *, Double_t *, Int_t *,
                           Double_t * );

  Double_t GetChi2() {
    return chi2_ran;
  };

  Double_t GetEner() {
    return ener_ran;
  };

  Double_t GetCumulus() {
    return cum_ran;
  };

  Int_t GetIerr() {
    return ier_ran;
  };

  virtual void ReadRangeTables();
  virtual TA2DataManager *CreateChild( const char *, Int_t ) {
    return NULL;
  }

  ClassDef( TA2RangeFit, 1 );
};

// -------------------------------------------------------------------------------------------------------

#endif
