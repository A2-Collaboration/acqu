//--Author      Erik Heid
//--Rev
//--Rev......
//--Update...
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
//
// User-defined physics class

#include "TA2RangeFit.h"
#include "TRandom.h"
#include "Math/DistFuncMathMore.h"
#include "Math/DistFuncMathCore.h"

ClassImp( TA2RangeFit )

//-----------------------------------------------------------------------------

  TA2RangeFit::TA2RangeFit( const char *name,
			    TA2Analysis * analysis ):TA2Physics( name,
								 analysis )
{
}

//-----------------------------------------------------------------------------

TA2RangeFit::~TA2RangeFit()
{

  // Free up allocated memory...after checking its allocated
  // detector and cuts lists

}

//---------------------------------------------------------------------------

void TA2RangeFit::PostInit()
{

  // Some further initialisation after all setup parameters read in
  // default Cut setup

  TA2Physics::PostInit();

  // RangeFit initializers

  itms   = 0;
  itimes = 0;

}

//-----------------------------------------------------------------------------

void TA2RangeFit::LoadVariable()
{

  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  TA2Physics::LoadVariable();

  return;

}

//-----------------------------------------------------------------------------

void TA2RangeFit::ReadRangeTables()
{

  Double_t  Params[11];
  char      line[80];

  FILE *RangeFile = fopen( "/home/susanna/acquroot/a2/acqu/data/RangeTables.dat", "r" );

  fscanf( RangeFile, "%s", line );

  for ( Int_t t = 0; t < 3; t++ ) {
    fscanf( RangeFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &Params[0], &Params[1], &Params[2], &Params[3], &Params[4],
            &Params[5], &Params[6], &Params[7], &Params[8], &Params[9],
            &Params[10] );
    for ( Int_t s = 0; s < 11; s++ )
      epar[t * 11 + s] = Params[s];
  }

  for ( Int_t s = 0; s < 13; s++ ) {
    fscanf( RangeFile, "%s", line );
    //  fprintf(fLogStream, " Material ", );

    for ( Int_t t = 0; t < 3; t++ ) {
      fscanf( RangeFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &Params[0], &Params[1], &Params[2], &Params[3], &Params[4],
              &Params[5], &Params[6], &Params[7], &Params[8], &Params[9],
              &Params[10] );
      for ( Int_t u = 0; u < 11; u++ )
        pr_range[s * 33 + t * 11 + u] = Params[u];
    }
  }

  for ( Int_t s = 0; s < 13; s++ ) {
    fscanf( RangeFile, "%s", line );

    for ( Int_t t = 0; t < 3; t++ ) {
      fscanf( RangeFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &Params[0], &Params[1], &Params[2], &Params[3], &Params[4],
              &Params[5], &Params[6], &Params[7], &Params[8], &Params[9],
              &Params[10] );
      for ( Int_t u = 0; u < 11; u++ ) {
        pi_range[s * 33 + t * 11 + u] = Params[u];
      }
    }
  }

  fclose( RangeFile );

}

//-----------------------------------------------------------------------------

// **************************************************************************
// *
// *   R A N G E F I T
// *   ===============
// *
// *  list of routines
// *            range_fit               main fitting routine
// *            range_los               entry for Energy losses in DAPHNE
// *            path_calc               entry calculate mean range
// *            frange                  function to fit
// *            quench                  function for quenching
// *            deloss                  function energy loss calculation
// *            givran                  function energy --> range
// *            givene                  function range --> energy
// *            range_dat               load range-energy data tables
// *            quad2                   parabola through 3 given points
// *            range_geo               calculate path, energy limits
// *            pb_geom                 geometry for lead layers
// *            pfit                    fitting routine
// *            gridls                  least square fit routine
// *            fchisq                  function chi2 calculation
// *            fbin                    function dicotomic search
// *            ekprot_max              function max proton kin energy
// *
// **************************************************************************

Int_t TA2RangeFit::range_fit( Int_t *ipart, Int_t *n, Double_t *xcos,
                              Double_t *y, Double_t *ey, Int_t *itarget )
{

  return range_fit_core( 0, ipart, n, xcos, y, ey, itarget, ( Double_t *) 0,
                         ( Int_t *) 0, ( Double_t *) 0 );
}

Int_t TA2RangeFit::range_los( Int_t *ipart, Double_t *e0_par, Double_t *xcos,
                              Int_t *itarget, Double_t *y )
{

  return range_fit_core( 1, ipart, ( Int_t *) 0, xcos, y, ( Double_t *) 0,
                         itarget, e0_par, ( Int_t *) 0, ( Double_t *) 0 );
}

Int_t TA2RangeFit::noquench_los( Int_t *ipart, Double_t *e0_par,
                                 Double_t *xcos, Int_t *itarget,
                                 Double_t *y )
{

  return range_fit_core( 2, ipart, ( Int_t *) 0, xcos, y, ( Double_t *) 0,
                         itarget, e0_par, ( Int_t *) 0, ( Double_t *) 0 );
}

Int_t TA2RangeFit::path_calc( Int_t *istop, Double_t *xcos, Int_t *itarget,
                              Double_t *xpath )
{

  return range_fit_core( 3, ( Int_t *) 0, ( Int_t *) 0, xcos,
                         ( Double_t *) 0, ( Double_t *) 0, itarget,
                         ( Double_t *) 0, istop, xpath );
}

//-----------------------------------------------------------------------------

Int_t TA2RangeFit::range_fit_core( Int_t level, Int_t *ipart, Int_t *n,
                                   Double_t *xcos, Double_t *y, Double_t *ey,
                                   Int_t *itarget, Double_t *e0_par,
                                   Int_t *istop, Double_t *xpath )
{

  /* Initialized data */

  static Double_t  r_targ_choice[7] = { 21.5, 21.5, 21.5, 21.5, 10., 10., 10. };
  static Double_t  x[6]             = { 1., 2., 3., 4., 5., 6. };

  Int_t     const1  = 1;
  Int_t     const3  = 3;
  Int_t     const10 = 10;
  Double_t  x0      = 0;

  Int_t  ii;
  Int_t  i, j;
  Int_t  n_dof;
  Int_t  n_tmp;
  Int_t  ng;
  Int_t  itr;

  Double_t  chisq, elast, theta;
  Double_t  dstep, dummy, cf;
  Double_t  yf[3];
  Double_t  xg, yg;
  Double_t  par[1], sig[1];
  Double_t  chitest[10];
  Double_t  min;

  TRandom *RAND = new TRandom();

  // *
  // ** input
  // *
  // *   ipart            particle type   (1 = proton / 2 = pion)
  // *   n                stop layer      (1:3)
  // *   xcos(6)          vertex coord + cosinus director
  // *   y(3)             experimental energies deposited in scint
  // *   ey(3)            sigma energies (MeV)
  // *   itarget          = 1 : Hydrogen
  // *                    = 2 : Deuterium
  // *                    = 3 : He-3
  // *                    = 4 : He-4
  // *                    = 5 : Carbon
  // *                    = 6 : CH2
  // *                    = 7 : Butanol
  // ** output common block
  // *
  // *
  // *   ener_ran         kinetic energy  (MeV)
  // *   sigm_ran         error of the fit on energy (MeV)
  // *   enran_min        min energy  for the stop layer
  // *   enran_max        max energy for the stop layer
  // *   chi2_ran         reduced chi-square of the fit
  // *   cum_ran          cumulative
  // *   nst_ran          no. of iterations
  // *   ier_ran          error flag 0:7, combination of:
  // *                     2^0 = chi-square > chi_cut
  // *                     2^1 = don't stop in n-th layer
  // *                     2^2 = for F energy > max
  // *                     2^3 = nst_ran > maxiter
  // *
  // * EXPERIMENT DEPENDENT PARAMETERS
  // *   thick_tube         material between target and A-layer (mwpc, pipe
  // *                    = line, ...) expressed as an equivalent thickness
  // *                      of scintillator (e.t.s.)
  // *   thick_dres         material between A and B layer (scint. dressing)
  // *                      in e.t.s.
  //
  // liquid targets: radius = 21.5 mm
  // solid targets : radius = 10.0 mm

  /* Parameter adjustments */

  if ( y )
    --y;

  if ( ey )
    --ey;

  /* Function Body */

  switch ( level ) {
  case 1:
    goto L_range_los;
  case 2:
    goto L_noquench_los;
  case 3:
    goto L_path_calc;
  }

  /* ************************************************************************** */

  ener_ran = 0.;
  sigm_ran = 0.;
  chi2_ran = 0.;
  nst_ran  = 0;
  ier_ran  = 0;

  /* --------------------------------------------------------------------------- */

  /*                                                     initialize range tables */

  if ( itms == 0 ) {
    range_dat();
    itms = 1;
  }

  np_ran     = *ipart;
  med_targ   = *itarget + 5;
  rad_targ   = r_targ_choice[*itarget - 1];
  thick_tube = 3.5;

  /* butanol target ? new value for the thickness before A-layer */

  if ( *itarget == 7 )
    thick_tube = 12.;

  for ( i = *n + 1; i <= 3; ++i ) {
    y[i]  = 0.;
    ey[i] = .001;
  }

  // ---------------------------------------------------------------------------
  //                                             evaluate guess value for energy
  // paths and upper/lower limit for energy

  range_geo( n, xcos, &enran_min, &enran_max );

  /* 1 more point to constraint stopping energy */

  n_tmp = *n;

  if ( *n < 2 ) {
    n_tmp = *n + 1;

    // take a normal distribution for the additional point
    // with mean_value=0 and stand_dev=0.001
    // this is important to have a correct chi-square distribution

    RAND->Rannor( xg, yg );

    y[n_tmp] = ey[n_tmp] * xg * 0.001;
  }

  /* build a grid of 10 values and find min chi2 */

  par[0] = enran_min;
  dstep  = ( enran_max - enran_min ) / 10.;

  for ( j = 1; j <= 10; ++j ) {
    for ( i = 1; i <= n_tmp; ++i )
      yf[i - 1] = frange( x, &i, &const1, par );

    ii             = n_tmp - 1;
    chitest[j - 1] = fchisq( &y[1], &ey[1], &n_tmp, &ii, &const1, yf );
    par[0]        += dstep;
  }

  ng  = 0;
  min = chitest[0];

  for ( i = 1; i < 10; i++ ) {
    if ( min > chitest[i] ) {
      ng  = i;
      min = chitest[i];
    }
  }

  /* starting energy correspond to min chi2 */

  par[0] = enran_min + ( ng - 1 ) * dstep;

  /* --------------------------------------------------------------------------- */

  /*                                                            Least Square Fit */

  itr    = maxiter;
  dstep /= 3.;

  pfit( x, &y[1], &ey[1], &n_tmp, &const1, par, &cf, sig, &itr, &dstep );

  // ---------------------------------------------------------------------------
  //                                                       set output parameters

  ener_ran = par[0];
  sigm_ran = sig[0];
  chi2_ran = cf;
  nst_ran  = itr;

  // ---------------------------------------------------------------------------
  //                                                        calculate cumulative

  n_dof = n_tmp - 1;
  chisq = chi2_ran * ( Double_t ) n_dof;

  cum_ran = ROOT::Math::chisquared_cdf( chisq, n_dof, x0 );

  // ---------------------------------------------------------------------------
  //                                                            test convergence
  // 1. test chi-square

  if ( chi2_ran > chi_cut[*n - 1] )
    ier_ran = 1;

  /* 2. test stopping layer */

  elast = edep[ispos[*n - 1] - 1];

  if ( elast == 0. ) {
    chi2_ran = 1.e3;
    ier_ran += 2;
  }

  /* 3. test max energy for protons stopping */

  if ( *n == 3 ) {
    if ( *ipart == 1 ) {

      /* proton: find & test upper limit */

      theta     = acos( xcos[5] );
      enran_max = ekprot_max( &theta );

      if ( ener_ran > enran_max )
        ier_ran += 4;

      /* ...	else				! pion: bad in any case */

      /* ...	  ier_ran = ier_ran + 4 */

    }
  }

  /* 4. test number of iterations */

  if ( nst_ran > maxiter )
    ier_ran += 8;

  goto L114;

  /* ************************************************************************** */

 L_range_los:

  // *
  // * utility to calculate energies losses in daphne
  // * input:
  // *    ipart           = particle type 1 = p / 2 = pi
  // *    e0_par          = initial kinetic energy
  // *    xcos(1:6)       = vertex coords + cosinus director
  // *    itarget         = target 1=H/2=D/3=He3/4=He4/5=C/6=CH2/7=but
  // * output:
  // *    y(1:3)          = energies deposited in the different layers (A:F)
  // *
  // **************************************************************************

  if ( itms == 0 ) {
    range_dat();
    itms = 1;
  }

  np_ran     = *ipart;
  med_targ   = *itarget + 5;
  rad_targ   = r_targ_choice[*itarget - 1];
  thick_tube = 3.5;

  /* butanol target ? new value for the thickness before A-layer */

  if ( *itarget == 7 )
    thick_tube = 12.;

  range_geo( &const3, xcos, &enran_min, &enran_max );

  //    fprintf(fLogStream, "const3 %d  xcos %f %f %f %f %f %f  min %f  max %f\n", 
  //            const3, xcos[0], xcos[1], xcos[2], xcos[3], xcos[4], xcos[5], enran_min, enran_max);

  par[0] = *e0_par;
  dummy  = frange( x, &const3, &const1, par );

  for ( i = 1; i <= 3; ++i ) {
    y[i] = edep[ispos[i - 1] - 1];

    //  fprintf(fLogStream, "los: i %d  y[i] %f  par[0] %f  np_ran %d\n", i, y[i], par[0], np_ran);
    //  for (Int_t j = 0; j<7; j++) 
    //    fprintf(fLogStream, " j %d  edep %f \n", j, edep[j]);
    //  fflush(fLogStream);

  }

  goto L114;

  /* ************************************************************************** */

 L_noquench_los:

  // *
  // * utility to calculate energies losses in daphne
  // * input:
  // *    ipart           = particle type 1 = p / 2 = pi
  // *    e0_par          = initial kinetic energy
  // *    xcos(1:6)       = vertex coords + cosinus director
  // *    itarget         = target 1=H/2=D/3=He3/4=He4/5=C/6=CH2/7=but
  // * output:
  // *    y(1:3)          = energies deposited in the different layers (A:F)
  // *
  // **************************************************************************

  if ( itms == 0 ) {
    range_dat();
    itms = 1;
  }

  np_ran     = *ipart;
  med_targ   = *itarget + 5;
  rad_targ   = r_targ_choice[*itarget - 1];
  thick_tube = 3.5;

  /* butanol target ? new value for the thickness before A-layer */

  if ( *itarget == 7 )
    thick_tube = 12.;

  range_geo( &const3, xcos, &enran_min, &enran_max );

  par[0] = *e0_par;
  dummy  = frnoque( x, &const3, &const1, par );

  for ( i = 1; i <= 3; ++i )
    y[i] = edep[ispos[i - 1] - 1];

  goto L114;

  /* ************************************************************************** */

 L_path_calc:

  // *
  // * calculate path corresponding to a given stopping layers
  // * input:
  // *    istop           = stopping layer
  // *    xcos(1:6)       = vertex coords + cosinus director
  // *    itarget         = target 1=H/2=D/3=He3/4=He4/5=C/6=CH2/7=but
  // * output:
  // *    xpath           = path to istop layer (mm of scintillator)
  // *
  // **************************************************************************

  if ( itms == 0 ) {
    range_dat();
    itms = 1;
  }

  np_ran = 1;

  /* p or pi does not change calculation ... */

  med_targ   = *itarget + 5;
  rad_targ   = r_targ_choice[*itarget - 1];
  thick_tube = 3.5;

  /* butanol target ? new value for the thickness before A-layer */

  if ( *itarget == 7 )
    thick_tube = 12.;

  range_geo( istop, xcos, &enran_min, &enran_max );
  *xpath = givran( &enran_min, &const1 );

 L114:
  delete RAND;

  return 0;

}                                                          /* range_fit_core */

//-----------------------------------------------------------------------------

/* ************************************************************************** */

double TA2RangeFit::frange( Double_t *xp, Int_t *kp, Int_t *na,
                            Double_t *par )
{

  Int_t     i, j;
  Double_t  e0, e1, e_save;

  // *
  // * Nov 97 V. 2.0 Ale
  // * Jan 98 V. 3.0 Ale
  // *
  // * calculates the energy deposited in the kp layer
  // * from the initial energy E0
  // *
  // *                   --->thick<--- (couche no. kp)
  // *                        ____
  // *                       |    |
  // *                       |    |
  // *                       |    |
  // *                       |    |
  // *                       |    |
  // *                       |____|
  // *                    E0  ---------------------> R0
  // *                          E1 ----------------> R1
  // *
  // *                    edep = E0-E1
  // *
  // **************************************************************************
  // set upper limit for energy
  
  /* Function Body */

  if ( par[0] > epar[30] )
    par[0] = epar[30];

  /* previous energy? short circuit */

  //    if (par[0] == e_save)
  //      goto L100;

  /* save current energy */

  e0     = par[0];
  e_save = par[0];

  // --------------------------------------------------------------------------
  // loop over all layers

  for ( i = 1; i <= ispos[2]; ++i ) {

    /* calculate energy loss */

    edep[i - 1] = deloss( &e0, &medium[i - 1], &totpath[i - 1] );
    e1          = e0 - edep[i - 1];

    //  fprintf(fLogStream, " frange : i %d  ispos[2] %d  edep[i-1] %f  e0 %f  medium[i-1] %d  totpath[i-1] %f  e1 %f\n",
    //          i, ispos[2], edep[i-1], e0, medium[i-1], totpath[i-1], e1 );
    //      fflush(fLogStream);
    
    if ( e0 < 0. )
      goto L100;

    /* the particle comes to rest */

    if ( e1 <= 0. ) {

      /* proton ? apply quenching only for scint (medium = 1) */

      if ( np_ran == 1 && medium[i - 1] == 1 )
        edep[i - 1] = quench( &e0 ) * e0;

      /* fill next positions with zero */

      for ( j = i + 1; j <= ispos[2]; ++j )
        edep[j - 1] = 0.;

      goto L100;
    }

    /* proton ? apply quenching only for scint */

    if ( np_ran == 1 && medium[i - 1] == 1 )
      edep[i - 1] = quench( &e0 ) * e0 - quench( &e1 ) * e1;

    /* prepare next step: input energy = current output energy */

    e0 = e1;
  }

 L100:
  return ( edep[ispos[*kp - 1] - 1] );

}                                                          /* frange */


/* ************************************************************************** */

double TA2RangeFit::frnoque( Double_t *xp, Int_t *kp, Int_t *na,
                             Double_t *par )
{

  Int_t     i, j;
  Double_t  e0, e1, e_save;

  // *
  // * Nov 97 V. 2.0 Ale
  // * Jan 98 V. 3.0 Ale
  // *
  // * calculates the energy deposited in the kp layer
  // * from the initial energy E0
  // *
  // *                   --->thick<--- (couche no. kp)
  // *                        ____
  // *                       |    |
  // *                       |    |
  // *                       |    |
  // *                       |    |
  // *                       |    |
  // *                       |____|
  // *                    E0  ---------------------> R0
  // *                          E1 ----------------> R1
  // *
  // *                    edep = E0-E1
  // *
  // **************************************************************************
  // set upper limit for energy
  
  /* Function Body */

  if ( par[0] > epar[30] )
    par[0] = epar[30];

  /* previous energy? short circuit */
  
  if ( par[0] == e_save )
    goto L100;

  /* save current energy */

  e0     = par[0];
  e_save = par[0];

  // --------------------------------------------------------------------------
  // loop over all layers

  for ( i = 1; i <= ispos[2]; ++i ) {

    /* calculate energy loss */

    edep[i - 1] = deloss( &e0, &medium[i - 1], &totpath[i - 1] );
    e1          = e0 - edep[i - 1];

    /* the particle comes to rest */

    if ( e1 <= 0. ) {

      // proton ? apply quenching
      //          if (np_ran.eq.1) edep(i) = quench(e0)*e0
      // fill next positions with zero

      for ( j = i + 1; j <= ispos[2]; ++j )
        edep[j - 1] = 0.;

      goto L100;
    }

    // proton ? apply quenching
    //          if (np_ran.eq.1) edep(i) = quench(e0)*e0 - quench(e1)*e1
    // prepare next step: input energy = current output energy

    e0 = e1;
  }

 L100:
  return ( edep[ispos[*kp - 1] - 1] );

}                                                          /* frnoque */


/* ************************************************************************** */

double TA2RangeFit::quench( Double_t *x )
{

  double  ret_val;

  // *
  // * for energies > 80 MeV no quenching (y = x + k)
  // *
  // **************************************************************************
  
  if ( *x < 0. ) {
    fprintf( fLogStream, "Error in function range_fit/quench. Energy < 0\n" );

    return 1.;
  }

  //if (*x < 880.)
  //  quench=0.95*x-8.*(1.-exp(-0.1*(x**0.9)))
  //else
  //  quench=x-11.95414
  //
  //-I ok
  //quench=1.01039*x-20.7547*(1.-exp(-0.0237568*(x**0.9)))
  
  if ( *x < 1.5 )
    ret_val = 0.;
  else {

    /*  quench=0.934893*x-9.506*(1-exp(-0.108431*(x**0.9))) */
    /* 	quench=1.14879*x-34.252*(1-exp(-0.0242684*(x**0.9))) */

    ret_val =
      *x * .978173 - ( 1 -
		       exp( pow( ( double ) ( *x ), .933797 ) *
			    -.0877829 ) ) * 9.63761;
  }

  ret_val /= *x;

  /* 	quench=0. */

  return ret_val;

}                                                          /* quench */

/* ***************************************************************************** */

double TA2RangeFit::deloss( Double_t *epkin, Int_t *imat, Double_t *path )
{

  Int_t     niter;
  Int_t     nx;

  Double_t  erel, rlen, step, step_max, a, b, c, e_res, r_min, r_res;
  Double_t  r_tot, deemax;
  Double_t  destep;
  Double_t  energy, dde, ekf;

  // Jan 98 V. 1.0 Ale
  //
  // Gives energy loss of a particle in a given thickness of material.
  // Divide the path in partial steps to minimize error (see Geant 3.21)
  // P2 interpolation of data table
  // Input:
  //      epkin [MeV]     kinetic energy of the particle
  //      imat            material type
  //                      1=sci/2=Fe/3=Pb/4=Al/5=Si/6=H/7=H2/8=He3/9=He4
  //                      10=C/11=CH2/12=but
  //      path [mm]       thickness to be traversed
  //
  // *****************************************************************************

  deemax = .05;

  /* max energy loss rate */

  erel = 0.;

  /* energy loss */

  rlen = *path;

  /* current distance to the end of the layer */

  energy = *epkin;

  /* current energy */

  niter = 0;

  while ( rlen > 0. ) {
    ++niter;

    /* particle below energy threshold? short circuit */

    if ( energy <= epar[0] + .001 ) {
      erel += energy;

      goto L114;
    }

    /* find energy bin */

    nx = Int_t( geka * log10( energy ) + gekb );

    if ( nx > 31 )
      nx = 31;

    if ( nx <= 0 )
      nx = 1;

    /* empirical test to minimize discontinuities */

    /* add by V.Lisin 25-jul-2005 */

    dde = ( energy - epar[nx - 1] ) / ( epar[nx] - epar[nx - 1] );

    if ( dde < .7 )

      /* Computing MAX */

      nx = max_( ( nx - 1 ), 1 );
    else
      nx = min_( nx, 31 );

    /* recall parameters fo P2 interpolation */

    a = a2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435];
    b = b2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] / a / 2.;
    c = c2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] / a;

    //  fprintf(fLogStream, " a %f  b %f  c %f  nx %d  imat %d\n", a, b, c, nx, *imat);

    /* total range of the particle */

    r_tot = -b + a / fabs( a ) * sqrt( b * b - c + energy / a );

    /* only first loop */

    if ( erel == 0. ) {

      /* the particle comes to rest? all is energy is released */

      if ( r_tot <= *path ) {
        erel = energy;
        goto L114;
      }

      /* the particle is near to be stopped? takes longer steps */

      if ( r_tot <= *path * 1.02 )
        deemax = .3;
    }

    /* -------------------------------------------------------------------------- */

    /*        divide path into partial steps according to deemax = max fractional */
    /*        energy loss in a single step. Fixed empirically by user. */
    /*        Phylosophy: */
    /*        1) det min energy after the step (ekf) */
    /*        2) calculate range correspondig to ekf */
    /*        3) max step = range(energy) - range(ekf) */
    /*        4) step = min(max step, distance to the end of the layer) */

    /* min energy after step */

    ekf = ( 1. - deemax ) * energy;

    /* test limits and correct energy */

    if ( ekf < epar[0] )
      ekf = epar[0];
    else if ( ekf >= epar[30] )
      ekf = epar[30] * .99;

    /* range corresponding to ekf */

    r_min = givran( &ekf, imat );

    /* max step possible according to deemax */

    step_max = r_tot - r_min;

    /* take the current step */

    step = dmin_( rlen, step_max );

    /* distance to the end of the layer */

    rlen -= step;

    /* residual range, residual energy of the particle */

    r_res = r_tot - step;

    if ( deemax < .1 )
      e_res = a * ( c + r_res * ( b * 2. + r_res ) );
    else
      e_res = givene( &r_res, imat );

    /* energy loss in the current step */

    destep = energy - e_res;

    /* total energy loss in the layer */

    erel += destep;

    /* prepare next step: current energy */

    energy = e_res;
  }

 L114:
  return erel;

}                                                          /* deloss */

/* ***************************************************************************** */

double TA2RangeFit::givran( Double_t *energy, Int_t *imat )
{

  Int_t  nx;

  Double_t  a, b, c;
  Double_t  dde;

  // give the range of particle np_ran with given energy in material imat
  // P2 interpolation of data table
  //
  // np_ran : 1 = proton / 2 = pion
  // imat   : 1=Sci / 2=Fe / 3=Pb / 4=Al / 5=Si / 6=H / 7=H2 / 8=He3 / 9=He4
  //          10=C / 11=CH2 / 12=But
  //
  // *****************************************************************************
  // find energy bin

  nx = Int_t( geka * log10( *energy ) + gekb );

  if ( nx > 31 )
    nx = 31;

  /* add by V.Lisin 25-jul-2005 */

  if ( nx <= 0 )
    nx = 1;

  // empirical test to minimize discontinuities
  // add by V.Lisin 25-jul-2005

  dde = ( *energy - epar[nx - 1] ) / ( epar[nx] - epar[nx - 1] );

  if ( dde < .7 )

    /* Computing MAX */

    nx = max_( ( nx - 1 ), 1 );
  else
    nx = min_( nx, 31 );

  /* recall parameters fo P2 interpolation */

  a = a2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435];
  b = b2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] / a / 2.;
  c = c2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] / a;

  /* range of the particle */

  return ( -b + a / fabs( a ) * sqrt( b * b - c + *energy / a ) );

}                                                          /* givran */

/* ***************************************************************************** */

double TA2RangeFit::givene( Double_t *range, Int_t *imat )
{

  Int_t  nx;
  Int_t  var;

  // give the energy of particle np_ran with given range in material imat
  // P2 interpolation of data table
  //
  // np_ran : 1 = proton / 2 = pion
  // imat   :  1=Sci /  2=Fe  /  3=Pb  /  4=Al / 5=Si / 6=H / 7=H2 / 8=He3 / 9=He4
  //          10=C   / 11=CH2 / 12=but / 13=NaI
  //
  // *****************************************************************************
  // find range bin

  var = 33;

  if ( np_ran == 1 )
    nx = Int_t( fbin( &pr_range[*imat * 33 - 33], range, &var ) );

  if ( np_ran == 2 )
    nx = Int_t( fbin( &pi_range[*imat * 33 - 33], range, &var ) );

  /* Computing 2nd power */

  return ( a2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] * ( *range * *range ) +
           b2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] * *range +
           c2cf[nx + ( *imat + np_ran * 13 ) * 31 - 435] );

}                                                          /* givene */

/* ***************************************************************************** */

Int_t TA2RangeFit::range_dat( void )
{

  Int_t  imat;
  Int_t  i1, i2, i3;

  Double_t  a, b, c;
  Double_t  x1, y1, y2, y3, x2, x3;

  // Jan 98 V. 2.0 Ale
  //
  // Range table: energy bins equally spaced with log(energy)
  // energy limits [1 MeV, 1000 MeV]
  // Initialize data table
  // calculate coefficients a,b,c for P2 interpolation: R = aE^2 +bE +c
  // materials :  1=Scint /  2=Fe  /  3=Pb  /  4=Al / 5=Si / 6=LH / 7=LH2 / 8=LHe3 / 9=LHe4
  //             10=C     / 11=CH2 / 12=But / 13=NaI
  //
  // *****************************************************************************
  
  fprintf( fLogStream, "\n ... Initialize Range tables ...\n\n" );

  ReadRangeTables();

  rad_targ   = 0.;
  med_targ   = 0;
  //thick_mwpc = 3.7;
  thick_tube = 3.5;
  thick_dres = 0.8;
  maxiter    = 30;

  chi_cut[0] = 0.;
  chi_cut[1] = 4.; //8.;
  chi_cut[2] = 10.;
  chi_cut[3] = 12.;
  chi_cut[4] = 12.;
  chi_cut[5] = 8.;

  for ( Int_t i = 0; i < 7; i++ ) {
    totpath[i] = 0.;
    edep[i]    = 0.;
  }

  medium[0] = 6;                                           // target (Deut)
  medium[1] = 1;                                           // equiv scint
  medium[2] = 1;                                           // equiv scint for PID
  medium[3] = 1;                                           // equiv scint
  medium[4] = 2;                                           // steel
  medium[5] = 13;                                          // NaI
  medium[6] = 1;                                           // dummy scint

  ispos[0] = 3;
  ispos[1] = 6;
  ispos[2] = 7;

  np_ran = 0;

  for ( imat = 1; imat <= 13; ++imat ) {

    /* built coefficient table */

    for ( i1 = 1; i1 <= 31; ++i1 ) {

      /* take i, i+1, i+2 */

      i2 = i1 + 1;
      i3 = i2 + 1;
      y1 = epar[i1 - 1];
      y2 = epar[i2 - 1];
      y3 = epar[i3 - 1];

      /* proton */

      x1 = pr_range[i1 + imat * 33 - 34];
      x2 = pr_range[i2 + imat * 33 - 34];
      x3 = pr_range[i3 + imat * 33 - 34];

      quad2( &x1, &x2, &x3, &y1, &y2, &y3, &a, &b, &c );

      a2cf[i1 + ( imat + 13 ) * 31 - 435] = a;
      b2cf[i1 + ( imat + 13 ) * 31 - 435] = b;
      c2cf[i1 + ( imat + 13 ) * 31 - 435] = c;

      /* pion */

      x1 = pi_range[i1 + imat * 33 - 34];
      x2 = pi_range[i2 + imat * 33 - 34];
      x3 = pi_range[i3 + imat * 33 - 34];

      quad2( &x1, &x2, &x3, &y1, &y2, &y3, &a, &b, &c );

      a2cf[i1 + ( imat + 26 ) * 31 - 435] = a;
      b2cf[i1 + ( imat + 26 ) * 31 - 435] = b;
      c2cf[i1 + ( imat + 26 ) * 31 - 435] = c;
    }
  }

  /* constants calculation */

  ekmin = epar[0];
  ekmax = epar[30];
  nebin = 30;

  geka = nebin / ( log10( ekmax ) - log10( ekmin ) );
  gekb = 1 - geka * log10( ekmin );

  return 0;

}                                                          /* range_dat */

/* ************************************************************************** */

Int_t TA2RangeFit::quad2( Double_t *x1, Double_t *x2, Double_t *x3,
                          Double_t *y1, Double_t *y2, Double_t *y3,
                          Double_t *a, Double_t *b, Double_t *c )
{

  Double_t  f1, f2, f3;

  // *
  // * Jan 98 V. 1.0 Ale
  // *
  // * coefficients of quadratic function through 3 given points
  // * equation form: y = a*x**2 + b*x + c
  // *
  // **************************************************************************

  f1 = *y1 / ( ( *x1 - *x2 ) * ( *x1 - *x3 ) );
  f2 = *y2 / ( ( *x2 - *x1 ) * ( *x2 - *x3 ) );
  f3 = *y3 / ( ( *x3 - *x1 ) * ( *x3 - *x2 ) );
  *a = f1 + f2 + f3;
  *b = -f1 * ( *x2 + *x3 ) - f2 * ( *x1 + *x3 ) - f3 * ( *x1 + *x2 );
  *c = f1 * *x2 * *x3 + f2 * *x1 * *x3 + f3 * *x1 * *x2;

  return 0;

}                                                          /* quad2 */

/* ************************************************************************** */

Int_t TA2RangeFit::range_geo( Int_t *istop, Double_t *x, Double_t *emin,
                              Double_t *emax )
{

  /* Initialized data */

  // mm
  static Double_t  r_bar_min[7] = { 0., 22.0, 47.0, 92.0, 240.0, 250.0, 660.0 };
  //   static Double_t  thick_mat[7] = { 21.5, 3.5, 2.0, 0.4, 1.6, 407.0, 10.0 };
  static Double_t  thick_mat[7] = { 21.5, 3.5, 4.0, 0.4, 1.6, 407.0, 10.0 }; // PID-II

  Int_t  isec;
  Int_t  imin, imax;
  Int_t  i, j;

  Double_t  phi_diff;
  Double_t  cos_phi_diff, rcur, rnex, c, d, e, f;
  Double_t  a_min, a_max, delta, coord_min[3], coord_max[3], dummy, z_diff;
  Double_t  phi_bar;

  // *
  // * materials:1=sci / 2=Fe / 3=Pb / 4=Al / 5=Si / 6=H / 7=H2 / 8=He3 / 9=He4
  // *
  // * Calculates:
  // * - The path totpath(i) of the particle through daphne
  // * - The energy limits emin, emax of the particle with stopping
  // *   channels = istop. Propagates backwards energy loss calculation and
  // *   takes a safety interval = +/- 10%
  // *
  // **************************************************************************

  // --------------------------------------------------------------------------

  if ( itimes == 0 ) {
    thick_mat[0] = rad_targ;
    thick_mat[1] = thick_tube;
    //  thick_mat[3] = thick_dres;

    for ( i = 0; i < 7; ++i ) {
      r_bar_max[i] = r_bar_min[i] + thick_mat[i];

      /* Computing 2nd power */

      r_bar_min_2[i] = r_bar_min[i] * r_bar_min[i];

      /* Computing 2nd power */

      r_bar_max_2[i] = r_bar_max[i] * r_bar_max[i];
    }

    itimes = 1;
  }

  /* set target material */

  medium[0] = med_targ;

  // --------------------------------------------------------------------------
  //                                 thickness of material seen by the particle
  // Computing 2nd power

  c = x[3] * x[3] + x[4] * x[4];
  d = x[0] * x[3] + x[1] * x[4];

  /* Computing 2nd power */

  e = x[0] * x[0] + x[1] * x[1];
  f = d / c;

  for ( j = 1; j <= ispos[2]; ++j ) {                      // for spherical layers

    if ( j >= ispos[1] - 1 ) {
      c = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
      d = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
      e = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      f = d / c;
    }

    /* special for target */

    if ( j == 1 ) {

      /* Computing 2nd power */

      delta = r_bar_max[0] * r_bar_max[0] - e;

      if ( delta < 0. ) {

	/* vertex outside target, take vertex x=y=0 */

        totpath[0] = r_bar_max[0] / sin( acos( x[5] ) );
        goto L1;
      }
      a_min = 0.;
    }
    else

      /* other layers */

      /* Computing 2nd power */

      a_min = sqrt( ( r_bar_min_2[j - 1] - e ) / c + f * f ) - f;

    /* Computing 2nd power */

    a_max = sqrt( ( r_bar_max_2[j - 1] - e ) / c + f * f ) - f;

    // --------------------------------------------------------------------------
    //                                         in/out coordinates of impact point

    for ( i = 0; i < 3; ++i ) {
      coord_min[i] = a_min * x[i + 3] + x[i];
      coord_max[i] = a_max * x[i + 3] + x[i];
    }

    // --------------------------------------------------------------------------
    //                                                    special for Pb geometry

    //      if (j == 8) {
    //          z_diff = fabs(coord_min[2]) - 304.;
    //
    //          if (z_diff > 0.)
    //              pb_geom(coord_min, &x[1], &z_diff);
    //      }
    //      if (j == 10) {
    //          z_diff = fabs(coord_min[2]) - 317.5;
    //
    //          if (z_diff > 0.)
    //              pb_geom(coord_min, &x[1], &z_diff);
    //      }

    // --------------------------------------------------------------------------
    //                                                             calculate path

    dummy = 0.;

    for ( i = 0; i < 3; ++i )

      /* Computing 2nd power */

      dummy +=
	( coord_max[i] - coord_min[i] ) * ( coord_max[i] - coord_min[i] );

    totpath[j - 1] = sqrt( dummy );

  L1:
    ;
  }

  // --------------------------------------------------------------------------
  //                                          correction for poligonal geometry

  //    phi_bar = atan2(coord_min[1], coord_min[0]);
  //
  //    if (phi_bar < 0.)
  //      phi_bar += 6.283185307;
  //
  //    isec         = (Int_t) (phi_bar / .392799);
  //    phi_diff     = (Double_t) isec * .392799 + .19634954 - phi_bar;
  //    cos_phi_diff = 1. / cos(phi_diff);
  //
  //    for (j = ispos[0]; j <= ispos[2]; ++j)
  //      if (j != 4)
  //          totpath[j - 1] *= cos_phi_diff;

  /* min energy: when particle stops before istop layer */

  imin  = ispos[*istop - 1] - 1;
  *emin = givene( &totpath[imin - 1], &medium[imin - 1] );

  for ( i = imin - 1; i >= 1; --i ) {
    rcur  = givran( emin, &medium[i - 1] );
    rnex  = rcur + totpath[i - 1];
    *emin = givene( &rnex, &medium[i - 1] );
  }

  /* safe value 10% lower */

  *emin *= .9;

  /* max energy for particles going through daphne */

  if ( *istop == 3 )
    *emax = 1e3;
  else {

    /* max energy: when particle dead before istop+1 layer */

    imax  = ispos[*istop] - 1;
    *emax = givene( &totpath[imax - 1], &medium[imax - 1] );

    for ( i = imax - 1; i >= 1; --i ) {
      rcur  = givran( emax, &medium[i - 1] );
      rnex  = rcur + totpath[i - 1];
      *emax = givene( &rnex, &medium[i - 1] );
    }

    /* safe value 10% greater */

    *emax *= 1.1;
  }

  return 0;

}                                                          /* range_geo */

/* ************************************************************************** */

Int_t TA2RangeFit::pb_geom( Double_t *coord_min, Double_t *x,
                            Double_t *z_diff )
{

  /* Initialized data */

  static Double_t  tg_alpha = .00485;

  Double_t  tg_theta, a, theta, z1;

  // **************************************************************************
  // *

  /* Function Body */

  theta    = acos( x[5] );
  tg_theta = fabs( tan( theta ) );

  if ( theta > 1.57 )
    *z_diff = -( *z_diff );

  z1            = tg_alpha * *z_diff / ( tg_theta - tg_alpha );
  coord_min[2] += z1;

  /* z */

  a            = ( coord_min[2] - x[2] ) / x[5];
  coord_min[0] = x[0] + a * x[3];

  /* x */

  coord_min[1] = x[1] + a * x[4];

  /* y */

  return 0;

}                                                          /* pb_geom */

/* ************************************************************************** */

Int_t TA2RangeFit::pfit( Double_t * x, Double_t * y, Double_t * dy, Int_t * n,
                         Int_t * na, Double_t * a, Double_t * cf, Double_t * da,
                         Int_t * itr, Double_t * dstep )
{

  Int_t  i;
  Int_t  nf;
  Int_t  ibound;
  Int_t  itr_max;

  Double_t  dlta[2], dltc;
  Double_t  c1;
  Double_t  yf[3];

  Int_t  const1 = 1;

  // *
  // *
  // * itr in input  = max nr. iter.
  // * itr in output = nr. of iter done
  // *
  // **************************************************************************

  /* Parameter adjustments */

  --dy;
  --y;
  --x;
  --da;
  --a;

  /* Function Body */

  for ( i = 0; i < 3; ++i )
    yf[i] = 0.;

  nf      = *n - *na;
  dltc    = .01;
  itr_max = *itr + 1;
  *itr    = 0;
  c1      = 0.;

  *cf = fchisq( &y[1], &dy[1], n, &nf, &const1, yf );

  for ( i = 0; i < *na; ++i )
    dlta[i] = *dstep;

  ibound = 0;

  /* inserted */

  while ( fabs( c1 - *cf ) >= dltc && *itr != itr_max &&
          ibound == 0 && dlta[0] > .001 ) {

    /* inserted */

    c1 = *cf;
    ++( *itr );

    gridls( &x[1], &y[1], &dy[1], n, na, &const1, dlta, &a[1], &da[1], yf, cf,
            &ibound );

    /* mod */

  }

  return 0;

}                                                          /* pfit */

// **************************************************************************
//
//      bevington program 11-1 gridls modified
//
//      69 bev 04-apr-88 jmh
//
//      *** modified feb/93 by ale to prevent infinite loops ***
//
//         *** modified feb/99 by pedro to prevent inifinite loops ***
//         *** modified jul/99 by pedro to prevent inifinite loops ***
//
// **************************************************************************

Int_t TA2RangeFit::gridls( Double_t *x, Double_t *y, Double_t *dy, Int_t *n,
                           Int_t *na, Int_t *md, Double_t *dlta,
                           Double_t *a, Double_t *da, Double_t *yf,
                           Double_t *cf, Int_t *ibound )
{

  Int_t  i, j;
  Int_t  nd, nf;

  Double_t  s, oldaj, c1, c2, c3;
  Double_t  dlt;

  /* Parameter adjustments */

  --yf;
  --da;
  --a;
  --dlta;
  --dy;
  --y;
  --x;

  /* Function Body */

  nf    = *n - *na;
  oldaj = 0.;

  /* inserted */

  if ( nf <= 0 ) {
    *cf = 0.;

    return 0;
  }

  for ( j = 1; j <= *na; ++j ) {
    for ( i = 1; i <= *n; ++i )
      yf[i] = frange( &x[1], &i, na, &a[1] );

    /* L1: */
    c1  = fchisq( &y[1], &dy[1], n, &nf, md, &yf[1] );
    nd  = 0;
    dlt = dlta[j];

  L2:
    a[j] += dlt;


    /* c old      if(a(j).eq.oldaj) then             ! inserted (feb 99) */

    if ( fabs( a[j] - oldaj ) < .01 ) {

      /* inserted (jul 99) */

      *ibound = 1;

      /* inserted */

      goto L7;

      /* inserted */

    }

    /* inserted */

    oldaj = a[j];

    /* inserted */

    for ( i = 1; i <= *n; ++i )
      yf[i] = frange( &x[1], &i, na, &a[1] );

    /* L3: */
    c2 = fchisq( &y[1], &dy[1], n, &nf, md, &yf[1] );

    if ( c1 < c2 ) {
      dlt   = -dlt;
      a[j] += dlt;
      s     = c1;
      c1    = c2;
      c2    = s;
    }
    else if ( c1 == c2 )
      goto L2;

  L5:
    ++nd;
    a[j] += dlt;

    for ( i = 1; i <= *n; ++i )
      yf[i] = frange( &x[1], &i, na, &a[1] );

    /* L6: */

    c3 = fchisq( &y[1], &dy[1], n, &nf, md, &yf[1] );

    if ( c3 < c2 ) {
      c1 = c2;
      c2 = c3;

      goto L5;
    }
    else if ( c3 == c2 )

      /* modified */

      dlt *= .5;

    /* inserted */

    else
      dlt *= 1. / ( ( c1 - c2 ) / ( c3 - c2 ) + 1. ) + .5;

    a[j]   -= dlt;
    da[j]   = dlta[j] * sqrt( 2. / ( ( Double_t ) nf * ( c3 - c2 * 2. + c1 ) ) );
    dlta[j] = dlta[j] * ( Double_t ) nd / 3.;

  L7:
    ;
  }

  for ( i = 1; i <= *n; ++i )
    yf[i] = frange( &x[1], &i, na, &a[1] );

  /* L8: */

  *cf = fchisq( &y[1], &dy[1], n, &nf, md, &yf[1] );

  return 0;

}                                                          /* gridls */

// **************************************************************************
//
//      bevington program 10-2 fchisq modified
//
//      69 bev 04-apr-88 jmh
//
// **************************************************************************

double TA2RangeFit::fchisq( Double_t *y, Double_t *dy, Int_t *n, Int_t *nf,
                            Int_t *md, Double_t *yf )
{

  Double_t  ret_val;

  Int_t  i;

  Double_t  c;
  Double_t  w;

  /* Parameter adjustments */

  --yf;
  --dy;
  --y;

  /* Function Body */

  c = 0.;

  if ( *nf <= 0 )
    ret_val = 0.;
  else {
    for ( i = 1; i <= *n; ++i ) {
      if ( *md < 0 )
        if ( y[i] != 0. )
          w = 1. / fabs( y[i] );
        else
          w = 1.;
      else if ( *md == 0 )
        w = 1.;
      else if ( dy[i] != 0. )

	/* inserted */

	/* Computing 2nd power */

        w = 1. / ( dy[i] * dy[i] );
      else
        w = 1.;

      /* Computing 2nd power */

      c += w * ( ( y[i] - yf[i] ) * ( y[i] - yf[i] ) );

      /* L1: */
    }
    ret_val = c / ( Double_t ) ( *nf );
  }

  return ret_val;

}                                                          /* fchisq */

/* ************************************************************************** */

double TA2RangeFit::fbin( Double_t *a, Double_t *x, Int_t *nmax )
{

  Int_t  ix, iy, mid;

  // *
  // * find i=1,...,nmax-1 with a(i) < x < a(i+1)
  // *
  // **************************************************************************

  /* Parameter adjustments */

  --a;

  /* Function Body */

  ix = 1;
  iy = *nmax;

  /* search increasing arguments */

 L1:
  mid = ( ix + iy ) / 2;

  if ( *x >= a[mid] )
    goto L2;

  iy = mid;

  goto L3;

 L2:
  ix = mid;

 L3:
  if ( iy - ix > 1 )
    goto L1;

  return ( Double_t ) ix;

}                                                          /* fbin */

/* ************************************************************************** */

double TA2RangeFit::ekprot_max( Double_t *theta )
{

  Double_t  thdeg, pm_max;

  // *
  // * calculate max proton kinetic energy to allow good identification
  // * input theta angle in radiants
  // *
  // radiants --> degrees

  thdeg = *theta * 57.2958;

  /* calculate max proton momentum */

  /* Computing 2nd power */

  pm_max = ( thdeg - 90. ) * ( thdeg - 90. ) * .025017 + 850.;

  /* momentum --> kin. energy */

  /* Computing 2nd power */

  return ( sqrt( pm_max * pm_max + MASS_PROTON * MASS_PROTON ) - MASS_PROTON );

}                                                          /* ekprot_max */
