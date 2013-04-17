#ifndef __TA2ExampleRangeFit_h__
#define __TA2ExampleRangeFit_h__

// AcquROOT includes
#include "TA2Physics.h"
#include "TAcquRoot.h"
#include "TA2Tagger.h"
#include "TA2Event.h" 
#include "TA2CrystalBall.h"
#include "TA2Analysis.h"
#include "TA2Physics.h"
#include "TA2Particle.h"
#include "TA2CrystalBall.h"
#include "TA2CB.h"
#include "TA2RangeFit.h"
// Defines the number of separate time windows (for prompt/random hits)
#define OPT_NWIND 3

// AcquROOT classes
class TA2Tagger;
class TA2KensTagger;
class TA2Ladder;
class TA2CentralApparatus;
class TA2CylMwpc;

class TA2ExampleRangeFit : public TA2Physics {
protected:
  
  static const Int_t kNch = 2;
  
  ////////////////////////////
  // Apparatuses & Detectrors
  ////////////////////////////
  

  TA2CrystalBall        *fCBall;        // CB apparatus
  TA2CentralApparatus	*fCB;	        // Central apparatus
  TA2CalArray		*fNaI;	        // NaI array from CB
  TA2CylMwpc		*fMwpc;		// Mwpc
  TA2PlasticPID		*fPid;	        // Pid
  TA2Tagger             *fTagger;       // Tagger apparatus
  TA2RangeFit           *fRangeFit;
  TA2Ladder* fFPD;                      // Focal Plane Detector apparatus
  TLorentzVector* fTAGGp4;

  // Incoming tagged photons
  Int_t nBeam;                             // Max. no. of tagged photons
  const Int_t *fNintersTrue;
  // Time windows (for prompt/random hits)
  Char_t WindowFilename[255];              // Filename to time window definition text file
  Double_t Window[OPT_NWIND+1][2];         // Lower [0] and upper [1] value for each time window
  Int_t evtDaq;

  ////////////////////////////
  // Histograms
  ////////////////////////////
  Int_t fNtracks;
  // Neutral
  Int_t fNne;
  Double_t *fThetaNe;
  Double_t *fPhiNe;
  Double_t *fEne;
  Double_t *fMclNe;
  Double_t *fTne;
  Double_t *fCentralIndexNaINe;
  Int_t controlcheck;
  Int_t gedanken;
  // Charged
  Double_t *fThetaCh;
  Double_t *fThetaChCB;
  Double_t *fPhiCh;
  Double_t *fEch;
  Double_t *fMclCh;
  Double_t *fTch;
  Double_t *fCentralIndexNaICh;
  Double_t *fEpid;
  Double_t *fIhitPid;
  Double_t *fPsVertexX;
  Double_t *fPsVertexY;
  Double_t *fPsVertexZ;
  Double_t *fPsVertexGoodZ;
  Double_t *fChTrackType;
  Double_t *fTotal_Energy;
  Double_t *fTotal_Energy_prot;
  Double_t *fTotal_Energy_pion;

  Int_t fNch;
  Int_t RanErrProt;
  Int_t RanErrPion;

  Double_t *fDzVertex;
  Double_t *fZVertex;
  Int_t fNvertex;
  Double_t *fEmwpc;
  Double_t *fDEmwpc;
  Double_t fM2g;
  Double_t fECB;
  Double_t fEPID;
  Double_t fEPIDs;
  Double_t EneRFProt;
  Double_t EneRFPion;
  Double_t chi2Prot;
  Double_t chi2Pion;
  Double_t chi2ProtPion;
  Double_t chi2PionProt;
  Double_t CumulProt;
  Double_t CumulPion;
  Double_t CumulAllProt;
  Double_t CumulAllPion;
  Double_t PiTheta, PiPhi, PiEnergy, MesonMinv;

  TA2Particle* Photon;                     // TA2Particle carrying all information about detected photons (CB/TAPS)
  TA2Particle Meson;                       // TA2Particle carrying all information about reconstructed mesons


  Double_t EgammaCal[353];


  // Beam Polarization
  static const Int_t fgBeamPolConst = 239;
  // Energy release in PID
  Double_t* fPID_dE[24];                      
  Double_t* fPID_dE_prot[24];                      
  Double_t* fPID_dE_pion[24];                      
  Double_t* fPID_dE_raw[24]; 

  // Multihit special 
  Double_t          *fFPDTimeOR;	   // FPD TimeOR for all the hits multiplicities
  Int_t              fFPDNHits;		   // # FPD hits for the hits multiplicities
  Int_t             *fFPDHits;		   // FPD hits for all the hits multiplicities
 Double_t* EventTime;                     // Event time as diff between tagged beam photon time and meson time
 Int_t nEventTime;                        // Counter for corresponding array

 // For Charged particle energies
 Int_t IndPID[24], IndPID_prot[24], IndPID_pion[24] ;
 Int_t nCB ;
 Int_t nCB_prot;
 Int_t nCB_pion;
 Int_t nCBnoPID ;
 Double_t* fDE_MWPC;         

 // General functions
  void VarInit();                          // Clear all used variables
  void TermArrays();                       // Terminate arrays with EBufferEnd markers
  void ParticleProperties();               // Demonstrate some TA2Particle properties
  Int_t* BeamPol10[OPT_NWIND+1];           // Beam pol (antiparallel) for each tagged beam photon for bgd
  Int_t nBeamPol10[OPT_NWIND+1];           // Counter for corresponding array
  Double_t* Egamma[OPT_NWIND+1];
  Int_t nEgamma[OPT_NWIND+1];
  Double_t* PhiPi[OPT_NWIND+1];
  Int_t nPhiPi[OPT_NWIND+1];
  Double_t* ThetaPi[OPT_NWIND+1];
  Int_t nThetaPi[OPT_NWIND+1];
  Double_t* PhiPi0[OPT_NWIND+1];
  Int_t nPhiPi0[OPT_NWIND+1];
  Double_t* ThetaPi0[OPT_NWIND+1];
  Int_t nThetaPi0[OPT_NWIND+1];
  Double_t* ThetaPi0Kin[OPT_NWIND+1];
  Int_t nThetaPi0Kin[OPT_NWIND+1];
  Double_t* ThetaNeuKin[OPT_NWIND+1];
  Int_t nThetaNeuKin[OPT_NWIND+1];
  Double_t* DiffThetaPi0[OPT_NWIND+1];
  Int_t nDiffThetaNeu[OPT_NWIND+1];
  Double_t* DiffThetaNeu[OPT_NWIND+1];
  Int_t nDiffThetaPi0[OPT_NWIND+1];
  Double_t* EnePi[OPT_NWIND+1];
  Int_t nEnePi[OPT_NWIND+1];
  Double_t* EnePiKin[OPT_NWIND+1];
  Int_t nEnePiKin[OPT_NWIND+1];
  Double_t* EnePi0[OPT_NWIND+1];
  Int_t nEnePi0[OPT_NWIND+1];
  Double_t* EnePi0Kin[OPT_NWIND+1];
  Int_t nEnePi0Kin[OPT_NWIND+1];
  Double_t* EneNeuKin[OPT_NWIND+1];
  Int_t nEneNeuKin[OPT_NWIND+1];
  Double_t* ThetaProton[OPT_NWIND+1];
  Int_t nThetaProton[OPT_NWIND+1];
  Double_t* PhiProton[OPT_NWIND+1];
  Int_t nPhiProton[OPT_NWIND+1];
  Double_t* EneProton[OPT_NWIND+1];
  Int_t nEneProton[OPT_NWIND+1];   
  Double_t* EneProtonKin[OPT_NWIND+1];
  Int_t nEneProtonKin[OPT_NWIND+1];  
  Double_t* KinProton[OPT_NWIND+1];
  Int_t nKinProton[OPT_NWIND+1];
  Double_t* KinPi0[OPT_NWIND+1];
  Int_t nKinPi0[OPT_NWIND+1];
  Double_t* KinNeu[OPT_NWIND+1];
  Int_t nKinNeu[OPT_NWIND+1];
  Double_t* KinPi[OPT_NWIND+1];
  Int_t nKinPi[OPT_NWIND+1];
  Double_t* deltaphi[OPT_NWIND+1];
  Int_t ndeltaphi[OPT_NWIND+1];
  Double_t* EneCBProt[OPT_NWIND+1];
  Int_t nEneCBProt[OPT_NWIND+1];
  Double_t* EneCBKinProt[OPT_NWIND+1];
  Int_t nEneCBKinProt[OPT_NWIND+1];
  Double_t* EnePIDProt[OPT_NWIND+1];
  Int_t nEnePIDProt[OPT_NWIND+1];
  Double_t* EnePIDKinProt[OPT_NWIND+1];
  Int_t nEnePIDKinProt[OPT_NWIND+1];
  Double_t* EnePIDNormProt[OPT_NWIND+1];
  Int_t nEnePIDNormProt[OPT_NWIND+1];
  Double_t* MissingEneProt[OPT_NWIND+1];
  Int_t nMissingEneProt[OPT_NWIND+1]; 
  Double_t* CheckMEProt[OPT_NWIND+1];
  Int_t nCheckMEProt[OPT_NWIND+1];
  Double_t* CheckMEPion[OPT_NWIND+1];
  Int_t nCheckMEPion[OPT_NWIND+1];
  Double_t* EneCBPion[OPT_NWIND+1];
  Int_t nEneCBPion[OPT_NWIND+1];
  Double_t* EneCBKinPion[OPT_NWIND+1];
  Int_t nEneCBKinPion[OPT_NWIND+1];
  Double_t* EneCBKinStopPion[OPT_NWIND+1];
  Int_t nEneCBKinStopPion[OPT_NWIND+1];
  Double_t* EneCBKinLowPion[OPT_NWIND+1];
  Int_t nEneCBKinLowPion[OPT_NWIND+1];
  Double_t* EnePIDPion[OPT_NWIND+1];
  Int_t nEnePIDPion[OPT_NWIND+1];
  Double_t* EnePIDKinPion[OPT_NWIND+1];
  Int_t nEnePIDKinPion[OPT_NWIND+1];
  Double_t* EnePIDKinStopPion[OPT_NWIND+1];
  Int_t nEnePIDKinStopPion[OPT_NWIND+1];
  Double_t* EnePIDNormPion[OPT_NWIND+1];
  Int_t nEnePIDNormPion[OPT_NWIND+1];
  Double_t* MissingEnePion[OPT_NWIND+1];
  Int_t nMissingEnePion[OPT_NWIND+1]; 
  Double_t* MissingEnePIDPion[OPT_NWIND+1];
  Int_t nMissingEnePIDPion[OPT_NWIND+1];
  Double_t* MissingEnePIDProt[OPT_NWIND+1];
  Int_t nMissingEnePIDProt[OPT_NWIND+1];
  Double_t* EneCBStopPion[OPT_NWIND+1];
  Int_t nEneCBStopPion[OPT_NWIND+1];
  Double_t* EnePIDStopPion[OPT_NWIND+1];
  Int_t nEnePIDStopPion[OPT_NWIND+1];
  Double_t* CheckEneCBProt[OPT_NWIND+1];
  Int_t nCheckEneCBProt[OPT_NWIND+1];
  Double_t* CheckEneCBPion[OPT_NWIND+1];
  Int_t nCheckEneCBPion[OPT_NWIND+1];
  Double_t* CheckEnePIDPion[OPT_NWIND+1];
  Int_t nCheckEnePIDPion[OPT_NWIND+1];
  Double_t* CheckEnePIDProt[OPT_NWIND+1];
  Int_t nCheckEnePIDProt[OPT_NWIND+1];
  Double_t * TrfBadChi2[OPT_NWIND+1];
  Int_t nTrfBadChi2[OPT_NWIND+1];
  Double_t* TkineBadChi2[OPT_NWIND+1];
  Int_t nTkineBadChi2[OPT_NWIND+1];
  Double_t* TrfOkChi2[OPT_NWIND+1];
  Int_t nTrfOkChi2[OPT_NWIND+1];
  Double_t* TkineOkChi2[OPT_NWIND+1];
  Int_t nTkineOkChi2[OPT_NWIND+1];
  Double_t * TrfBadChi2Pion[OPT_NWIND+1];
  Int_t nTrfBadChi2Pion[OPT_NWIND+1];
  Double_t* TkineBadChi2Pion[OPT_NWIND+1];
  Int_t nTkineBadChi2Pion[OPT_NWIND+1];
  Double_t* TrfOkChi2Pion[OPT_NWIND+1];
  Int_t nTrfOkChi2Pion[OPT_NWIND+1];
  Double_t* TkineOkChi2Pion[OPT_NWIND+1];
  Int_t nTkineOkChi2Pion[OPT_NWIND+1];
  Double_t* CheckEneCBStopPion[OPT_NWIND+1];
  Int_t nCheckEneCBStopPion[OPT_NWIND+1];
  Double_t* OkEnePIDNormProt[OPT_NWIND+1];
  Int_t nOkEnePIDNormProt[OPT_NWIND+1];
  Double_t* OkEneCBProt[OPT_NWIND+1];
  Int_t nOkEneCBProt[OPT_NWIND+1];
  Double_t* OkEnePIDNormPion[OPT_NWIND+1];
  Int_t nOkEnePIDNormPion[OPT_NWIND+1];
  Double_t* OkEneCBPion[OPT_NWIND+1];
  Int_t nOkEneCBPion[OPT_NWIND+1];

 public:
  TA2ExampleRangeFit( const char*, TA2Analysis* );
  virtual ~TA2ExampleRangeFit();
  virtual void 		  LoadVariable();            // variables for display/cuts
  virtual void 		  PostInit();                // init after parameter input
  virtual void 		  Reconstruct();             // reconstruct detector info
  virtual TA2DataManager *CreateChild( const char*, Int_t ) { return NULL; }
  virtual void		  MarkEndBuffers();	     // Mark EndBuffer for the output arrays
  virtual void            ParseMisc(char* line);     // Parses additional information from cobnfiguration file
  virtual void	          ResetEvent();
  virtual void            Read_Cal();
  virtual void            ComputeKineticEne(Double_t Egamma, Double_t th_prot, 
					    Double_t &p_prot, Double_t &E_prot, Double_t &T_prot,
					    Double_t &th_pi, Double_t &p_pi, Double_t &E_pi, Double_t &T_pion, 
					    Int_t &ierr);
  virtual void            ComputeKineticEneNPi(Double_t Egamma, Double_t th_pi, 
					       Double_t &p_pi, Double_t &E_pi, Double_t &T_pi,
					       Double_t &th_n, Double_t &p_n, Double_t &E_n, Double_t &T_n, 
					       Int_t &ierr2);

  Double_t Read_CBProton(Double_t );
  Double_t Read_CBPion(Double_t);

  Double_t Read_PIDProton(Double_t);
  Double_t Read_PIDPion(Double_t);

  void Speed();                            // Do speed measurement (events/sec)

  ClassDef(TA2ExampleRangeFit,1)
};

#endif
