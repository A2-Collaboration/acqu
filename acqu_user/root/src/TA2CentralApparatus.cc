// SVN info: $Id: TA2CentralApparatus.cc 64 2011-09-21 22:28:32Z mushkar $
// AcquRoot
#include "TA2CentralApparatus.h"
#include "TA2PlasticPID.h"
#include "TA2CylMwpc.h"
#include "TA2CalArray.h"
// #include "HitCluster_t.h"
#include "TA2CentralTrack.h"
// #include "TA2MwpcTrack.h"
#include "TA2Particle.h"

// Default list of detector classes that the TA2CentralApparatus
// apparatus may contain
enum { ECB_PlasticPID, ECB_CalArray, ECB_CylMWPC, ECB_TA2CylMwpc};
static Map_t kValidDetectors[] =
{
  {"TA2PlasticPID", ECB_PlasticPID},
  {"TA2CylMwpc",    ECB_TA2CylMwpc},
  {"TA2CalArray",   ECB_CalArray},
  {NULL,            -1}
};

// Valid Keywords for command-line setup of CB apparatus
enum { ECBAngleLimits = 1000, ECBParticleCuts,
       ECBUseTracksBestMwpc, ECBDroop, ECBTrackLimits };
static const Map_t kCBKeys[] =
{
  {"AngleLimits:",       ECBAngleLimits},
  {"ParticleID-Cut:",    ECBParticleCuts},
  {"UseBestMwpcTracks:", ECBUseTracksBestMwpc},
  {"Droop:",             ECBDroop}, //To-do: It would be better to move it to TA2PlasticPID
  {"TrackLimits:",       ECBTrackLimits},
  {NULL,                 -1}
};

ClassImp(TA2CentralApparatus)
//_________________________________________________________________________________________
//  TA2CentralApparatus
//
//  The TA2CentralApparatus class describes the PID, cylindrical MWPC and Crystall Ball as a single apparatus and 
//  provides methods for geometrical assotiation of their signals combining them into the particles tracks.
//

//_______________________________________________________________________________
TA2CentralApparatus::TA2CentralApparatus( const char* name, TA2System* analysis  ) : TA2Apparatus( name, analysis, kValidDetectors )
{
  // Zero pointers and counters, add local SetConfig command list
  // to list of other apparatus commands

  // Test
  fWait = kFALSE;

  // Parameters/Options
  fUseTracksBestMwpc = kFALSE;
  // Angle limits
  fMaxPhiMwpcNaI  = 1000.;
  fMaxPhiTrackPid = 1000.;
  fMaxPhiPidNaI   = 1000.;
  fMaxPhiMwpcInterNaI = 1000.;
  fMaxPhiMwpcInterPid = 1000.;
  // Track limits
  fMaxTrack = 10;
  fLimitsPsVertex[0] =  1000.;
  fLimitsPsVertex[1] =  1000.;
  fLimitsPsVertex[2] = -1000.;
  // NaI 
  fFactorTrackLengthNaI = .9; // TEST option being tested
  // Pid droop
  fLengthPid = 30.;
  fNpointsDroopPid = 0;
  fLostCorrPosPid = NULL;
  fLostCorrFacPid = NULL;

  // Detectors
  // PID
  fPid = NULL;
  fNhitsPid = 0;
  fHitsPid  = NULL;
  fPositionsPid = NULL;
  fIhitsPidUsed = NULL;
  fIhitsPidNotUsed = NULL;
  // MWPC
  fMwpc = NULL;
  fNchambers = 0;
  fTracksMwpc = NULL;
  fTracksTrueMwpc = NULL;
  fTracksTrueCandidateMwpc = NULL;
  fTracksCandidateMwpc = NULL;
  fNintersTrueMwpc = NULL;
  fIintersTrueMwpc = NULL;
  fNintersCandidateMwpc = NULL;
  fIintersCandidateMwpc = NULL;
  // NaI
  fNaI = NULL;
  fNclNaI = 0;
  fClNaI  = NULL;
  fIdClNaI = NULL;
  fUsedClNaI = NULL;
  fPositionsNaI = NULL;

  // Tracks
  fNtracks = 0;
  fTracks  = NULL;
  fNchTracks = 0;
  fIchTracks  = NULL;
  fNneTracks = 0;
  fIneTracks = NULL;
  // Single MWPC tracks
  fNtracksSingleMwpc = 0;
  fTracksSingleMwpc = NULL;
  fNtracksSingleMwpcBest = 0;
  // Not used MWPC intersections
  fNintersTrueMwpcNotUsed = NULL;
  fNintersCandidateMwpcNotUsed = NULL;
  // For histograms
  fPhiMwpcTrackNaI = NULL;
  fPhiTrackPid = NULL;
  fPhiMwpcTrackPid = NULL;
  fPhiInterSingleMwpcNaI = NULL;
  fPhiTrackSingleMwpcPid = NULL;  
  fPhiMwpcInterNaI = NULL;
  fPhiMwpcInterPid = NULL;
  fPhiPidNaI = NULL;
  fTrackType = NULL;
  fTrackTheta = NULL;
  fTrackPhi = NULL;
  fEclNaI = NULL;
  fEhitPid = NULL;
  fEtrackMwpc = NULL;
  fPsVertex[0] = fPsVertex[1] = fPsVertex[2] = fPsVertexR = NULL;

  // Vertexes
  fNvertexes = 0;
  fVertexes = NULL;
  fIvertexes = NULL;
  // For histograms
  for(Int_t i=0; i<3; ++i) fVert[i] = NULL;
  fVertR = NULL;

  // ParticleInfo
  fDet     = NULL;
  fSize    = NULL;
  fCentral = NULL;
  fTime    = NULL;
  fType    = NULL;

  //
  AddCmdList(kCBKeys); // for SetConfig()
}

//_______________________________________________________________________________
TA2CentralApparatus::~TA2CentralApparatus()
{
  // Free up allocated memory
  DeleteArrays();
}

//_______________________________________________________________________________
void TA2CentralApparatus::DeleteArrays()
{
  // Free up all allocated memory

  // Parameters
  if (fLostCorrFacPid) delete [] fLostCorrFacPid;
  if (fLostCorrPosPid) delete [] fLostCorrPosPid;

  // Detectors
  if (fPid) {
    delete fPid; // PID
    delete [] fPositionsPid;
    delete [] fIhitsPidUsed;
    delete [] fIhitsPidNotUsed;
  }
  if (fMwpc) {
    delete fMwpc;// MWPC
    delete [] fIintersTrueMwpc;
    delete [] fIintersCandidateMwpc;
    delete [] fNintersTrueMwpcNotUsed;
    delete [] fNintersCandidateMwpcNotUsed;
  }
  if (fNaI) {
    delete fNaI; // NaI
    delete [] fPositionsNaI;
    delete [] fUsedClNaI;
  }
  
  // Tracks
  delete [] fTracks;
  delete [] fIchTracks;
  delete [] fIneTracks;
  //
  delete [] fTracksSingleMwpc;
  // For histograms
  delete [] fPhiMwpcTrackNaI;
  delete [] fPhiTrackPid;
  delete [] fPhiMwpcTrackPid;
  delete [] fPhiInterSingleMwpcNaI;
  delete [] fPhiTrackSingleMwpcPid;  
  delete [] fPhiMwpcInterNaI;
  delete [] fPhiMwpcInterPid;
  delete [] fPhiPidNaI;
  delete [] fTrackType;
  delete [] fTrackTheta;
  delete [] fTrackPhi;
  delete [] fEclNaI;
  delete [] fEhitPid;
  delete [] fEtrackMwpc;
  for (Int_t i=0;i<3;++i) delete [] fPsVertex[i];
  delete [] fPsVertexR;
  
  // Vertexes
  delete [] fVertexes;
  for (Int_t i=0; i<fMaxTrack; ++i) delete [] fIvertexes[i];
  delete [] fIvertexes;
  //
  for (Int_t i=0;i<3;++i) delete [] fVert[i];
  delete [] fVertR;
  delete [] fTrackDist;
  
  // Particles
  delete [] fParticleInfo;
  delete [] fType;
  delete [] fDet;
  delete [] fSize;
  delete [] fCentral;
  delete [] fTime;
}

//_______________________________________________________________________________
TA2DataManager*  TA2CentralApparatus::CreateChild(const char* name, Int_t dclass)
{
  // Create a TA2Detector class for use with the TA2CentralApparatus
  // Valid detector types are...
  // 1. TA2PlasticPID
  // 2. TA2CalArray
  // 3. TA2CylMwpc
  
  if ( !name ) name = Map2String( dclass );
  switch (dclass) {
    default:
      return NULL;
    case ECB_PlasticPID:
      fPid = new TA2PlasticPID( name, this );
      return fPid;
    case ECB_CalArray:
      fNaI = new TA2CalArray( name, this );
      return fNaI;
    case ECB_TA2CylMwpc:
      fMwpc =  new TA2CylMwpc( name, this );
      return fMwpc;
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::LoadVariable( )
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup.
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable
  
  TA2DataManager::LoadVariable("Ntracks",	&fNtracks,	EISingleX);
  TA2DataManager::LoadVariable("NchTracks",	&fNchTracks,	EISingleX);
  TA2DataManager::LoadVariable("NneTracks",	&fNneTracks,	EISingleX);
  
  TA2DataManager::LoadVariable("PhiMwpcTrackNaI",	fPhiMwpcTrackNaI,	EDMultiX);
  TA2DataManager::LoadVariable("PhiTrackPid",		fPhiTrackPid,		EDMultiX);
  TA2DataManager::LoadVariable("PhiMwpcTrackPid",	fPhiMwpcTrackPid,	EDMultiX);
  TA2DataManager::LoadVariable("PhiInterSingleMwpcNaI",	fPhiInterSingleMwpcNaI,	EDMultiX);
  TA2DataManager::LoadVariable("PhiTrackSingleMwpcPid",	fPhiTrackSingleMwpcPid,	EDMultiX);  
  TA2DataManager::LoadVariable("PhiMwpcInterNaI",	fPhiMwpcInterNaI,	EDMultiX);
  TA2DataManager::LoadVariable("PhiMwpcInterPid",	fPhiMwpcInterPid,	EDMultiX);
  TA2DataManager::LoadVariable("PhiPidNaI",		fPhiPidNaI,		EDMultiX);
  
  TA2DataManager::LoadVariable("TrackTypes",	 fTrackType,	EDMultiX);
  TA2DataManager::LoadVariable("TrackTheta",	 fTrackTheta,	EDMultiX);
  TA2DataManager::LoadVariable("TrackPhi",	 fTrackPhi,	EDMultiX);
  TA2DataManager::LoadVariable("EclNaI",	 fEclNaI,	EDMultiX);
  TA2DataManager::LoadVariable("EhitPid",	 fEhitPid,	EDMultiX);
  TA2DataManager::LoadVariable("EtrackMwpc",	 fEtrackMwpc,	EDMultiX);
  TA2DataManager::LoadVariable("PsVertexX",	 fPsVertex[0],	EDMultiX);
  TA2DataManager::LoadVariable("PsVertexY",	 fPsVertex[1],	EDMultiX);
  TA2DataManager::LoadVariable("PsVertexZ",	 fPsVertex[2],	EDMultiX);
  TA2DataManager::LoadVariable("PsVertexR",	 fPsVertexR,	EDMultiX);
  
  // Vertexes
  TA2DataManager::LoadVariable("VertexN",	&fNvertexes,	EISingleX);
  TA2DataManager::LoadVariable("VertexX",	 fVert[0],	EDMultiX);
  TA2DataManager::LoadVariable("VertexY",	 fVert[1],	EDMultiX);
  TA2DataManager::LoadVariable("VertexZ",	 fVert[2],	EDMultiX);
  TA2DataManager::LoadVariable("VertexR",	 fVertR,	EDMultiX);
  TA2DataManager::LoadVariable("TrackDist",	 fTrackDist,	EDMultiX);
  
  // ParticleInfo
  TA2DataManager::LoadVariable("Type",	 fType,	EIMultiX);
  
  TA2Apparatus::LoadVariable();
}

//_______________________________________________________________________________
void TA2CentralApparatus::SetConfig(Char_t* line, Int_t key)
{
  // Any special command-line input for Crystal Ball apparatus
  
  Double_t dparm[20];
  
  switch (key){
    case ECBUseTracksBestMwpc:
      // Use best MWPC tracks
      fUseTracksBestMwpc = kTRUE;
      cout << "TA2CentralApparatus: Best MWPC tracks are used" << endl;
      break;
    case ECBDroop:
      // Droop correction for light attenuation inside PID strips
      if( sscanf( line, "%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	          &fLengthPid, &fNpointsDroopPid,
	          dparm,    dparm+1,  dparm+2,  dparm+3,  dparm+4,  dparm+5,  dparm+6,  dparm+7,  dparm+8,  dparm+9,
		  dparm+10, dparm+11, dparm+12, dparm+13, dparm+14, dparm+15, dparm+16, dparm+17, dparm+18, dparm+19 ) < 2*fNpointsDroopPid + 2 )
      {
	PrintError(line,"<Droop correction for light attenuation inside PID strips>");
	return;
      }
      fLostCorrPosPid = new Double_t[fNpointsDroopPid];
      fLostCorrFacPid = new Double_t[fNpointsDroopPid];
      for(Int_t i=0; i<fNpointsDroopPid; ++i)
      {
	fLostCorrPosPid[i] = dparm[i];
	fLostCorrFacPid[i] = dparm[fNpointsDroopPid+i];
      }
      break;
    case ECBAngleLimits:
      // Angle difference limits
      if ( sscanf(line,"%lf%lf%lf%lf%lf",&fMaxPhiMwpcNaI,&fMaxPhiTrackPid,&fMaxPhiPidNaI,&fMaxPhiMwpcInterNaI,&fMaxPhiMwpcInterPid) != 5 )
      {
	PrintError(line,"<CB angle difference limits>");
	return;
      }
      fMaxPhiMwpcNaI  *= kDegToRad;
      fMaxPhiTrackPid *= kDegToRad;
      fMaxPhiPidNaI   *= kDegToRad;
      fMaxPhiMwpcInterNaI *= kDegToRad;
      fMaxPhiMwpcInterPid *= kDegToRad;
      break;
    case ECBTrackLimits:
      // Track limits
      if( sscanf(line,"%d%lf%lf%lf",&fMaxTrack,&fLimitsPsVertex[0],&fLimitsPsVertex[1],&fLimitsPsVertex[2]) != 4 )
      {
	PrintError(line,"<Track limits>");
	return;
      }
      break;
    default:
      // default main apparatus SetConfig()
      TA2Apparatus::SetConfig( line, key );
      break;
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::PostInit()
{
  // Initialise arrays used to correlate hits in NaI, PiD and MWPC detectors.
  // Load 2D cuts file and get the contained cuts classes
  // Demand particle ID class...if not there self destruct
  // Does not come back if called with EErrFatal
  
  // Pointer to PID hits, elements, paddle positions
  if (fPid)
  {
    fHitsPid = fPid->GetHits();
    fPositionsPid = new TVector2[fPid->GetNelement()];
    fIhitsPidUsed = new Int_t[fPid->GetNelement()];
    fIhitsPidNotUsed = new Bool_t[fPid->GetNelement()];
  }
  
  // Init the pointers to MWPC tracks
  if (fMwpc)
  {
    fNchambers = fMwpc->GetNchamber();
    //
    fRmwpc = fMwpc->GetR();
    fTracksMwpc = fMwpc->GetTracks();
    if (fUseTracksBestMwpc)
    {
      fTracksTrueMwpc          = fMwpc->GetTracksTrueBest();
      fTracksTrueCandidateMwpc = fMwpc->GetTracksTrueCandidateBest();
      fTracksCandidateMwpc     = fMwpc->GetTracksCandidateBest();
    }
    else
    {
      fTracksTrueMwpc          = fMwpc->GetTracksTrue();
      fTracksTrueCandidateMwpc = fMwpc->GetTracksTrueCandidate();
      fTracksCandidateMwpc     = fMwpc->GetTracksCandidate();
    }
    fNintersTrueMwpc = fMwpc->GetNintersTrue();
    fIintersTrueMwpc = new const Int_t*[fNchambers];
    fNintersCandidateMwpc = fMwpc->GetNintersCandidate();
    fIintersCandidateMwpc = new const Int_t*[fNchambers];
    fNintersTrueMwpcNotUsed = new Int_t[fNchambers];
    fNintersCandidateMwpcNotUsed = new Int_t[fNchambers];
    for (Int_t i=0; i<fNchambers; ++i)
    {
      fIintersTrueMwpc[i] = fMwpc->GetIintersTrue(i);
      fIintersCandidateMwpc[i] = fMwpc->GetIintersCandidate(i);
      fNintersTrueMwpcNotUsed[i] = 0;
      fNintersCandidateMwpcNotUsed[i] = 0;
    }
  }
  
  // Pointer to NaI clasters
  if (fNaI)
  {
    if (((TA2CalArray*)fNaI)->IsUCLAClustering())
      fClNaI = ((TA2CalArray*)fNaI)->GetClusterUCLA(); // UCLA claster analysis
    else
      fClNaI = fNaI->GetCluster(); // Standard claster analysis
    fIdClNaI = fNaI->GetClustHit();
    fPositionsNaI = new TVector3[fNaI->GetMaxCluster()];
    fUsedClNaI = new Bool_t[fNaI->GetMaxCluster()];
    for (UInt_t i=0; i<fNaI->GetMaxCluster(); ++i) fUsedClNaI[i] = kFALSE;
  }
  
  // Tracks
  fTracks = new TA2CentralTrack[fMaxTrack];
  fIchTracks = new Int_t[fMaxTrack];
  fIneTracks = new Int_t[fMaxTrack];
  fTracksSingleMwpc = new TA2CentralTrack[fMaxTrack];
  // For histograms
  fPhiMwpcTrackNaI = new Double_t[fMaxTrack+1];		fPhiMwpcTrackNaI[fMaxTrack] = EBufferEnd;
  fPhiTrackPid     = new Double_t[fMaxTrack+1];		fPhiTrackPid[fMaxTrack]     = EBufferEnd;
  fPhiMwpcTrackPid = new Double_t[fMaxTrack+1];		fPhiMwpcTrackPid[fMaxTrack] = EBufferEnd;
  fPhiInterSingleMwpcNaI = new Double_t[fMaxTrack+1];	fPhiInterSingleMwpcNaI[fMaxTrack] = EBufferEnd;
  fPhiTrackSingleMwpcPid = new Double_t[fMaxTrack+1];	fPhiTrackSingleMwpcPid[fMaxTrack] = EBufferEnd;
  fPhiMwpcInterNaI = new Double_t[fMaxTrack+1];		fPhiMwpcInterNaI[fMaxTrack] = EBufferEnd;
  fPhiMwpcInterPid = new Double_t[fMaxTrack+1];		fPhiMwpcInterPid[fMaxTrack] = EBufferEnd;
  fPhiPidNaI = new Double_t[fMaxTrack+1];		fPhiPidNaI[fMaxTrack] = EBufferEnd;
  fTrackType = new Double_t[fMaxTrack+1];
  fTrackTheta = new Double_t[fMaxTrack+1];
  fTrackPhi   = new Double_t[fMaxTrack+1];
  fEclNaI     = new Double_t[fMaxTrack+1];
  fEhitPid    = new Double_t[fMaxTrack+1];
  fEtrackMwpc = new Double_t[fMaxTrack+1];
  for (Int_t i=0;i<3;++i) fPsVertex[i] = new Double_t[fMaxTrack+1];
  fPsVertexR = new Double_t[fMaxTrack+1];
  
  // Vertexes
  Int_t fMaxVertex = fMaxTrack*(fMaxTrack-1)>>1;
  fVertexes = new TVector3[fMaxVertex];
  fIvertexes = new Int_t*[fMaxTrack];
  for (Int_t i=0; i<fMaxTrack; ++i)
  {
    fIvertexes[i] = new Int_t[fMaxTrack];
    fIvertexes[i][i] = -1;
  }
  // for histograming
  for (Int_t i=0;i<3;++i) fVert[i] = new Double_t[fMaxVertex+1];
  fVertR = new Double_t[fMaxVertex+1];
  fTrackDist = new Double_t[fMaxVertex+1];
  
  // ParticleInfo
  fMaxParticle = fMaxTrack;
  fParticleInfo = new TA2Particle[fMaxParticle];
  fType = new Int_t[fMaxParticle];
  fDet = new Int_t[fMaxParticle];
  fSize = new Int_t[fMaxParticle];
  fCentral = new Int_t[fMaxParticle];
  fTime = new Double_t[fMaxParticle];
  
  // Finally call the default apparatus post initialise
  TA2Apparatus::PostInit();
}

//_______________________________________________________________________________
void TA2CentralApparatus::Cleanup()
{
  // Reset event
  
  //
  TA2Apparatus::Cleanup();
  
  //
  for (Int_t i=0; i<fNparticle; ++i) fP4[i].SetXYZT(0.,0.,0.,0.);
  fNparticle = 0;
  
  // Reset Pid
  if (fPid) {
    fNhitsPidNotUsed = 0;
    for(UInt_t i=0; i<fPid->GetNelement(); ++i) {
      fIhitsPidUsed[i] = 0;
      fIhitsPidNotUsed[i] = kFALSE;
    }
  }
  
  // Reset Mwpc
  if (fMwpc) {
    for(Int_t i=0; i<fNchambers; ++i) {
      fNintersTrueMwpcNotUsed[i] = fNintersCandidateMwpcNotUsed[i] = 0;
    }
  }
  
  // Reset NaI
  if (fNaI) {
    for(UInt_t i=0; i<fNaI->GetMaxCluster(); ++i) {
      fUsedClNaI[i] = kFALSE;
    }
  }
  
  // Tracks
  fNtracks   = 0;
  fNchTracks = 0;
  fNneTracks = 0;
  // Tracks with 2 MWPC
  fMapMwpcTrueTrackNaI.clear();
  fMapMwpcTrueCandidateTrackNaI.clear();
  fMapMwpcCandidateTrackNaI.clear();
  // Tracks with 1 MWPC
  fNtracksSingleMwpc = 0;
  fMapMwpcTrueInterNaI.clear();
  fMapMwpcCandidateInterNaI.clear();
  fNtracksSingleMwpcBest = 0;
  // For histograms
  for (Int_t i=0; i<fMaxTrack; ++i)
  {
    fPhiMwpcTrackNaI[i] = kNullFloat;
    fPhiTrackPid[i]     = kNullFloat;
    fPhiMwpcTrackPid[i] = kNullFloat;
    fPhiInterSingleMwpcNaI[i] = kNullFloat;
    fPhiTrackSingleMwpcPid[i] = kNullFloat;
    fPhiMwpcInterNaI[i] = kNullFloat;
    fPhiMwpcInterPid[i] = kNullFloat;
    fPhiPidNaI[i] = kNullFloat;
    for (Int_t j=0; j<3; ++j) fPsVertex[j][i] = kNullFloat;
    fPsVertexR[i] = kNullFloat;
  }
  
  // Vertexes
  fNvertexes = 0;
}

//_______________________________________________________________________________
void TA2CentralApparatus::Reconstruct()
{
  // Main reconstruction function
  
  // # of hits in PID
  if (fPid)
  {
    Double_t rPid, phiPid;
    TVector3 **pos = fPid->GetPosition();
    fNhitsPid = fPid->GetNhits();
    for (Int_t i=0; i<fNhitsPid; ++i)
    {
      rPid   = pos[fHitsPid[i]]->X();		// Pid radius in mm
      phiPid = pos[fHitsPid[i]]->Z()*kDegToRad;	// Pid phi in radians
      fPositionsPid[i] = TVector2(rPid*TMath::Cos(phiPid),rPid*TMath::Sin(phiPid));
    }
  }
  else
    fNhitsPid = 0;
  
  // # of clastrers in NaI
  if (fNaI)
  {
    fNclNaI = fNaI->GetNCluster();
    for (Int_t i=0; i<fNclNaI; ++i)
    {
      fPositionsNaI[i] = *(fClNaI[fIdClNaI[i]]->GetMeanPosition())*fFactorTrackLengthNaI*10.;// TODO: init fFactorPenetroNaI in cfg
    }
  }
  else
    fNclNaI = 0;
  
  // PID-MWPC-CB geometrical association
  if (fMwpc) fMwpc->ResetWiresAndStripsUsage();
  MakeTracks();
  
  // Reconstruct vertexes
  MakeVertexes();
  
  // Mark end buffers
  MarkEndBuffers();

  // Test
//   Test();
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracks()
{
  // Make tracks of particles comming from the target
  
  if (fMwpc)
  {
    // Tracks with 2 MWPC
    if (fMwpc->GetNtracks() > 0)
    {
      // PID & MWPC(true) & NaI
      MakeTracksTrue(fTracksTrueMwpc,fMapMwpcTrueTrackNaI);
      // PID & MWPC(true-candidate) & NaI
      MakeTracksTrue(fTracksTrueCandidateMwpc,fMapMwpcTrueCandidateTrackNaI);
      // PID & MWPC(true)
      MakeTracksPidMwpc(fTracksTrueMwpc);
      // PID & MWPC(true-candidate)
      MakeTracksPidMwpc(fTracksTrueCandidateMwpc);
      // PID & MWPC(candidate) & NaI
      MakeTracksTrue(fTracksCandidateMwpc,fMapMwpcCandidateTrackNaI);
      // PID & MWPC(candidate)
      MakeTracksPidMwpc(fTracksCandidateMwpc);
      // MWPC(true)
      MakeTracksMwpc(fTracksTrueMwpc);
      // MWPC(true-candidate)
      MakeTracksMwpc(fTracksTrueCandidateMwpc);
    }
    
    // Tracks with 1 MWPC
    // PID & MWPC(true intersections) & NaI
    MakeTracksTrue(fIintersTrueMwpc,fNintersTrueMwpc,fMapMwpcTrueInterNaI);
    // PID & MWPC(true intersections)
    MakeTracksPidMwpc(fIintersTrueMwpc,fNintersTrueMwpc);
    // PID & MWPC(condidate intersections) & NaI
    MakeTracksTrue(fIintersCandidateMwpc,fNintersCandidateMwpc,fMapMwpcCandidateInterNaI);
    // PID & MWPC(condidate intersections)
    MakeTracksPidMwpc(fIintersCandidateMwpc,fNintersCandidateMwpc);
    
    // Tracks with 2 MWPC
    // MWPC(candidate tracks)
    if (fMwpc->GetNtracksCandidate() > 0)
    {
      MakeTracksMwpc(fTracksCandidateMwpc);
    }
    
    // Not used MWPC true intersections
    CountNotUsedIntersMwpc(fIintersTrueMwpc,fNintersTrueMwpc,fNintersTrueMwpcNotUsed);
    // Not used MWPC candidate intersections
    CountNotUsedIntersMwpc(fIintersCandidateMwpc,fNintersCandidateMwpc,fNintersCandidateMwpcNotUsed);
  }
  
  // Tracks without MWPC
  // PID & NaI
  MakeTracksPidNaI();
  // PID
  CountNotUsedPid();
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracksTrue(const map<Double_t,Int_t> *tracksMwpc, map<Double_t,pair<Int_t,Int_t> > &mapPhiMwpcNaI)
{
  // Find CBcluster & MWPCtrack & PIDpaddle associated events and
  // combine them into common tracks
  
  // All possible MWPC track & NaI combinations
  Double_t phi;
  Int_t iTrMwpc;
  for (map<Double_t,Int_t>::const_iterator iterTracksMwpc=tracksMwpc->begin(); iterTracksMwpc!=tracksMwpc->end(); ++iterTracksMwpc)
  {
    iTrMwpc = iterTracksMwpc->second;
    if (fMwpc->IsUsedWIE(iTrMwpc)) continue; // check if the strips or wires were already used
    //
    for (Int_t iClNaI=0; iClNaI<fNclNaI; ++iClNaI)
    {
      if (IsUsedClNaI(iClNaI)) continue; // check if the clNaI was used already
      phi = fTracksMwpc[iTrMwpc].Angle(fPositionsNaI[iClNaI]);
      if (phi > fMaxPhiMwpcNaI) continue;
      mapPhiMwpcNaI[phi] = make_pair(iTrMwpc,iClNaI);
    }
  }
  
  // The best MWPC track & NaI & PID and MWPC track & NaI pairs
  Int_t iClNaI, iPid;
  const TVector3 *r[2];
  TA2TrackLine track;
  for (map<Double_t,pair<Int_t,Int_t> >::iterator iterPhi = mapPhiMwpcNaI.begin(); iterPhi != mapPhiMwpcNaI.end(); ++iterPhi)
  {
    if (TooManyTracks()) break;
    //
    iTrMwpc = iterPhi->second.first;
    if (fMwpc->IsUsedWIE(iTrMwpc)) continue; // check if the strips or wires were already used
    iClNaI = iterPhi->second.second;
    if (IsUsedClNaI(iClNaI)) continue; // check if the clNaI was used already
    //
    r[0] = &(fTracksMwpc[iTrMwpc].GetOrigin());
    r[1] =  fPositionsNaI + iClNaI;
    if (track.BuildTrack(*r[0],*r[1],-1.,fLimitsPsVertex))
    {
      // Test the track & PID combinations and find the best one
      fPhiTrackPid[fNtracks] = fMaxPhiTrackPid;
      iPid = kNullHit;
      for (Int_t i=0; i<fNhitsPid; ++i)
      {
	phi = track.Angle(fPositionsPid[i]);
	if (TMath::Abs(phi) > TMath::Abs(fPhiTrackPid[fNtracks])) continue;
	fPhiTrackPid[fNtracks] = phi;
	iPid = i;
      }
      fPhiTrackPid[fNtracks] *= (iPid == kNullHit) ? kNullFloat/fMaxPhiTrackPid : kRadToDeg;
      fPhiMwpcTrackNaI[fNtracks] = iterPhi->first*kRadToDeg;
      // Add a new track
      AddTrack(iPid,
	       fTracksMwpc[iTrMwpc].GetIinter(0),
	       fTracksMwpc[iTrMwpc].GetIinter(1),
	       iClNaI,
	       track.GetOrigin(),
	       track.GetDirCos(),
	       track.GetPsVertex());
    }
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracksPidMwpc(const map<Double_t,Int_t> *tracksMwpc)
{
  // MWPC track & PID
  
  Int_t iTrMwpc, iPid;
  Double_t phi;
  for (map<Double_t,Int_t>::const_iterator iterTrackMwpc=tracksMwpc->begin(); iterTrackMwpc!=tracksMwpc->end(); ++iterTrackMwpc)
  {
    if (TooManyTracks()) break;
    //
    iTrMwpc = iterTrackMwpc->second;
    // Check if the strips or wires were already used
    if (fMwpc->IsUsedWIE(iTrMwpc)) continue;
    // Test the MWPC track & PID combinations and find the best one
    fPhiMwpcTrackPid[fNtracks] = fMaxPhiTrackPid;
    iPid = kNullHit;
    for (Int_t i=0; i<fNhitsPid; ++i)
    {
      phi = fTracksMwpc[iTrMwpc].Angle(fPositionsPid[i]);
      if (TMath::Abs(phi) > TMath::Abs(fPhiMwpcTrackPid[fNtracks])) continue;
      fPhiMwpcTrackPid[fNtracks] = phi;
      iPid = i;
    }
    if (iPid == kNullHit)
    {
      fPhiMwpcTrackPid[fNtracks] = kNullFloat;
      continue;
    }
    else
    {
      fPhiMwpcTrackPid[fNtracks] *= kRadToDeg;
      // Add a new track
      AddTrack(iPid,
	       fTracksMwpc[iTrMwpc].GetIinter(0),
	       fTracksMwpc[iTrMwpc].GetIinter(1),
	       kNullHit,
	       fTracksMwpc[iTrMwpc].GetOrigin(),
	       fTracksMwpc[iTrMwpc].GetDirCos(),
	       fTracksMwpc[iTrMwpc].GetPsVertex());
    }
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracksMwpc(const map<Double_t,Int_t> *tracksMwpc)
{
  // MWPC alone
  
  Int_t iTrMwpc;
  for (map<Double_t,Int_t>::const_iterator iterTracksMwpc=tracksMwpc->begin(); iterTracksMwpc!=tracksMwpc->end(); ++iterTracksMwpc)
  {
    if (TooManyTracks()) break;
    //
    iTrMwpc = iterTracksMwpc->second;
    if ( fMwpc->IsUsedWIE(iTrMwpc)) continue; // check if the strips or wires were already used
    // Add a new track
    AddTrack(kNullHit,
	     fTracksMwpc[iTrMwpc].GetIinter(0),
	     fTracksMwpc[iTrMwpc].GetIinter(1),
	     kNullHit,
	     fTracksMwpc[iTrMwpc].GetOrigin(),
	     fTracksMwpc[iTrMwpc].GetDirCos(),
	     fTracksMwpc[iTrMwpc].GetPsVertex());
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracksTrue(const Int_t **iInters, const Int_t *nInters, map<Double_t,Int_t> &mapMwpcInterNaI)
{
  // Find all possible NaI & MWPC1(2) intersection & (PID) combinations
  // and add the best ones to the fTracks array
  
  // Possible tracks
  const TA2MwpcIntersection *inter;
  const TVector3 *r[2];
  Double_t phi;
  TA2CentralTrack *track = fTracksSingleMwpc + fNtracksSingleMwpc;
  for (Int_t iCh=0; iCh<fNchambers; ++iCh)
  {
    for (Int_t iInter=0; iInter<nInters[iCh]; ++iInter)
    {
      if (TooManyTracksSingleMwpc()) break;
      //
      inter = fMwpc->GetInters(iCh) + iInters[iCh][iInter];
      if (fMwpc->IsUsedWIE(iCh,inter->GetIclI(),inter->GetIclW(),inter->GetIclE())) continue; // check if the strips or wires were already used
      r[0] = inter->GetPosition();
      for (Int_t iNaI=0; iNaI<fNclNaI; ++iNaI)
      {
	if (TooManyTracksSingleMwpc()) break;
	//
	if (IsUsedClNaI(iNaI)) continue; // Check if the clNaI was used already
	//
	r[1] = fPositionsNaI + iNaI;
	phi = r[0]->Angle(*r[1]);
	if (phi > fMaxPhiMwpcInterNaI) continue;
	if (track->BuildTrack(*r[0],*r[1],-1.,fLimitsPsVertex))
	{
	  track->SetIintersMwpc(kNullHit,kNullHit); // Reset IintersMwpc
	  track->SetIinterMwpc(iCh,iInters[iCh][iInter]);
	  track->SetIclNaI(iNaI);
	  mapMwpcInterNaI[phi] = fNtracksSingleMwpc;
	  fPhiInterSingleMwpcNaI[fNtracksSingleMwpc] = phi*kRadToDeg;
	  // Find the best PID association with the track
	  fPhiTrackSingleMwpcPid[fNtracksSingleMwpc] = fMaxPhiTrackPid;
	  track->SetIhitPid(kNullHit); // Reset IhitPid
	  for (Int_t iPid=0; iPid<fNhitsPid; ++iPid)
	  {
	    phi = track->Angle(fPositionsPid[iPid]);
	    if (TMath::Abs(phi) > TMath::Abs(fPhiTrackSingleMwpcPid[fNtracksSingleMwpc])) continue;
	    fPhiTrackSingleMwpcPid[fNtracksSingleMwpc] = phi;
	    track->SetIhitPid(iPid);
	  }
	  fPhiTrackSingleMwpcPid[fNtracksSingleMwpc] *= (track->GetIhitPid() == kNullHit) ? kNullFloat/fMaxPhiTrackPid : kRadToDeg;
	  ++fNtracksSingleMwpc;
	  ++track;
	}
      }
    }
  }
  
  // Best tracks
  Int_t iTrack, iInterMwpc, iCh, iNaI;
  for (map<Double_t,Int_t>::iterator iterTrack = mapMwpcInterNaI.begin(); iterTrack != mapMwpcInterNaI.end(); ++iterTrack) {
    if (TooManyTracks()) break;
    //
    iTrack = iterTrack->second;
    iCh = fTracksSingleMwpc[iTrack].HasSingleMwpc()-1;
    iInterMwpc = fTracksSingleMwpc[iTrack].GetIinterMwpc(iCh);
    if (fMwpc->IsUsedWIE(iCh,iInterMwpc)) continue; // Check if the strips or wires were already used
    iNaI = fTracksSingleMwpc[iTrack].GetIclNaI();
    if (IsUsedClNaI(iNaI)) continue; // Check if the clNaI was used already
    // Add track
    fPhiMwpcInterNaI[fNtracks] = fPhiInterSingleMwpcNaI[iTrack];
    fPhiTrackPid[fNtracks]     = fPhiTrackSingleMwpcPid[iTrack];
    AddTrack(fTracksSingleMwpc[iTrack]);
    ++fNtracksSingleMwpcBest;
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracksPidMwpc(const Int_t **iInters, const Int_t *nInters)
{
  // Single Mwpc & Pid
  
  const TA2MwpcIntersection *inter;
  Double_t phi;
  Int_t iPid;
  TA2CentralTrack *track = fTracksSingleMwpc + fNtracksSingleMwpc;
  for (Int_t iCh=0; iCh<fNchambers; ++iCh)
  {
    for (Int_t iInter=0; iInter<nInters[iCh]; ++iInter)
    {
      if (TooManyTracksSingleMwpc()) break;
      if (TooManyTracks()) break;
      //
      inter = fMwpc->GetInters(iCh) + iInters[iCh][iInter];
      if (fMwpc->IsUsedWIE(iCh,iInters[iCh][iInter])) continue; // check if the strips or wires were already used
      // Find iPid correlated with the intersection
      fPhiMwpcInterPid[fNtracks] = fMaxPhiMwpcInterPid;
      iPid = kNullHit;
      for (Int_t i=0; i<fNhitsPid; ++i)
      {
	phi = TVector2::Phi_mpi_pi( inter->GetPhi() - fPositionsPid[i].Phi() );
	if (TMath::Abs(phi) > TMath::Abs(fPhiMwpcInterPid[fNtracks])) continue;
	fPhiMwpcInterPid[fNtracks] = phi;
	iPid = i;
      }
      if (iPid == kNullHit)
      {
	// There is no Pid hit correlated with the intersection
	fPhiMwpcInterPid[fNtracks] = kNullFloat;
	continue;
      }
      else
      {
	// There is a Pid hit correlated with the intersection, biuld a new track
	track->BuildTrack(TVector3(0,0,0),*(inter->GetPosition()));
	track->SetIhitPid(iPid);
	track->SetIintersMwpc(kNullHit,kNullHit); // Reset IintersMwpc
	track->SetIinterMwpc(iCh,iInters[iCh][iInter]);
	track->SetIclNaI(kNullHit);
	fPhiMwpcInterPid[fNtracks] *= kRadToDeg;
	// Add a new track
	AddTrack(*track);
	++track;
	++fNtracksSingleMwpc;
	++fNtracksSingleMwpcBest;
      }
    }
  }
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeTracksPidNaI()
{
  // NaI & Pid or NaI alone (pure neutral claster)
  
  Double_t phi;
  Int_t iHitPid;
  for (Int_t iNaI=0; iNaI<fNclNaI; ++iNaI)
  {
    if (TooManyTracks()) break;
    if (IsUsedClNaI(iNaI)) continue; //Check if the clNaI was used already
    // Find iHitPid correlated with the NaI cluster
    fPhiPidNaI[fNtracks] = fMaxPhiPidNaI;
    iHitPid = kNullHit;
    for (Int_t iPid=0; iPid<fNhitsPid; ++iPid)
    {
      phi = TVector2::Phi_mpi_pi( fPositionsPid[iPid].Phi() - fPositionsNaI[iNaI].Phi() );
      if (TMath::Abs(phi) > TMath::Abs(fPhiPidNaI[fNtracks])) continue;
      fPhiPidNaI[fNtracks] = phi;
      iHitPid = iPid;
    }
    fPhiPidNaI[fNtracks] *= (iHitPid == kNullHit) ? kNullFloat/fMaxPhiPidNaI : kRadToDeg;
    AddTrack(iHitPid,
	     kNullHit,
	     kNullHit,
	     iNaI,
	     TVector3(0,0,0),
	     fPositionsNaI[iNaI],
	     TVector3(0,0,0));
  }
  
}

//_______________________________________________________________________________
void TA2CentralApparatus::CountNotUsedPid()
{
  // Count not used PID
  
  for (Int_t i=0; i<fNhitsPid; ++i)
  {
    if (IsUsedHitPid(i)) continue;
    fIhitsPidNotUsed[i]++;
    ++fNhitsPidNotUsed;
  }
  
}

//_______________________________________________________________________________
void TA2CentralApparatus::CountNotUsedIntersMwpc(const Int_t **iInters, const Int_t *nInters, Int_t *nIntersMwpcNotUsed)
{
  // Count not used Mwpc intersections
  
  const TA2MwpcIntersection *inter;
  for (Int_t iCh=0; iCh<fNchambers; ++iCh)
  {
    for (Int_t iInter=0; iInter<nInters[iCh]; ++iInter)
    {
      inter = fMwpc->GetInters(iCh) + iInters[iCh][iInter];
      if (fMwpc->IsUsedWIE(iCh,inter->GetIclI(),inter->GetIclW(),inter->GetIclE())) continue; // check if the strips or wires were already used
      nIntersMwpcNotUsed[iCh]++;
    }
  }
  
}

//_______________________________________________________________________________
void TA2CentralApparatus::MakeVertexes()
{
  // Reconstruct all possible vertexes (based on 2 tracks)
  
  const TVector3 *o1, *o2, *d1, *d2;
  TVector3 r1, r2;
  for (Int_t ii, i=0; i<fNchTracks-1; ++i)
  {
    ii = fIchTracks[i];
    if ( !fTracks[ii].HasMwpc() || !fTracks[ii].HasNaI() ) continue;
    o1 = &(fTracks[ii].GetOrigin());
    d1 = &(fTracks[ii].GetDirCos());
    for (Int_t jj, j=i+1; j<fNchTracks; ++j)
    {
      jj = fIchTracks[j];
      if ( !fTracks[jj].HasMwpc() || !fTracks[jj].HasNaI() ) continue;
      o2 = &(fTracks[jj].GetOrigin());
      d2 = &(fTracks[jj].GetDirCos());
      //
      fTrackDist[fNvertexes] = TA2Math::DistanceBetweenTwoLines(*o1,*d1,*o2,*d2,r1,r2);
      fVertexes[fNvertexes] = 0.5*(r1 + r2);
      //
      fVert[0][fNvertexes] = fVertexes[fNvertexes].X();
      fVert[1][fNvertexes] = fVertexes[fNvertexes].Y();
      fVert[2][fNvertexes] = fVertexes[fNvertexes].Z();
      fVertR[fNvertexes]   = fVertexes[fNvertexes].Perp();
      //
      fIvertexes[ii][jj] = fIvertexes[jj][ii] = fNvertexes++;
    }
  }
  
}

//_______________________________________________________________________________
void TA2CentralApparatus::AddTrack(const Int_t iHitPid, const Int_t iInterMwpc0, const Int_t iInterMwpc1, const Int_t iClNaI, const TVector3 &orig, const TVector3 &dircos, const TVector3 &psvertex)
{
  // Add a new track to fTracks array
  // NOTE: No check for value of fNtracks is done
  
  // Set track
  TA2CentralTrack &track = fTracks[fNtracks];
  track.Reset();
  track.SetTrack(iHitPid,iInterMwpc0,iInterMwpc1,iClNaI,orig,dircos,psvertex);
  // Set Pid info
  if (track.HasPid())
  {
    track.SetEhitPid(CalcEhitPid(iHitPid,orig,dircos));
    if (fPid->IsTime()) track.SetThitPid(fPid->GetTime(fHitsPid[iHitPid]));
  }
  // Set MWPC info
  if (track.HasMwpc())
  {
    track.SetEtrackMwpc(CalcEtrackMwpc(iInterMwpc0,iInterMwpc1));
    //   track.SetTtrackMwpc(); // TODO
  }
  // Set NaI info
  if (track.HasNaI())
  {
    track.SetEclNaI(fClNaI[fIdClNaI[iClNaI]]->GetEnergy());
    track.SetTclNaI(fClNaI[fIdClNaI[iClNaI]]->GetTime());
    track.SetMclNaI(fClNaI[fIdClNaI[iClNaI]]->GetNhits());
    track.SetCentralIndexNaI(fClNaI[fIdClNaI[iClNaI]]->GetIndex());
    track.SetSigmaTimeNaI(fNaI->GetSigmaTime());
    track.SetSigmaEnergyNaI(fNaI->GetSigmaEnergy(track.GetEclNaI()));
    track.SetSigmaPhiNaI(fNaI->GetSigmaPhi(fPositionsNaI[iClNaI].Theta()));
    track.SetSigmaThetaNaI(fNaI->GetSigmaTheta());
  }
  
  // Mark the PID hit as "used"
  if (iHitPid != kNullHit) SetUsedHitPid(iHitPid);
  // Mark the MWPC strips and wires clusters as "used"
  if (iInterMwpc0 != kNullHit) fMwpc->SetUsedWIE(0,iInterMwpc0);
  if (iInterMwpc1 != kNullHit) fMwpc->SetUsedWIE(1,iInterMwpc1);
  // Mark the NaI cluster as "used"
  if (iClNaI != kNullHit) SetUsedClNaI(iClNaI);
  
  // Neutral or charged
  if (track.IsNeutral())
    fIneTracks[fNneTracks++] = fNtracks;
  else
    fIchTracks[fNchTracks++] = fNtracks;
  
  // For histograms
  fTrackType[fNtracks]  = track.GetType();
  fTrackTheta[fNtracks] = track.GetTheta()*kRadToDeg;
  fTrackPhi[fNtracks]   = TVector2::Phi_0_2pi(track.GetPhi())*kRadToDeg;
  fEclNaI[fNtracks]     = track.GetEclNaI();
  fEhitPid[fNtracks]    = track.GetEhitPid();
  fEtrackMwpc[fNtracks] = track.GetEtrackMwpc();
  if ( track.HasMwpc() && track.HasNaI() )
  {
    // To plot PsVertex only for NaI & Mwpc tracks
    fPsVertex[0][fNtracks] = track.GetPsVertex().X();
    fPsVertex[1][fNtracks] = track.GetPsVertex().Y();
    fPsVertex[2][fNtracks] = track.GetPsVertex().Z();
    fPsVertexR[fNtracks]   = track.GetPsVertex().Perp();
  }
  
  // Add a new TA2Particle to the fParticleInfo array
  AddParticleInfo(track);
  
  // Increment the tracks counter
  ++fNtracks;
}

//_______________________________________________________________________________
void TA2CentralApparatus::AddTrack(const TA2CentralTrack &track)
{
  // Add the track to fTracks
  AddTrack(track.GetIhitPid(),
	   track.GetIinterMwpc(0),
	   track.GetIinterMwpc(1),
	   track.GetIclNaI(),
	   track.GetOrigin(),
	   track.GetDirCos(),
	   track.GetPsVertex());
}

//_______________________________________________________________________________
void TA2CentralApparatus::AddParticleInfo(const TA2CentralTrack &track)
{
  // Add a new particle to the fParticleInfo array
  
  fParticleInfo[fNparticle].Reset();
  
  fType[fNparticle] = kRootino;
  Double_t mass = fParticleID->GetMassMeV(kRootino);
  Int_t m, iPid;
  //Check dE/E cuts
  if ( fPCut && track.HasPid() && track.HasNaI() )
  {
    iPid = fHitsPid[track.GetIhitPid()];
    m = fPCutStart[iPid];					// Start cut index for pid[iPid]
    for (Int_t i=0; i<fNSectorCut[iPid]; ++i, ++m)		// Loop over specified cuts
    {
      if ( static_cast<TA2Cut2DM<Double_t>* >(fPCut[m])->TestXY(track.GetEclNaI(),track.GetEhitPid()) )	// Condition met?
      {
	fType[fNparticle] = GetCutPDGIndex(m);			//Set particle ID code...
	mass = fParticleID->GetMassMeV(fType[fNparticle]);	//...and particle mass
	break;							// cut OK so exit loop
      }
    }
  }
  else if ( track.IsNeutral() )
  {
    fType[fNparticle] = kGamma;
  }
  
  //
  fP4[fNparticle]      = TA2Math::LorentzVector(track.GetDirCos().Unit(), track.GetEclNaI(), mass);
  fTime[fNparticle]    = track.GetTclNaI();
  fSize[fNparticle]    = track.GetMclNaI();
  fCentral[fNparticle] = track.GetCentralIndexNaI();
  fDet[fNparticle] = EDetUndef;
  if (track.HasNaI())  fDet[fNparticle] += EDetNaI;
  if (track.HasMwpc()) fDet[fNparticle] += EDetMWPC;
  if (track.HasPid())  fDet[fNparticle] += EDetPID;
  //
  fParticleInfo[fNparticle].SetP4(fP4[fNparticle]);
  fParticleInfo[fNparticle].SetTime(fTime[fNparticle]);
  fParticleInfo[fNparticle].SetClusterSize(fSize[fNparticle]);
  fParticleInfo[fNparticle].SetCentralIndex(fCentral[fNparticle]);
  fParticleInfo[fNparticle].SetApparatus(EAppCB);
  fParticleInfo[fNparticle].SetDetector(fDet[fNparticle]);
  fParticleInfo[fNparticle].SetSigmaTime(track.GetSigmaTimeNaI());
  fParticleInfo[fNparticle].SetSigmaE(track.GetSigmaEnergyNaI());
  fParticleInfo[fNparticle].SetSigmaPhi(track.GetSigmaPhiNaI());
  fParticleInfo[fNparticle].SetSigmaTheta(track.GetSigmaThetaNaI());
  fParticleInfo[fNparticle].SetParticleIDA(fType[fNparticle]);
  fParticleInfo[fNparticle].SetVetoTime(track.GetThitPid());
  fParticleInfo[fNparticle].SetVetoEnergy(track.GetEhitPid());
  fParticleInfo[fNparticle].SetVetoIndex(track.GetIhitPid());
  fParticleInfo[fNparticle].SetDetectorA( fDet[fNparticle] - fParticleInfo[fNparticle].GetDetectors() );
  fParticleInfo[fNparticle].SetTrackIntersects(track.GetIinterMwpc(0),track.GetIinterMwpc(1));
  fParticleInfo[fNparticle].SetTrackEnergy(track.GetEtrackMwpc());
  fParticleInfo[fNparticle].SetTrackTime(track.GetTtrackMwpc());
  
  ++fNparticle;
  
}

//////////////////////////////////////////////////////////////////////
////////// TEST ////////////// TEST ///////////// TEST ///////////////
//////////////////////////////////////////////////////////////////////
//_______________________________________________________________________________
void TA2CentralApparatus::Test()
{
  // Test output
  
  // Print file name & event number
  cout << "____________________________________ TA2CentralApparatus _______________________________"<< endl;
  cout << "Event#: " << gAN->GetNEvent() << "\tAccepted#: " << gAN->GetNEvAnalysed() << "\tDAQ#: " << gAN->GetNDAQEvent() << endl;
  
  // NaI
  cout << "> NaI" << endl;
  cout << "Ncl: " << fNclNaI << endl;
  for (Int_t i=0; i<fNclNaI; ++i)
  {
    cout << i << "\tPhi: " << TVector2::Phi_0_2pi(fPositionsNaI[i].Phi())*kRadToDeg << "\tTheta: " << fPositionsNaI[i].Theta()*kRadToDeg << endl;
  }
  
  // Pid
  cout << "> Pid" << endl;
  cout << "Npid: " << fNhitsPid << "\tNpidNotUsed: " << fNhitsPidNotUsed << endl;
  cout << "NotUsedPidHits:\t";
  for (Int_t i=0; i<fNhitsPid; ++i)
  {
    cout << fIhitsPidNotUsed[i] << " ";
  }
  cout << endl;
  cout << "UsedPidHits:\t";
  for (Int_t i=0; i<fNhitsPid; ++i)
  {
    cout << fIhitsPidUsed[i] << " ";
  }
  cout << endl;
  for (Int_t i=0; i<fNhitsPid; ++i)
  {
    cout << i << "\tPhi: "<< TVector2::Phi_0_2pi(fPositionsPid[i].Phi())*kRadToDeg << endl;
  }
  
  // Posslible tracks with 2 MWPC
  // True
  cout << "> fMapMwpcTrueTrackNaI" << endl;
  for (map<Double_t,pair<Int_t,Int_t> >::const_iterator iterPhi = fMapMwpcTrueTrackNaI.begin(); iterPhi != fMapMwpcTrueTrackNaI.end(); ++iterPhi)
  {
    cout << "dPhi: " << iterPhi->first*kRadToDeg << "\tiTrMwpc: " << iterPhi->second.first << "\tiClNaI: " << iterPhi->second.second << endl;
  }
  
  // TrueCandidate
  cout << "> fMapMwpcTrueCandidateTrackNaI" << endl;
  for (map<Double_t,pair<Int_t,Int_t> >::const_iterator iterPhi = fMapMwpcTrueCandidateTrackNaI.begin(); iterPhi != fMapMwpcTrueCandidateTrackNaI.end(); ++iterPhi)
  {
    cout << "dPhi: " << iterPhi->first*kRadToDeg << "\tiTrMwpc: " << iterPhi->second.first << "\tiClNaI: " << iterPhi->second.second << endl;
  }
  
  // Candidate
  cout << "> fMapMwpcCandidateTrackNaI" << endl;
  for (map<Double_t,pair<Int_t,Int_t> >::const_iterator iterPhi = fMapMwpcCandidateTrackNaI.begin(); iterPhi != fMapMwpcCandidateTrackNaI.end(); ++iterPhi)
  {
    cout << "dPhi: " << iterPhi->first*kRadToDeg << "\tiTrMwpc: " << iterPhi->second.first << "\tiClNaI: " << iterPhi->second.second << endl;
  }
  
  // Tracks with 1 MWPC
  cout << "> Possible tracks with SingleMwpc" << endl;
  cout << "NtracksSingleMwpc: " << fNtracksSingleMwpc << endl;
  for (Int_t i=0; i<fNtracksSingleMwpc; ++i)
  {
    cout << i;
    cout << "\tiHitPid: " << fTracksSingleMwpc[i].GetIhitPid();
    cout << "\tiInterMwpc0: " << fTracksSingleMwpc[i].GetIinterMwpc(0) << "\tiInterMwpc1: " << fTracksSingleMwpc[i].GetIinterMwpc(1);
    cout << "\tiClNaI: " << fTracksSingleMwpc[i].GetIclNaI();
    cout << "\tType: " << fTracksSingleMwpc[i].GetType();
    cout << "\tPhi: " << TVector2::Phi_0_2pi(fTracksSingleMwpc[i].GetPhi())*kRadToDeg << "\tTheta: " << fTracksSingleMwpc[i].GetTheta()*kRadToDeg;
    cout << "\tPhiInterMwpcClNaI: " << fPhiInterSingleMwpcNaI[i];
    cout << "\tPhiTrSingleMwpcPid: " << fPhiTrackSingleMwpcPid[i] << endl;
  }
  cout << "NtracksTrueSingleMwpc: " << fMapMwpcTrueInterNaI.size() << endl;
  for (map<Double_t,Int_t>::const_iterator iter=fMapMwpcTrueInterNaI.begin(); iter!=fMapMwpcTrueInterNaI.end(); ++iter)
  {
    cout << "PhiInterSingleMwpcNaI: " << iter->first << "\tiTracksSingleMwpc: " << iter->second << endl;
  }
  cout << "NtracksCandidatesSingleMwpc: " << fMapMwpcCandidateInterNaI.size() << endl;
  for (map<Double_t,Int_t>::const_iterator iter=fMapMwpcCandidateInterNaI.begin(); iter!=fMapMwpcCandidateInterNaI.end(); ++iter)
  {
    cout << "PhiInterSingleMwpcNaI: " << iter->first << "\tiTracksSingleMwpc: " << iter->second << endl;
  }
  cout << "NtracksSingleMwpcBest: " << fNtracksSingleMwpcBest << endl;
  for (Int_t i=0; i<fNchambers; ++i)
  {
    cout << "iCh" << i << ":" << "\tNintersTrueMwpcNotUsed: " << fNintersTrueMwpcNotUsed[i] << "\tNintersCandidateMwpcNotUsed: " << fNintersCandidateMwpcNotUsed[i] << endl;
  }
  
  // Tracks
  cout << "===================== Tracks =======================" << endl;
  cout << "Ntracks: " << fNtracks << "\tNch: " << fNchTracks << "\tNne: " << fNneTracks << endl;
  for (Int_t i=0; i<fNtracks; ++i)
  {
    cout << "\tiHitPid: " << fTracks[i].GetIhitPid();
    cout << "\tiInterMwpc0: " << fTracks[i].GetIinterMwpc(0) << "\tiInterMwpc1: " << fTracks[i].GetIinterMwpc(1);
    cout << "\tiClNaI: " << fTracks[i].GetIclNaI();
    cout << "\tType: " << fTracks[i].GetType() << "\tPhi: " << TVector2::Phi_0_2pi(fTracks[i].GetPhi())*kRadToDeg << "\tTheta: " << fTracks[i].GetTheta()*kRadToDeg << "\tiCh: " << fTracks[i].HasSingleMwpc()-1 << endl;
  }
  
//   fWait = fNhitsPid ==2 && fNclNaI == 2 && fMwpc->GetNtracks() == 2;//fTracksTrueMwpc->size() > 0;
//   fWait = fMwpc->GetNtracks() == 2;
//      fWait = fNintersTrueMwpc[1] == 2 && fNtracks == 2;
//   fWait = gAN->GetNDAQEvent() == 3000;
//   fWait = fNtracksSingleMwpc > 0;

  // Wait?
  if (fWait)
  {
    if (cin.get()=='q') gROOT->ProcessLine(".q");
  }
  
  // Reset fWait
  fWait = kFALSE;
  
  cout << "________________________________________________________________________________________"<< endl;
}
