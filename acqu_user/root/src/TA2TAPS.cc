//--Author      R Gregor   3th Jul 2005
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2TAPS
//
// User coded version of TAPS apparatus for identification
// of charged particles.

#include "TA2TAPS.h"
#include "TA2KensTagger.h"
#include "TA2KensLadder.h"
#include "TA2Analysis.h"

// Default list of detector classes that the TA2TAPS
// apparatus may contain
enum { ETAPS_BaF2=1501, ETAPSVeto, ETAPSVetoInit, ETAPSPulse, ETAPSTof, ETAPSPsa,
       ETAPSFactor, ETAPSSigma, ETAPSVetodet,
       ETAPSPosShift, ETAPSTimeShift};

static Map_t kValidDetectors[] = {
  {"TA2TAPS_BaF2", ETAPS_BaF2},
  {NULL,                -1}
};

static const Map_t kTAPSKeys[] = {
  {"TAPSVeto:",      ETAPSVeto},
  {"TAPSVetoInit:",  ETAPSVetoInit},
  {"TAPS-Factor:",   ETAPSFactor},
  {"PSA-Sigma:",     ETAPSSigma},
  {"TAPSTOF:",       ETAPSTof},
  {"TAPSPSA:",       ETAPSPsa},
  {"Vetodetection:", ETAPSVetodet},
  {"Pulseshape:",    ETAPSPulse},
  {"PosShift:",      ETAPSPosShift},
  {"TimeShift:",     ETAPSTimeShift},
  {NULL,            -1}
};

ClassImp(TA2TAPS)

//-----------------------------------------------------------------------------
TA2TAPS::TA2TAPS(const char* name, TA2System* analysis)
        :TA2Apparatus(name, analysis, kValidDetectors)
{
  // Zero pointers and counters, add local SetConfig command list
  // to list of other apparatus commands
  fBaF2 = NULL;
  //fNWVetosAmount = 520;
  fTapsFudge = 1.0;
  particle = 0;
  crystal = 0;
  fMaxParticle = 0;
  fClusterTime = NULL;
  for(Int_t t=1; t<4; t++)
    PosShift[t] = 0.0;
  TimeShift = 0.0;
  AddCmdList( kTAPSKeys );                  // for SetConfig()
}


//-----------------------------------------------------------------------------
TA2TAPS::~TA2TAPS()
{
  // Free up allocated memory
  delete RealVetoHits;
  delete psaInfo;
  delete isVCharged;
  delete IsProton;
  if (fParticleID) delete fParticleID;

}

//-----------------------------------------------------------------------------
TA2DataManager*  TA2TAPS::CreateChild(const char* name, Int_t dclass)
{
  // Create a TA2Detector class for use with the TA2TAPS
  // Valid detector types are...
  // 1. TA2TAPS_BaF2

  if(!name) name = Map2String(dclass);
  switch (dclass)
  {
   case ETAPS_BaF2:
    fBaF2 = new TA2TAPS_BaF2(name, this);
    return fBaF2;
   default:
    return NULL;
  }
}

//---------------------------------------------------------------------------
void TA2TAPS::SetConfig(char* line, int key)
{
  // Load TAPS parameters from file or command line
  // TAPS specific configuration
  TCutG* CutPtr;

  switch( key )
  {
   case ETAPSFactor:
    //TAPS Fudge Factor
    if(sscanf(line, "%lf", &fTapsFudge) < 1 ) goto error;
    break;
   case ETAPSSigma:
    //TAPS Fudge Factor
    if(sscanf(line, "%lf%lf", &fTapsSigma1, &fTapsSigma2) < 1 ) goto error;
    break;
   case ETAPSVetodet:
    //TAPS Fudge Factor
    if(sscanf(line, "%lf", &fTapsVetodet) < 1 ) goto error;
    break;
   case ETAPSTof:
    //TAPS Fudge Factor
    if(sscanf( line, "%d", &fTapsTof) < 1) goto error;
        //fCutFile=new TFile("rew.root","READ");
    if(fTapsTof > 0)
    {
        fCutFile = new TFile("CutTOF.root","READ");
        CutPtr =(TCutG*)fCutFile->Get("CUTG");
        fProtonCut = new TCutG;
        *fProtonCut = *CutPtr;
        fCutFile->Close();
    }
    break;
  case ETAPSPsa:
    //TAPS Fudge Factor
    if(sscanf(line, "%d", &fTapsPsa) < 1) goto error;
    break;
  case ETAPSVetoInit:
    //Max number of Vetos
    if(sscanf(line, "%d", &fNWVetosAmount) < 1) goto error;
    VetoMap = new Int_t[fNWVetosAmount];
    for (UInt_t v=0; v<fNWVetosAmount; v++)
      VetoMap[v] = 0;
    break;
  case ETAPSVeto:
    //Misaligned or broken channels
    if(sscanf(line, "%d%d", &fNWVetos, &fVetoStatus) < 1) goto error;
    VetoMap[fNWVetos] = fVetoStatus;
    break;
  case ETAPSPulse:
    //Crystal settings for pulseshape
    if( sscanf( line, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d", &minRad[crystal], &angle1[crystal], &sigma1[crystal] , &angle2[crystal], &sigma2[crystal] , &angle3[crystal], &sigma3[crystal] , &angle4[crystal], &sigma4[crystal] , &angle5[crystal], &sigma5[crystal] , &angle6[crystal], &sigma6 [crystal], &angle7[crystal], &sigma7[crystal] , &angle8[crystal], &sigma8[crystal] , &angle9[crystal], &sigma9[crystal] , &angle10[crystal], &sigma10[crystal] , &angle11[crystal], &sigma11[crystal] ,&maxRad[crystal] ) != 24 ) goto error;
    crystal++;
    break;
  case ETAPSPosShift:
    if( sscanf(line, "%lf %lf %lf", &PosShift[1], &PosShift[2], &PosShift[3])!=3) goto error;
    printf("TAPS Shift for target position set to (%f, %f, %f)\n",
           PosShift[1], PosShift[2], PosShift[3]);
    break;
  case ETAPSTimeShift:
    if(sscanf(line, "%lf", &TimeShift)!=1) goto error;
    printf("TAPS Time Shift set to %f\n", TimeShift);
    break;
  default:
    // Command not found...possible pass to next config
    TA2Apparatus::SetConfig(line, key);
    break;;
  }
  return;
 error: PrintError(line);
  return;
}

//-----------------------------------------------------------------------------

void TA2TAPS::LoadVariable()
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

  //                            name        pointer          type-spec
  TA2DataManager::LoadVariable("NCharged",  &fNCharged,      EISingleX);
  TA2DataManager::LoadVariable("NProtons",  &fNtapsproton,      EISingleX);
  TA2DataManager::LoadVariable("NNeutrons",  &fNtapsneutron,      EISingleX);
  TA2DataManager::LoadVariable("NElectrons",  &fNtapselectron,      EISingleX);
  TA2DataManager::LoadVariable("NGammas",  &fNtapsgamma,      EISingleX);
  TA2DataManager::LoadVariable("NRootinos",  &fNtapsrootino,      EISingleX);
  TA2DataManager::LoadVariable("PSASG", fTestSG, EDMultiX);
  TA2DataManager::LoadVariable("PSALG", fTestLG, EDMultiX);
  TA2DataManager::LoadVariable("PSARad", psarad, EDMultiX);
  TA2DataManager::LoadVariable("PSAAng", psaang, EDMultiX);
  TA2DataManager::LoadVariable("STesting0", &STesting[0], EDSingleX);
  TA2DataManager::LoadVariable("STesting1", &STesting[1], EDSingleX);
  TA2DataManager::LoadVariable("STesting2", &STesting[2], EDSingleX);
  TA2DataManager::LoadVariable("STesting3", &STesting[3], EDSingleX);
  //
  TA2Apparatus::LoadVariable();
}


//-----------------------------------------------------------------------------

void TA2TAPS::PostInit( )
{
  // Initialise arrays used to correlate hits in BaF2 and Veto detectors.
  // Load 2D cuts file and get the contained cuts classes
  // Demand particle ID class...if not there self destruct
  // Does not come back if called with EErrFatal
  if( !fParticleID )
    PrintError("",
               "<Configuration aborted: ParticleID class MUST be specified>",
               EErrFatal);
  fMaxParticle = fBaF2->GetMaxCluster();
  RealVetoHits = new Int_t[520];
  IsProton = new Bool_t[510];
  minRad = new Int_t[510];
  maxRad = new Int_t[510];
  angle1 = new Double_t[510];
  angle2 = new Double_t[510];
  angle3 = new Double_t[510];
  angle4 = new Double_t[510];
  angle5 = new Double_t[510];
  angle6 = new Double_t[510];
  angle7 = new Double_t[510];
  angle8 = new Double_t[510];
  angle9 = new Double_t[510];
  angle10 = new Double_t[510];
  angle11 = new Double_t[510];
  sigma1 = new Double_t[510];
  sigma2 = new Double_t[510];
  sigma3 = new Double_t[510];
  sigma4 = new Double_t[510];
  sigma5 = new Double_t[510];
  sigma6 = new Double_t[510];
  sigma7 = new Double_t[510];
  sigma8 = new Double_t[510];
  sigma9 = new Double_t[510];
  sigma10 = new Double_t[510];
  sigma11 = new Double_t[510];
  psaInfo = new Int_t[fMaxParticle+1];
  isVCharged = new Bool_t[fMaxParticle+1];
  fClusterTime = new Double_t[fMaxParticle+1];
  fTestLG = new Double_t[fMaxParticle+1];
  fTestSG = new Double_t[fMaxParticle+1];
  psarad = new Double_t[fMaxParticle+1];
  psaang = new Double_t[fMaxParticle+1];
  fSize = new Int_t[fMaxParticle+1];
  fParticles = new TA2Particle[fMaxParticle+1];
  fVetoIndex = new UInt_t[fMaxParticle+1];

  // Finally call the default apparatus post initialise
  fTAGG = (TA2KensTagger*)((TA2Analysis*)fParent)->GetChild("TAGG");   //tagger
  fLadder = NULL;
  if(!fLadder) fLadder = (TA2KensLadder*)((TA2Analysis*)fParent)->GetGrandChild("Ladder");
  if(!fLadder) fLadder = (TA2KensLadder*)((TA2Analysis*)fParent)->GetGrandChild("FPD");
  if(!fLadder) printf("Error: Could not find any Ladder/FPD class\n");
  TA2Apparatus::PostInit();
}

//-----------------------------------------------------------------------------

void TA2TAPS::Reconstruct()
{
  UInt_t* id_clBaF2 = fBaF2->GetClustHit();                //central detector array
  HitCluster_t** clBaF2 = fBaF2->GetCluster();             //cluster array
  UInt_t nBaF2 = fBaF2->GetNCluster();                     //nr. of clusters
  fNparticle = nBaF2;                                      //particles=nr.ofclusters
  RealVetoHits[0] = EBufferEnd;                            //starting value->no Vetohits

  fTestLG[0] = EBufferEnd;
  fTestSG[0] = EBufferEnd;
  psarad[0] = EBufferEnd;
  psaang[0] = EBufferEnd;
  for(Int_t t=0; t<3; t++)
    STesting[t] = EBufferEnd;

  if (nBaF2==0) return;                                    // no Hits in TAPS -> exit
  UInt_t startID;
  TVector3 p_;

  fBitPattern = fBaF2->GetBitPattern();                    //tapspattern
  Int_t* vetopattern = fBitPattern->GetHits(0);            //vetopattern from tapspattern
  UInt_t* vetohits = fBitPattern->GetNHits();              //nr. of vetohits
  Int_t ntagg = fLadder->GetNhits(); //fTAGG->GetNparticle();                     //nr. of tagger hits
  Double_t* TaggTime = fLadder->GetTimeOR();               //taggertime
  Double_t* BaF2Time = fBaF2->GetTime();                   //baf2time
  fNCharged = 0;                                           //nr. of charged particles
  fNtapsproton = 0;
  fNtapsneutron = 0;
  fNtapselectron = 0;
  fNtapsgamma = 0;
  fNtapsrootino = 0;
  for (Int_t a0=0;a0<fMaxParticle;a0++)
  {
    isVCharged[a0]=false;  // charged signature from vetos
    psaInfo[a0]=0;  // charged signature from pulseshape
    IsProton[a0] =false;
  }

//Vetocorrection
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UInt_t VVV=0;
  UInt_t* vetoNeighbourArray;
  UInt_t vetoNeighbours;

  for(UInt_t VV = 0; VV < vetohits[0]; VV++)
  {
    if(vetopattern[VV] == 127 || vetopattern[VV] == 383) continue;         //vetos but missing baf2-> sort out
    if(VetoMap[vetopattern[VV]]==-1) continue;                             //broken veto, defined in dat file
    if(VetoMap[vetopattern[VV]]==0)                                        // normal vetos
    {
      if(vetopattern[VV] > 127)
      {
        if(vetopattern[VV] > 383)
          RealVetoHits[VVV] = vetopattern[VV]-2;          //move vetos down, because of missing baf2
        else
          RealVetoHits[VVV] = vetopattern[VV]-1;
      }
      else RealVetoHits[VVV] = vetopattern[VV];
    }
    else  // correct wrong cabled vetos, defined in dat file
    {
      if(vetopattern[VV] > 127)
      {
        if(vetopattern[VV] > 383)
          RealVetoHits[VVV] = VetoMap[vetopattern[VV]]-2;
        else
          RealVetoHits[VVV] = VetoMap[vetopattern[VV]]-1;
      }
      else RealVetoHits[VVV]=VetoMap[vetopattern[VV]];
    }
    VVV++;
  }
  RealVetoHits[VVV] = EBufferEnd;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for(UInt_t i=0; i<nBaF2; i++) // first all particles are rootinos (or gammas?)
  {
    p_ = *((clBaF2[id_clBaF2[i]])->GetMeanPosition());
    p_.SetX(p_.X() + PosShift[1]);
    p_.SetY(p_.Y() + PosShift[2]);
    p_.SetZ(p_.Z() + PosShift[3]);
    fParticleID->SetP4(&fP4[i], kRootino, (fTapsFudge*(clBaF2[id_clBaF2[i]])->GetEnergy()), &p_);
    fClusterTime[i] = (clBaF2[id_clBaF2[i]])->GetTime() + TimeShift;  // cluster time
    fSize[i] = (clBaF2[id_clBaF2[i]])->GetNhits();
    fVetoIndex[i] = EBufferEnd;
  }

//Particle ID
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //Time of Flight
  if(fTapsTof > 0)
  {
    for(UInt_t i=0; i<nBaF2; i++)                     //loop on all TAPS hits
    {
      IsProton[i] = false;
      for ( Int_t j = 0; j < ntagg; j++ )                    // loop on all Tagger hits
      {
        Double_t TimeDifference = TaggTime[j] -  BaF2Time[id_clBaF2[i]];
        if(fProtonCut->IsInside((fP4+i)->E(), TimeDifference))
        {
          IsProton[i] = true;                                //time of flight -proton
          break;
        }
      }
    }
  }
  //Veto
  for(UInt_t i=0; i<nBaF2; i++)                      //loop on all TAPS hits
  {
    vetoNeighbourArray = (clBaF2[id_clBaF2[i]])->GetNeighbour();
    vetoNeighbours = (clBaF2[id_clBaF2[i]])->GetNNeighbour();
    for (UInt_t i2=0; i2<fNWVetosAmount; i2++)
    {
      if ((UInt_t)RealVetoHits[i2] == id_clBaF2[i])
      {
        isVCharged[i] = true;                //central detector has Veto-Hit
        fVetoIndex[i] = RealVetoHits[i2];
        fNCharged++;
      }
      else if(fTapsVetodet==2)               //1st ring used when fTapsVetodet==2, defined in dat file
      {
        for(UInt_t v1=0; v1<vetoNeighbours; v1++)
        {
          if((UInt_t)RealVetoHits[i2]==vetoNeighbourArray[v1])
          {
            isVCharged[i] = true; //1st ring of central detector has veto-hit
            fVetoIndex[i] = RealVetoHits[i2];
            fNCharged++;
          }
        }
      }
      if (RealVetoHits[i2]==(Int_t)EBufferEnd) break;
    }
  }

  //Pulseshape
  for (UInt_t a1=0; a1 < nBaF2; a1++)
  {
    fTestLG[a1] = fBaF2->GetLGEnergy(id_clBaF2[a1]);
    fTestSG[a1] = fBaF2->GetSGEnergy(id_clBaF2[a1]);
    psaang[a1] = (atan2(fTestSG[a1],fTestLG[a1]) * 180 / TMath::Pi());
    psarad[a1] = (sqrt(fTestLG[a1]*fTestLG[a1]+fTestSG[a1]*fTestSG[a1]));
    if(psarad[a1] > minRad[id_clBaF2[a1]] && psarad[a1] < maxRad[id_clBaF2[a1]])
    {
      if(psarad[a1]<=20)
      {
        if(psaang[a1] > (angle1[id_clBaF2[a1]]-fTapsSigma1*sigma1[id_clBaF2[a1]])) psaInfo[a1]=1; //photonarea
        if(psaang[a1] < (angle1[id_clBaF2[a1]]-fTapsSigma2*sigma1[id_clBaF2[a1]])) psaInfo[a1]=2; //protonarea
      }
      else if(psarad[a1]<=30)
      {
        if(psaang[a1] > (angle2[id_clBaF2[a1]]-fTapsSigma1*sigma2[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle2[id_clBaF2[a1]]-fTapsSigma2*sigma2[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=40)
      {
        if(psaang[a1] > (angle3[id_clBaF2[a1]]-fTapsSigma1*sigma3[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle3[id_clBaF2[a1]]-fTapsSigma2*sigma3[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=60)
{
        if(psaang[a1] > (angle4[id_clBaF2[a1]]-fTapsSigma1*sigma4[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle4[id_clBaF2[a1]]-fTapsSigma2*sigma4[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=90)
      {
        if(psaang[a1] > (angle5[id_clBaF2[a1]]-fTapsSigma1*sigma5[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle5[id_clBaF2[a1]]-fTapsSigma2*sigma5[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=130)
      {
        if(psaang[a1] > (angle6[id_clBaF2[a1]]-fTapsSigma1*sigma6[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle6[id_clBaF2[a1]]-fTapsSigma2*sigma6[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=200)
      {
        if(psaang[a1] > (angle7[id_clBaF2[a1]]-fTapsSigma1*sigma7[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle7[id_clBaF2[a1]]-fTapsSigma2*sigma7[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=280)
      {
        if(psaang[a1] > (angle8[id_clBaF2[a1]]-fTapsSigma1*sigma8[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle8[id_clBaF2[a1]]-fTapsSigma2*sigma8[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=360)
      {
        if(psaang[a1] > (angle9[id_clBaF2[a1]]-fTapsSigma1*sigma9[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle9[id_clBaF2[a1]]-fTapsSigma2*sigma9[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=420)
      {
        if(psaang[a1] > (angle10[id_clBaF2[a1]]-fTapsSigma1*sigma10[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle10[id_clBaF2[a1]]-fTapsSigma2*sigma10[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
      else if(psarad[a1]<=535)
      {
        if(psaang[a1] > (angle11[id_clBaF2[a1]]-fTapsSigma1*sigma11[id_clBaF2[a1]])) psaInfo[a1]=1;
        if(psaang[a1] < (angle11[id_clBaF2[a1]]-fTapsSigma2*sigma11[id_clBaF2[a1]])) psaInfo[a1]=2;
      }
    }
    else psaInfo[a1] = 9; //out of range
    if(id_clBaF2[a1]==32)
    {
      STesting[0] = fTestLG[a1];
      STesting[1] = fTestSG[a1];
      STesting[2] = psaang[a1];
      STesting[3] = psarad[a1];
    }
  }
  fTestLG[nBaF2] = EBufferEnd;
  fTestSG[nBaF2] = EBufferEnd;
  psaang[nBaF2] = EBufferEnd;
  psarad[nBaF2] = EBufferEnd;

  for(UInt_t ap=0;ap<nBaF2;ap++)
  {
    startID=0xf0;  //binary for :       e(+/-) gamma proton neutron | pi(+/-) muon deuteron notavailable
                   //                      1      1     1      1          0     0       0         0
    if(fTapsVetodet > 0)
    {
      if (isVCharged[ap]==true) { startID = startID & 0xa0; }
                    //                     1      0     1      0    |     0     0       0         0
      else { startID = startID & 0x50; }
                    //                     0      1     0      1    |     0     0       0         0
    }
    //printf("StartID1=%i\n",startID);
    if(fTapsTof > 0)
    {
      if (IsProton[ap]==true) { startID = startID & 0x30; }
                    //                     0      0     1      1    |     0     0       0         0
      else { startID = startID & 0xc0; }
                    //                     1      1     0      0    |     0     0       0         0
    }
    //printf("StartID2=%i, psa used=, %i\n",startID,fTapsPsa);
    if(fTapsPsa > 0)
    {
      if(psaInfo[ap]==1) { startID = startID & 0xc0; }           //gammaarea
              //                     1      1     0      0    |     0     0       0         0
      else if (psaInfo[ap]==2){ startID = startID & 0x30; }       //protonarea
              //                     0      0     1      1    |     0     0       0         0
      else if (psaInfo[ap]==9){ startID = startID & 0xf1; }       //not available
                    //                     1      1     1      1    |     0     0       0         1
    }
    else startID = startID & 0xc0; //Assume particles in gammaarea, if PSA turned off
    //printf("StartID3=%i\n",startID);
    if (startID==0x20)         //proton
    {
      fParticleID->SetP4(fP4+ap, kProton,(fTapsFudge*(clBaF2[id_clBaF2[ap]])->GetEnergy()), (clBaF2[id_clBaF2[ap]])->GetMeanPosition());
      fPDG_ID[ap] = kProton;
      fNtapsproton++;
    }
    else if (startID==0x40)  //gamma
    {
      fParticleID->SetP4(fP4+ap, kGamma,(fTapsFudge*(clBaF2[id_clBaF2[ap]])->GetEnergy()), (clBaF2[id_clBaF2[ap]])->GetMeanPosition());
      fPDG_ID[ap] = kGamma;
      fNtapsgamma++;
    }
    else if (startID==0x80)  //e+/-
    {
      fParticleID->SetP4(fP4+ap, kElectron,(fTapsFudge*(clBaF2[id_clBaF2[ap]])->GetEnergy()), (clBaF2[id_clBaF2[ap]])->GetMeanPosition());
      fPDG_ID[ap] = kElectron;
      fNtapselectron++;
    }
    else if (startID==0x10)  //neutron
    {
      fParticleID->SetP4(fP4+ap, kNeutron,(fTapsFudge*(clBaF2[id_clBaF2[ap]])->GetEnergy()), (clBaF2[id_clBaF2[ap]])->GetMeanPosition());
      fPDG_ID[ap] = kNeutron;
      fNtapsneutron++;
    } // else if (startID== ) { //user dependend
    else fNtapsrootino++;
    //if no condition matches it is the particle from the start, so in original code a rootino

    SetParticleInfo(ap);
  }
}

//-----------------------------------------------------------------------------

inline void TA2TAPS::Cleanup()
{
  // Clear any arrays holding variables
 TA2DataManager::Cleanup();

 if(!fNparticle) return;

 for(Int_t i=0; i<fNparticle; i++)
 {
   fP4[i].SetXYZT(0.,0.,0.,0.);
   fPDG_ID[i] = kRootino;
   fClusterTime[i] = (Double_t)EBufferEnd;
   fSize[i] = 0;
   fVetoIndex[i] = EBufferEnd;
 }
 fP4tot.SetXYZT(0.,0.,0.,0.);
 fNparticle = 0;
}

//-----------------------------------------------------------------------------
