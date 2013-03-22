//--Author	JRM Annand   27th Jan 2004...Apapt A.Starostin code
//--Rev         JRM Annand   12th May 2004...TA2Calorimeter etc. added
//--Rev         JRM Annand   21st Nov 2004...TA2CosmicCal added
//--Rev         JRM Annand   15th Jul 2005...TA2CrystalBall, TA2TAPS
//--Update      JRM Annand   15th Jul 2005...Use TA2Analysis::LoadVariable
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2UserAnalysis
//
// User analysis recognises user-written apparatus and physics classes
// New apparati and physics should be entered in Map_t kKnownChild
// and also in the switch statement of CreateChild() where the apparatus
// constructers are called

#include "TA2UserAnalysis.h"
#include "TA2Calorimeter.h"
#include "TA2KensTagger.h"
#include "TA2Tagger.h"
#include "TA2CosmicCal.h"
#include "TA2CrystalBall.h"
#include "TA2CB.h"
#include "TA2TAPS.h"
#include "TA2Taps.h"
#include "TA2CentralApparatus.h"
#include "TA2GenericApp.h"
#include "TA2UserPhysics.h"
#include "TA2MesonPhysics.h"
#include "TA2BasePhysics.h"
#include "TA2TriggerPhysics.h"

// Recognised apparatus classes.
// The "standard" set is held in TA2Analysis
enum { EA2Calorimeter,
       EA2KensTagger, EA2Tagger, EA2CosmicCal,
       EA2CrystalBall, EA2CB, EA2TAPS, EA2Taps,
       EA2CentralApparatus, EA2GenericApp,
       EA2Physics, EA2UserPhysics, EA2MesonPhysics, EA2BasePhysics, EA2TriggerPhysics };

static const Map_t kKnownChild[] =
{
  //Apparati
  {"TA2Calorimeter",      EA2Calorimeter},
  {"TA2KensTagger",       EA2KensTagger},
  {"TA2Tagger",           EA2Tagger},
  {"TA2CosmicCal",        EA2CosmicCal},
  {"TA2CrystalBall",      EA2CrystalBall},
  {"TA2CB",               EA2CB},
  {"TA2TAPS",             EA2TAPS},
  {"TA2Taps",             EA2Taps},
  {"TA2CentralApparatus", EA2CentralApparatus},
  {"TA2GenericApp",       EA2GenericApp},
  //Physics
  {"TA2Physics",          EA2Physics},
  {"TA2UserPhysics",      EA2UserPhysics},
  {"TA2MesonPhysics",     EA2MesonPhysics},
  {"TA2BasePhysics",      EA2BasePhysics},
  {"TA2TriggerPhysics",   EA2TriggerPhysics},
  {NULL,                  -1}
};

//-----------------------------------------------------------------------------

TA2UserAnalysis::TA2UserAnalysis(const char* name)
                :TA2Analysis(name)
{
  kValidChild = kKnownChild; //Valid apparatus list
  PrivateInit();
}

//---------------------------------------------------------------------------

TA2UserAnalysis::~TA2UserAnalysis()
{
  //Free up allocated memory
}

//---------------------------------------------------------------------------

TA2DataManager* TA2UserAnalysis::CreateChild(const char* name, Int_t a)
{
  //Add creation of user-defined apparati here
  switch(a)
  {
   //Apparatus stuff:
   case EA2KensTagger:
    //Ken's tagger
    return new TA2KensTagger(name, this);
   case EA2Tagger:
    //Standard tagger
    return new TA2Tagger(name, this);
   case EA2Calorimeter:
    //General calorimeter... TAPS, CB etc.
    return new TA2Calorimeter(name, this);
   case EA2CosmicCal:
    //Cosmic calibration array of long bars
    return new TA2CosmicCal(name, this);
   case EA2CrystalBall:
    //Modded CB stuff
    return new TA2CrystalBall(name, this);
   case EA2CB:
    //Other modded CB stuff
    return new TA2CB(name, this);
   case EA2TAPS:
    //Modded TAPS stuff (Ralf)
    return new TA2TAPS(name, this);
   case EA2Taps:
    //Henning's TAPS class
    return new TA2Taps(name, this);
   case EA2CentralApparatus:
    //Alexander Mushkarenkov's MWPC/Tracking systems
    return new TA2CentralApparatus(name, this);
   case EA2GenericApp:
    //Generic
    return new TA2GenericApp(name, this);

   //Physics stuff:
   case EA2Physics:
    //Default (dummy physics)
    return new TA2Physics(name, this);
   case EA2UserPhysics:
    //User defined Physics
    return new TA2UserPhysics(name, this);
   case EA2MesonPhysics:
    //Demonstration physics class for pi0/eta -> 2gamma
    return new TA2MesonPhysics(name, this);
   case EA2BasePhysics:
    //Base physics class with particle collection and trigger
    return new TA2BasePhysics(name, this);
   case EA2TriggerPhysics:
    //Trigger threshold physics class
    return new TA2TriggerPhysics(name, this);

   default:
    PrintError((char*)name, "<Unknown apparatus..cannot continue>", EErrFatal);
  }
  return NULL;
}

//-----------------------------------------------------------------------------

