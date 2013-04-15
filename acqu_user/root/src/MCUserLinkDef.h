//--Author	JRM Annand    9th Sep 2004
//--Rev
//--Update
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// ROOT dictionary specification
// All classes linked in libUserRoot.so must included in the ROOT
// dictionary so that they are recognised by the CINT command processor

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Monte Carlo classes
#pragma link C++ class TMCUserGenerator+;
#pragma link C++ class TMCUserParticle+;
#pragma link C++ class TMCHe4Phase+;

#endif
