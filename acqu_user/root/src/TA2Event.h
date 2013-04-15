#ifndef __TA2Event_h__
#define __TA2Event_h__

#include "TA2Particle.h"

#define MAXBEAM 32
#define MAXPART 16

//-----------------------------------------------------------------------------

class TA2Event : public TObject
{
 protected:
  TA2Particle Beam[MAXBEAM];
  TA2Particle Particle[MAXPART];
  Int_t nBeam;
  Int_t nParticle;
  Int_t nEvent;
  UInt_t PatternL1;
  UInt_t PatternL2;
 public:
  TA2Event();
  virtual ~TA2Event();
  TA2Particle* GetBeam(){ return Beam; }
  TA2Particle GetBeam(Int_t pIndex){ return Beam[pIndex]; }
  TA2Particle* GetParticle(){ return Particle; }
  TA2Particle GetParticle(Int_t pIndex){ return Particle[pIndex]; }
  Int_t GetNBeam(){ return nBeam; }
  Int_t GetNParticle(){ return nParticle; }
  Int_t GetEventNumber(){ return nEvent; }
  UInt_t GetL1Pattern(){ return PatternL1; }
  UInt_t GetL2Pattern(){ return PatternL2; }
  void AddBeam(TA2Particle);
  void AddParticle(TA2Particle);
  void Clear();
  void SetEventNumber(Int_t pEvent){ nEvent = pEvent; }
  void SetEventNumber(){ nEvent++; }
  void SetL1Pattern(Int_t pPatternL1){ PatternL1 = pPatternL1; }
  void SetL2Pattern(Int_t pPatternL2){ PatternL2 = pPatternL2; }

  ClassDef(TA2Event, 1)
};

//-----------------------------------------------------------------------------

#endif
