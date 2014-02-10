#ifndef __TA2HistManager_h__
#define __TA2HistManager_h__

//--Author	JRM Annand    9th Jan 2003
//--Rev		JRM Annand...26th Feb 2003...1st production version
//--Rev 	JRM Annand...28th Feb 2003...Separate from TA2System
//--Rev 	JRM Annand...15th Oct 2003...Upgrade Setup1D methods
//--Rev 	JRM Annand... 4th Feb 2004...General display parse, data-cuts
//--Rev 	JRM Annand...30th Mar 2004...2D names scaler plots
//--Rev 	JRM Annand...27th Jun 2004...Rate/multiscaler stuff
//--Rev 	JRM Annand... 5th Oct 2004...DecodeSaved()
//--Rev 	JRM Annand...13th Apr 2005...Cuts,Weights,3D hist,Tidy code
//--Rev 	JRM Annand...10th May 2005...'.' separator cut, weight names
//--Rev 	JRM Annand....6th Jul 2006...ZeroAll() also children
//--Update	JRM Annand....8th Nov 2008...Non-const channel widths
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2HistManager
//
// Contains methods for...
// Setup, filling and plotting of histograms
// All Acqu-Root analysis classes requiring histograms inherit from this one.
// It cannot itself be instantiated...pure virtual SetConfig().

#include "TList.h"
#include "TA2DataManager.h"
#include "TA2H1.h"                     // AcquRoot wrapper for TH1F
#include "TA2H2.h"                     // AcquRoot wrapper for TH2F
#include "TA2H3.h"                     // AcquRoot wrapper for TH3f

enum { EHist1D = 1, EHist2D = 2 };     // 1D, 2D histograms
enum { EHistSingleX=ESingleX, EHistMultiX=EMultiX, EHistSingleY=ESingleY, 
       EHistMultiY=EMultiY, EHistMultiXY=EMultiXY, EHistShapeX=EShapeX,
       EHistScalerX=EScalerX, EHistRateX=ERateX};

class TA2HistManager : public TA2DataManager {
private:
  const Map_t* fHistList;              // list of histogram names
protected:
  TList* f1Dhist;		       // list of 1D histograms
  TList* f2Dhist;		       // list of 2D histograms
  TList* f3Dhist;		       // list of 3D histograms
  UInt_t fN1Dhist;                     // no. 1D histograms
  UInt_t fN2Dhist;                     // no. 2D histograms
  UInt_t fN3Dhist;                     // no. 2D histograms
  Bool_t fIsDisplay;                   // histogramming on?
public:
  TA2HistManager( const Char_t*, const Map_t*, const Map_t* );
  virtual ~TA2HistManager();
  virtual void Decode( );
  virtual void DecodeSaved( );
  virtual void Reconstruct();
  virtual void SetConfig( Char_t*, int ) = 0;   // class-dependent setup
  virtual void ParseDisplay(const Char_t * );         // histogram setup
  virtual void ReadChanLim(Int_t, Double_t**, Char_t*);
  template<typename T>
    void Setup1D( Char_t*, T*, Int_t = EHistSingleX, Double_t* = NULL,
		 TA2Cut* = NULL, Stat_t* = NULL );
  template<typename T>
    void SetupMulti1D( Char_t*, T*, Int_t = EHistSingleX, Double_t* = NULL,
		      TA2Cut* = NULL, Stat_t* = NULL );
  template<typename Tx, typename Ty>
    void Setup2D( Char_t*, Tx*, Ty*, Int_t = EHistSingleX, Double_t* = NULL,
		 Double_t* = NULL, TA2Cut* = NULL, Stat_t* = NULL );
  template<typename Tx, typename Ty, typename Tz>
    void Setup3D( Char_t*, Tx*, Ty*, Tz*, Int_t = EHistSingleX,
		Double_t* = NULL, Double_t* = NULL, Double_t* = NULL,
		TA2Cut* = NULL, Stat_t* = NULL );
  void FillHist( );
  void FillAllHist();
  void Hist( Char_t*,Int_t=EHist1D,Axis_t=0,Axis_t=0,Axis_t=0,Axis_t=0 );
  void ZeroHist( Char_t* );
  void ZeroAll();
  void SaveHist();
  void ListHist();

  void H2D( Char_t* n,Int_t m=EHist2D,Axis_t xl=0,Axis_t xh=0,
	    Axis_t yl=0,Axis_t yh=0 ){
    Hist( n, m, xl, xh, yl, yh );
  }
  Bool_t IsDisplay(){ return fIsDisplay; }    // Is histograming turned on
  UInt_t GetN1Dhist(){ return fN1Dhist; }     // # 1D histograms
  UInt_t GetN2Dhist(){ return fN2Dhist; }     // # 2D histograms
  UInt_t GetN3Dhist(){ return fN3Dhist; }     // # 3D histograms
  TList* Get1Dhist(){ return f1Dhist; }       // list of 1D histograms
  TList* Get2Dhist(){ return f2Dhist; }       // list of 2D histograms
  TList* Get3Dhist(){ return f3Dhist; }       // list of 3D histograms
   // 
  ClassDef(TA2HistManager,1)   
};

//-----------------------------------------------------------------------------
inline void TA2HistManager::FillHist()
{
  // If display is enabled...
  // Step through lists of 1D and 2D histograms
  // Call the fill procedure for those found

  if( !fIsDisplay ) return;

  // 1D histogram list
  if( f1Dhist ){                               // if there is a 1D list
    TA2H1<void>* h1;                           // generic 1D histogram
    TIter nexthist( f1Dhist );                 // list iterator
    while( (h1 = (TA2H1<void>*)nexthist()) ){  // step thro' list
      ((TA2H*)h1)->Fill();
    }
  }
  // 2D histogram list
  if( f2Dhist ){                               // if there is a 2D list
    TA2H2<void,void>* h2;                      // generic 2D histogram
    TIter nexthist( f2Dhist );                 // list iterator
    while( (h2 = (TA2H2<void,void>*)nexthist()) ){
      ((TA2H*)h2)->Fill();
    }
  }
  // 3D histogram list
  if( f3Dhist ){                               // if there is a 3D list
    TA2H3<void,void,void>* h3;                 // generic 3D histogram
    TIter nexthist( f3Dhist );                 // list iterator
    while( (h3 = (TA2H3<void,void,void>*)nexthist()) ){ 
      ((TA2H*)h3)->Fill();
    }
  }
}

//-----------------------------------------------------------------------------
inline void TA2HistManager::Decode( )
{
  // Step through list of child analysis entities
  // and call the child decode method for each
  // If the child is flagged as a specific decoder
  // Evaluate and test data-cut conditions and conditionally fill histograms
  TIter nextchild(fChildren);
  TA2HistManager* child;
  fIsDecodeOK = ETrue;
  while( ( child = (TA2HistManager*)nextchild() ) ){
    child->Decode();
    if( child->IsDecode() ){ 
      if( !(child->GetDataCuts())->Test() ){
	fIsDecodeOK = EFalse;
	return;
      }
      child->FillHist();
    }
  }
}

//-----------------------------------------------------------------------------
inline void TA2HistManager::DecodeSaved( )
{
  // Step through list of child analysis entities
  // and call the child decode method for each
  // If the child is flagged as a specific decoder
  // Evaluate and test data-cut conditions and conditionally fill histograms
  TIter nextchild(fChildren);
  TA2HistManager* child;
  fIsDecodeOK = ETrue;
  while( ( child = (TA2HistManager*)nextchild() ) ){
    child->DecodeSaved();
    if( child->IsDecode() ){ 
      if( !(child->GetDataCuts())->Test() ){
	fIsDecodeOK = EFalse;
	return;
      }
      child->FillHist();
    }
  }
}

//-----------------------------------------------------------------------------
inline void TA2HistManager::Reconstruct( )
{
  // Step through list of child analysis entities
  // and call the child reconstruct method for each
  // If the child is flagged as a specific reconstructer
  // Evaluate and test data-cut conditions and conditionally fill histograms
  TIter nextchild(fChildren);
  TA2HistManager* child;
  fIsReconstructOK = ETrue;
  while( ( child = (TA2HistManager*)nextchild() ) ){
    child->Reconstruct();
    if( child->IsReconstruct() ){ 
      if( !(child->GetDataCuts())->Test() ){
	fIsReconstructOK = EFalse;
	return;
      }
      child->FillHist();
    }
  }
}
//-----------------------------------------------------------------------------
template<typename T>
void TA2HistManager::Setup1D( Char_t* line, T* xhb, Int_t mode, 
			      Double_t* axisbuff, TA2Cut* cut, Stat_t* wgt )
{
  // Read parameters for setting up a 1D histogram from file or command line
  // Expect parameters provided on input line:
  // name, number of channels, lower limit, upper limit
  // Create histogram and add to list of histograms to be filled

  Char_t hname[EMaxName];                  // histogram name
  Char_t* name;
  Int_t nparm;
  Int_t xChan;                          // x channel numbers
  Double_t xLow,xHigh;                  // x axis range
  const Char_t* title = this->ClassName();// name of class = hist title

  if( !f1Dhist ) f1Dhist = new TList(); // check list initialised
  if( !xhb ){                           // check a data buffer is provided
    PrintError(line,"<1D histogram: null ptr to variable>");
    return;
  }
  TObject* h1;

  // Name of histogram is preceded by the name of the calling class
  // Except in the case of TA2Analysis or derived classes
  if( this->InheritsFrom("TA2Analysis") ) name = hname;
  else{
    strcpy( hname, GetName() );
    strcat( hname, "_" );
    name = hname + strlen(hname);
  }
  nparm = sscanf( line, "%s%d%lf%lf", name,&xChan,&xLow,&xHigh );
  if( !strrchr( name,'_') ) name = hname;

  switch ( nparm ){
  case 2:
    if( !axisbuff ){
      PrintError(line,"<1D histogram: format error or null axis buffer>");
      return;
    }
  case 4:
    // 1D histograms
    switch( mode ){
    case EHistSingleX:
      // normal single value
      if( axisbuff )
      	h1 = new TA2H1S<T>(name,title,xChan,axisbuff,xhb,cut,wgt);
      else
      	h1 = new TA2H1S<T>(name,title,xChan,xLow,xHigh,xhb,cut,wgt);
      break;
    case EHistMultiX:
      // multi-value buffer
      if( axisbuff )
	h1 =new TA2H1M<T>(name,title,xChan,axisbuff,xhb,cut,wgt);
      else
	h1 =new TA2H1M<T>(name,title,xChan,xLow,xHigh,xhb,cut,wgt);
      break;
    case EHistShapeX:
      // flash ADC
      h1 = new TA2H1MF<T>(name,title,xChan,xhb);
      break;
    case EHistScalerX:
      // scaler spectrum
      h1 = new TA2H1MS<T>(name,title,xChan,xhb,
			  gAR->GetScalerRead());
      break;
    case EHistRateX:
      // rate or multiscaler (periodic fill of value)
      // here xlow is used to the determine the fill frequency
      // ie filled every xlow scaler-reads
      h1 = new TA2H1SR<T>(name,title,xChan,xhb,(Int_t)xLow,
			  gAR->GetScalerRead());
      break;
    default:
      PrintError(line,"<1D histogram: unknown spectrum type>");
      return;
    }
    f1Dhist->AddLast( h1 );
    fN1Dhist++;
    if( !fIsDisplay ) fIsDisplay = ETrue;
    return;
  default:
    PrintError(line,"<1D histogram: setup format error>");
    return;
  }
  return;
}

//-----------------------------------------------------------------------------
template<typename T>
void TA2HistManager::SetupMulti1D( Char_t* line, T* xhb, Int_t mode, 
				   Double_t* axisbuff,TA2Cut* cut,Stat_t* wgt )
{
  // Create multiple 1D histograms from command line
  // Expect parameters provided on input line:
  // name, lower histogram index, upper histogram index,
  // number of channels, lower histogram limit, upper histogram limit
  // Use Setup1D() to create individual histograms

  Char_t hname[EMaxName];               // name of histogram
  Int_t i,iLow,iHigh,xChan;           // multi-indices & x channel number
  Double_t xLow,xHigh;                // x axis range
  Char_t histline[256];                 // temporary line buffer

  UInt_t nparm = sscanf( line, "%s%d%d%d%lf%lf",
			hname,&iLow,&iHigh,&xChan,&xLow,&xHigh );
  if( nparm != 6 ){
    printf(" Error 1D histogram for class %s, at input line\n%s\n",
	   this->ClassName(), line);
    return;
  }
  for( i=iLow; i<=iHigh; i++ ){
    sprintf( histline,"%s%d %d %f %f", hname, i, xChan, xLow, xHigh );
    Setup1D( histline, xhb+i, mode, axisbuff, cut, wgt );
  }
  return;
}

//-----------------------------------------------------------------------------
template<typename Tx, typename Ty>
void TA2HistManager::Setup2D( Char_t* line, Tx* xhb, Ty* yhb, Int_t mode,
			      Double_t* xaxisbuff, Double_t* yaxisbuff,
			      TA2Cut* cut, Stat_t* wgt )
{
  // Read parameters for setting up a 2D histogram from file or command line
  // Create histogram and add to list

  Char_t hname[EMaxName];                 // histogram name
  Char_t* name;
  Int_t nparm;
  Int_t xChan, yChan;                	// x,y channel numbers
  Double_t xLow,xHigh,yLow,yHigh;   	// x,y axis range
  const Char_t* title = this->ClassName();// name of class = hist title

  if( !f2Dhist ) f2Dhist = new TList(); // check list initialised
  if( !xhb || !yhb ){            	// check a data buffer is provided
    PrintError(line,"<2D histogram: null ptr to x or y variable>");
    return;
  }
  TObject* h2;

  // Name of histogram is preceded by the name of the calling class
  // Except in the case of TA2Analysis or derived classes
  if( this->InheritsFrom("TA2Analysis") ) name = hname;
  else{
    strcpy( hname, GetName() );
    strcat( hname, "_" );
    name = hname + strlen(hname);
  }
  nparm = sscanf( line, "%s%d%lf%lf%d%lf%lf",
		  name,&xChan,&xLow,&xHigh,&yChan,&yLow,&yHigh );
  if( nparm != 7 ){
    PrintError(line,"<2D histogram: setup format error>");
    return;
  }
  switch( mode ){
  default:
    // default single-valued x & y
    if( xaxisbuff && yaxisbuff )
      h2 = new TA2H2SS<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else if( xaxisbuff )
      h2 = new TA2H2SS<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    else if( yaxisbuff )
      h2 = new TA2H2SS<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else
      h2 = new TA2H2SS<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    break;
  case EHistMultiX:
    // multi-valued x, single-valued y
    if( xaxisbuff && yaxisbuff )
      h2 = new TA2H2MS<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else if( xaxisbuff )
      h2 = new TA2H2MS<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    else if( yaxisbuff )
      h2 = new TA2H2MS<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else
      h2 = new TA2H2MS<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    break;
  case EHistMultiY:
    // single-valued x, multi-valued y
    if( xaxisbuff && yaxisbuff )
      h2 = new TA2H2SM<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else if( xaxisbuff )
      h2 = new TA2H2SM<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    else if( yaxisbuff )
      h2 = new TA2H2SM<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else
      h2 = new TA2H2SM<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    break;
  case EHistMultiXY:
    // multi-valued x and y
    if( xaxisbuff && yaxisbuff )
      h2 = new TA2H2MM<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else if( xaxisbuff )
      h2 = new TA2H2MM<Tx,Ty>(hname,title,xChan,xaxisbuff,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    else if( yaxisbuff )
      h2 = new TA2H2MM<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yaxisbuff,xhb,yhb,cut,wgt);
    else
      h2 = new TA2H2MM<Tx,Ty>(hname,title,xChan,xLow,xHigh,
			      yChan,yLow,yHigh,xhb,yhb,cut,wgt);
    break;
  }
  f2Dhist->AddLast( h2 );
  fN2Dhist++;
  if( !fIsDisplay ) fIsDisplay = ETrue;
  return;
}

//-----------------------------------------------------------------------------
template<typename Tx, typename Ty, typename Tz>
void TA2HistManager::Setup3D( Char_t* line, Tx* xhb, Ty* yhb, Tz* zhb,
			      Int_t mode, Double_t* xaxisbuff,
			      Double_t* yaxisbuff, Double_t* zaxisbuff,
			      TA2Cut* cut, Stat_t* wgt )
{
  // Read parameters for setting up a 3D histogram from file or command line
  // Create histogram and add to list

  Char_t hname[EMaxName];
  Char_t* name;
  TObject* h3;
  Int_t nparm;
  Int_t xChan, yChan, zChan;                	// x,y channel numbers
  Double_t xLow,xHigh,yLow,yHigh,zLow,zHigh;   	// x,y axis range
  const Char_t* title = this->ClassName();      // name of class = hist title

  if( !f3Dhist ) f3Dhist = new TList(); // check list initialised
  if( !xhb || !yhb || !zhb ){
    PrintError(line,"<3D Hist: null ptr to x,y or z variable>");
    return;
  }
  strcpy( hname, GetName() );
  strcat( hname, "_" );
  name = hname + strlen(hname);
  nparm = sscanf( line, "%s%d%lf%lf%d%lf%lf%d%lf%lf",
		  name,&xChan,&xLow,&xHigh,&yChan,&yLow,&yHigh,
		  &zChan, &zLow, &zHigh );
  if( nparm != 10 ){
    PrintError(line,"<3D histogram: bad setup line format>");
    return;
  }
  // If array to specify x-axis bin sizes provided...equivalent
  // y and z axis arrays must be provided
  if( xaxisbuff ){
    if( !yaxisbuff || !zaxisbuff ){
      PrintError(line,"<3D histogram:axis bin-size buffer for all or none>");
      return;
    }
  }
  switch( mode ){
  default:
    // default single-valued x,y,z
    if( xaxisbuff )
      h3 = new TA2H3S<Tx,Ty,Tz>(name,title,xChan,xaxisbuff,yChan,yaxisbuff,
				zChan,zaxisbuff,xhb,yhb,zhb,cut,wgt);
    else
      h3 = new TA2H3S<Tx,Ty,Tz>(name,title,xChan,xLow,xHigh,yChan,yLow,yHigh,
				zChan,zLow,zHigh,xhb,yhb,zhb,cut,wgt);
    break;
  case EHistMultiX:
    // multi-valued x,y,z
    if( xaxisbuff )
      h3 = new TA2H3M<Tx,Ty,Tz>(name,title,xChan,xaxisbuff,yChan,yaxisbuff,
				zChan,zaxisbuff,xhb,yhb,zhb,cut,wgt);
    else
      h3 = new TA2H3M<Tx,Ty,Tz>(name,title,xChan,xLow,xHigh,yChan,yLow,yHigh,
				zChan,zLow,zHigh,xhb,yhb,zhb,cut,wgt);
    break;
  }
  f3Dhist->AddLast( h3 );
  fN3Dhist++;
  if( !fIsDisplay ) fIsDisplay = ETrue;
  return;
  
}

#endif

