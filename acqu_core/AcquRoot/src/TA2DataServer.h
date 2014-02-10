//--Author	JRM Annand	12th Jan 2004   Convert Acqu++ -> AcquRoot
//--Rev 	JRM Annand
//--Rev 	JRM Annand      24th Mar 2004   Data stream merging
//--Rev 	JRM Annand      19th Apr 2004   Write data file
//--Rev 	JRM Annand       3rd Mar 2005   Batch-mode log-file directory
//--Rev 	JRM Annand       5th Nov 2006   Flexi start offset GetEvent()
//--Rev 	JRM Annand      13th Nov 2006   Delay in multi-source start
//--Rev 	JRM Annand      10th Dec 2006   FlushBuffers() limit on loop
//--Rev 	JRM Annand       8th Jan 2007   Attempt re-align during merge
//--Rev 	JRM Annand...   22nd Jan 2007   4v0 changes ARFile_t
//--Rev 	JRM Annand...   27th Apr 2007   ProcessHeader format iterface
//--Rev 	JRM Annand...   20th Nov 2007   Local DAQ coupling
//--Rev 	JRM Annand...   22nd Mar 2008   ProcessHeader error non-fatal
//--Rev 	JRM Annand...   21st Apr 2008   Call UpdateInfo (header merge)
//--Rev 	JRM Annand...    1st May 2008   Process() no end record fix
//--Rev 	JRM Annand...    3rd Jun 2008...const Char_t*...gcc 4.3
//--Rev 	JRM Annand...    1st Sep 2009   delete[]
//--Rev 	JRM Annand...   31st Aug 2012   Mk2 data buffer recognised
//                                              in data merging
//--Rev 	JRM Annand...   21st Nov 2012   More Mk2 data buffer recognised
//                                              in data merging
//--Rev 	JRM Annand...    1st Mar 2013   More Mk2 merge debugging
//--Rev 	JRM Annand...    6th Mar 2013   Add TA2TAPSMk2Format
//--Update	JRM Annand...   18th Sep 2013   Start sources inc usleep 100-200
//				      
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2DataServer
// Read data from file or network and make available to AcquRoot sorting
// process. Data may optionally be written to file
//
//---------------------------------------------------------------------------

#ifndef _TA2DataServer_h_
#define _TA2DataServer_h_

#include "TA2System.h"
#include "TAcquRoot.h"
#include "TA2DataSource.h"
#include "TA2DataFormat.h"
#include "TA2RingBuffer.h"
#include "TThread.h"                   // for threaded event loop
#include "SwapBytes.h"

// Subevent status markers
enum{ ESubEventOK, ESubEventLow, ESubEventHigh };

class ARFile_t;

class TA2DataServer : public TA2System {
protected:
  TAcquRoot* fAcquRoot;         // ptr to main Acqu - Root class
  ARFile_t* fDataOutFile;       // file to write output
  TA2DataSource** fDataSource;  // ptrs to data source classes
  TA2DataFormat** fDataFormat;  // ptrs to data format classes
  TA2RingBuffer** fSourceBuff;  // data buffers for each source
  TA2RingBuffer* fSortBuff;     // data passed to sort process
  TThread* fServerThread;       // To run data server as separate thread
  void* fHeaderBuff;            // Header information (e.g. in ACQU format)
  Char_t* fFileDir;             // directory info
  void** fEvStart;              // pointers to starts of subevents
  void** fEvTempStart;          // for sorting ID mismatch
  UInt_t* fEvLength;            // array of subevent lengths
  Int_t* fEvOffset;             // subevent # offsets
  Int_t* fNEvSearch;            // # subevents to search to find event-ID match
  Int_t* fEvMaxDiff;            // ignore subevent differences > fEvMaxDiff
  Int_t* fEvMax;                // ignore subevent diff if ev-ID > fEvMax
  Int_t* fEvStatus;             // status of subevents
  UInt_t* fOutBuff;             // output data buffer
  void* fEndBuff;               // end of buffer pointer
  Char_t* fEvTempBuff;          // temporary event store sorting mismatch
  Int_t fNDataSource;           // # of input data streams
  Int_t fNDataBuffer;           // # data buffers
  Int_t fNeventError;           // merge error counter
  Int_t fRecLen;                // record length of any output file
  Bool_t fIsStore;              // flag data storage
  Bool_t fIsSortLock;           // true if sort analyses all data
  Bool_t fIsHeaderInit;         // 1st header read?

public:
  TA2DataServer( const Char_t*, TAcquRoot* );
  // sets up inter-process communications
  virtual ~TA2DataServer();
  virtual void SetConfig(Char_t*, Int_t);
  virtual void Start();
  virtual void Run();
  virtual void StartSources();
  virtual void Process();
  virtual void MultiProcess();
  virtual Bool_t ProcessHeader();
  virtual void GetEvent(Int_t);
  virtual void CompareEvent();
  virtual void MergeBuffers(UInt_t = EDataBuff);
  virtual void FlushBuffers();

  // Inlined getters
  TA2DataSource** GetDataSource(){ return fDataSource; }
  TA2DataSource* GetDataSource(Int_t i){ return fDataSource[i]; }
  TA2DataFormat** GetDataFormat(){ return fDataFormat; }
  TA2DataFormat* GetDataFormat(Int_t i){ return fDataFormat[i]; }
  TA2RingBuffer** GetSourceBuff(){ return fSourceBuff; }
  TA2RingBuffer* GetSourceBuff(Int_t i){ return fSourceBuff[i]; }
  TA2RingBuffer* GetSortBuff(){ return fSortBuff; }
  void* GetHeaderBuff(){ return fHeaderBuff; }
  Char_t* GetFileDir(){ return fFileDir; }
  void** GetEvStart(){ return fEvStart; }
  void** GetEvTempStart(){ return fEvTempStart; }
  UInt_t* GetEvLength(){ return fEvLength; }
  Int_t* GetEvOffset(){ return fEvOffset; }
  Int_t* GetNEvSearch(){ return fNEvSearch; }
  Int_t* GetEvStatus(){ return fEvStatus; }
  UInt_t* GetOutBuff(){ return fOutBuff; }
  void* GetEndBuff(){ return fEndBuff; }
  Char_t* GetEvTempBuff(){ return fEvTempBuff; }
  Int_t GetNDataSource(){ return fNDataSource; }
  Int_t GetNDataBuffer(){ return fNDataBuffer; }
  Int_t GetNeventError(){ return fNeventError; }
  Int_t GetRecLen(){ return fRecLen; }
  void SetIsStore( Bool_t test ){ fIsStore = test; }
  Bool_t IsStore(){ return fIsStore; }
  Bool_t IsSortLock(){ return fIsSortLock; }
  Bool_t IsHeaderInit(){ return fIsHeaderInit; }
 

  ClassDef(TA2DataServer,1)     // for ROOT dictionary
};

//---------------------------------------------------------------------------
inline void TA2DataServer::GetEvent( Int_t n )
{
  // Mark start of event and determine length in bytes of event
  // Not counting the end-of-event marker

  TA2DataSource* ds = fDataSource[n];        // data source handler
  TA2RingBuffer* rb = fSourceBuff[n];        // data buffer handler
  TA2DataFormat* df = fDataFormat[n];        // data format handler
  Int_t length = 0;
  UInt_t* header;
  //  UInt_t* start;

 newbuffer:
  if( !fEvStart[n] ){
    rb->WaitFull();
    if( ds->IsSwap() )Swap4ByteBuff((UInt_t*)rb->GetStore(),ds->GetInRecLen());
    header = (UInt_t*)(rb->GetStore());
    //    if( *header != EDataBuff ) goto cleanup;
    if( *header != df->GetDataBuffHead() ) goto cleanup;
    fEvStart[n] = (Char_t*)header + df->GetStartOffset();
  }
  // Data format class determines length of event (bytes)
  length = df->GetEventLength(fEvStart[n]);
 cleanup:
  fEvLength[n] = length;
  if( length || !n ) return;
  rb->CNext();
  fEvStart[n] = NULL;
  goto newbuffer;
}

//---------------------------------------------------------------------------
inline void TA2DataServer::CompareEvent()
{
  // Check synchronisation via event tag or ID number
  // ie do the sub-events belonging to the same event.
  // If ID's don't match, optionally look for match in up to n (fNEvSearch[i])
  // subsequent events in data stream
  // Mark the status of each subevent....any bad ones will not be included
  // in the merged event.

  Int_t i,j,k;
  if(  fDataFormat[0]->GetEventIndex() != 0xffff ){
    UShort_t ev0 = fDataFormat[0]->GetEventID();
    for( i=1; i<fNDataSource; i++ ){
      k = fDataFormat[i]->GetEventID() + fEvOffset[i];
      fEvTempStart[i] = NULL;
      if( k == ev0 ) fEvStatus[i] = ESubEventOK;
      else{
	fprintf( fLogStream,
		 "Event ID Mismatch Error: Event# %d; ID0 %d ID%d %d\n",
		 fDataFormat[0]->GetEventNumber(),ev0,i,k );
	j = fNEvSearch[i];                 // max number of event scans
	fNeventError++;                    // increment error counter
	//	void* evst = fEvStart[i];          // save event start
	//	UInt_t evlen = fEvLength[i];       // save event length
	void* evst;
	UInt_t evlen;
	while( (k < ev0) && (j > 0) ){
	  //	  memcpy(fEvTempBuff, fEvStart[i], fEvLength[i]);
	  evst = fEvStart[i];
	  evlen = fEvLength[i];
	  fEvStart[i] = (Char_t*)fEvStart[i] + fEvLength[i] + 
	  fDataFormat[i]->GetStartOffset();
	  GetEvent(i);
	  k = fDataFormat[i]->GetEventID() + fEvOffset[i];
	  // Have you scanned too far in the event stream
	  if( k > ev0 ){
	    // Yes too far....rewind back one event and assume its OK
	    fEvTempStart[i] = fEvStart[i];
	    fEvStart[i] = evst;
	    fEvLength[i] = evlen;
	    k = ev0;
	    break;
	  }
	  j--;
	}
	// Mark the status of the subevent...it will only be included
	// in the main event if its OK
	if( (k == ev0) || (ev0 > fEvMax[i]) ) fEvStatus[i] = ESubEventOK;
	else if( k < ev0 ) fEvStatus[i] = ESubEventLow;
	else if( (k-ev0) < fEvMaxDiff[i] ) fEvStatus[i] = ESubEventHigh;
	else fEvStatus[i] = ESubEventOK;
      }
    }
  }
  else{
    UInt_t ev0 = *(UInt_t*)fEvStart[0];
    for( i=1; i<fNDataSource; i++ ){
      if( *(UInt_t*)fEvStart[i] != ev0 ){
	fNeventError++;
      }
    }
  }
  //  fEvStart[i] = (UInt_t*)fEvStart[i] + 1;
  //  fEvLength[i] -= sizeof(UInt_t);
}

//---------------------------------------------------------------------------
inline void TA2DataServer::FlushBuffers()
{
  // Data source 0, the controlling data source has detected an end-of-file
  // marker in the data...look for end markers in sources 1-n
  // limit flushes to the number of buffers in the ring buffer

  Int_t i;
  for( i=1; i<fNDataSource; i++ ){
    TA2RingBuffer* rb = fSourceBuff[i];        // data buffer handler
    UInt_t nb = rb->GetNbuff();                // total number of buffers
    UInt_t header, j;
    for(j=0; j<nb; j++){
      header = *(UInt_t*)(rb->GetStore());
      if( header == EEndBuff ){
	fprintf(fLogStream,"End Marker found for data source %d\n", i);
	rb->CNext();
	fEvStart[i] = NULL;                   // reset the event pointer
	return;
      }
      else rb->CNext();
    }
  }
}


#endif
