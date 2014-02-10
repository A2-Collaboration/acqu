//--Author      JRM Annand    19th Jun 2003
//--Rev         JRM Annand...
//--Rev         ss,dk         When ??
//--Rev         JRM Annand     7th Dec 2003
//--Rev         D. Krambrich  17th Jan 2004
//--Rev         JRM Annand    30th Apr 2011
//--Rev         B. Oussena    10th Jun 2011  removed htons I2cWriteChk()
//                                           usleep() after a GeSiCA write
//--Rev         B. Oussena    13th Jun 2011  split (if Read(GeSiCA-Reg) & Value)
//--Rev         B. Oussena    27th Jun 2011  sleep(2) in TVME_GeSiCA::ProgSADC
//--Rev         B. Oussena    4th  Jul 2011  usleep(2) in TVME_GeSiCA::ProgSADC
//--Rev         B. Oussena    5th  Jul 2011  split (if Read(GeSiCA-Reg) & Value)
//                                           TVME_GeSiCA::SpyRead()
//--Rev         JRM Annand    9th  Jul 2011  Tidy up and further delay fiddle
//--Rev         JRM Annand   10th  Jul 2011  Use EIBase register Spy reset
//--Rev         JRM Annand   11th  Jul 2011  Try to improve speed
//--Rev         JRM Annand   24th  Jan 2012  Save TCS event ID
//--Rev         JRM Annand   25th  Jan 2012  Constructor set bits = 13
//--Rev         B. Oussena   22nd  Nov 2012  Add Send Event ID in SpyRead()
//--Rev         JRM Annand   17th  Sep 2013  Spy buff timeout 200 us
//--Update      JRM Annand   27th  Sep 2013  Try pause() spy buff wait
//
//
//--Description
//
//                *** Acqu++  ***
// Data Acquisition for Sub-Atomic Physics Experiments.
//
// GeSiCA
// Procedures for GeSiCA driver boards for SG-ADC
// Converted from ACQU 3v7 to AcquDAQ format JRMA 30th April 2011
// Initialisation Debugging by B. Oussena June-July 2011
//

/*
  Cleanups to do:
  - remove dummy-variable
  - i2c broadcast addr. was changed from 0 to 0xff -> check!
  - i2c addr. now starts at 0 to port 0 -> check if there is still a i+1 
  somewhere
  - readback of software set addr. might only look at 7 bits to compensate for 
  harware failure! chnge 0x7f to 0xff!
  - make setting thresholds only for init-level > 0 

  Changelog:
  19.01.04: -dk-
  - clear read-buffer data[]usleep(1);// added  -----------Baya
 in every event to have better debugging
  possibilities
  - check buffersize slink-header vs. hardware information
  - inplemented readout looking al #wrds in register +0x24. "Old" readout 
  got conserved. Please switch back as soon as possible-
  20.01.04: -dk-
  - little cleanup. Switched back to "old" readout procedure, but 
  implemnted errorblock on inconsistent eventsize.
  05.05.04: -dk-
  - added directive to log setting of thresholds to file.
  - minor cleanups. There is still a lot of work ro do!
  22.07.04: -dk-
  - added ivmeTVME_GeSiCA-i2cwait.cc~-output to "inconsistent buffersize"-errormessage.

  26.07.04: jrma, dk:
  - added "ESG_SINGLE-mode in ADC-sparse-mode.

*/

//------------------------------------------------------------------------------

#include "TVME_GeSiCA.h"
#include <arpa/inet.h>

//------------------------------------------------------------------------------
Map_t kGeSiCA_Keys[] = {
  {"NADCBoards:",    ENADCBoards},
  {"GeSiCAFile:",    EGeSiCAFile},
  {"GeSiCAReg:",     EGeSiCAReg},
  {"GeSiCAMode:",    EGeSiCAMode},
  {"SADCReg:",       ESADCReg},
  {"SADCFile:",      ESADCFile},
  {"SADCThreshold:", ESADCThresh},
  {"SADCSum:",       ESADCSum},
  {"GeSiCALogFile:",  EGeSiCALog},
};

//
// GeSiCA register address
//
VMEreg_t GeSiCAReg[] = {
  // offset,    data,  type,  repeat
  // Onboard GeSiCA registers
  {0x0,         0x1,    'l',    0},  // EIBase    0 Base
  {0x14,        0x0,    'l',    0},  // EIFPGAReg 1 FPGA download
  {0x20,        0x0,    'l',    0},  // EIMStatus 2 Module Status
  {0x24,	0x0,	'l',	0},  // EIDStatus 3 Data-Buffer Status
  {0x28,	0x0,	'l',	0},  // EIDatum   4 Data-Buffer
  {0x2c,        0x0,    'l',    0},  // EIHlPort  5 Hot-Link Port selection
  // I2C-registers
  {0x48,        0x0,	'l',	0},  // EI2CAddrCtl 6 Address and Control
  {0x4c,        0x0,	'l',	0},  // EI2CStatus  7 Status
  {0x40,        0x0,	'l',	0},  // EI2CdwLow   8 Data Write Low
  {0x54,        0x0,	'l',	0},  // EI2CdwHigh  9 Data Write High 3/11/03
  {0x44,        0x0,	'l',	0},  // EI2CdrLow  10 Data Read Low
  {0x58,        0x0,	'l',	0},  // EI2CdrHigh 12 Data Read High 3/11/03
  {0x50,        0x0,	'l',	0},  // EI2CPortSelect  13 Select i2c-port
  {0x48,        0x0,	'w',	0},  // genLD8 'w' = short
  {0xffffffff, 0, 0, 0}              // terminator
};


//------------------------------------------------------------------------------
TVME_GeSiCA::TVME_GeSiCA( Char_t* name, Char_t* file, FILE* log,
			  Char_t* line ):
  TVME_CATCH( name, file, log, line )
{
  // Basic initialisation 
  AddCmdList( kGeSiCA_Keys );                // CATCH-specific setup commands
  fNSample = fNSampleWd = fADCmode = fLatency = fOpMode = fNSADC = fNChip = 0;
  fSADCfile = NULL;
  fSADCport = fSADCchan = fSADCthresh = NULL;
  fNthresh=0;
  fNBits = 13;                             // 13-bit ADC
  fPedFactor = 0.0;
  for(Int_t i=0; i<3; i++){ fSamStart[i] = fSamWidth[i] = 0; }
}

//-----------------------------------------------------------------------------
TVME_GeSiCA::~TVME_GeSiCA( )
{
  // Clean up
}

//---------------------------------------------------------------------------
void TVME_GeSiCA::ReadIRQ( void** outBuffer )
{
  // Read GeSiCA (SG-ADC) Data Buffer via Spy Buffer (VMEbus)
  // This is a cut-down version just to try
  // Full functionality is in ReadIRQorig
  //
  UInt_t* pdatum;
  Int_t adcMode;                // 0: latch-all, 1: sparsified
  Int_t chipID;                 // at the moment 1: top / 0: bottom chip
  Int_t adcBlock;
  UShort_t index, value;

  UInt_t nword = SpyRead(outBuffer);
  if( !nword ) return;
  // Now decode the contents of the local buffer previously filled
  // from GeSiCA. Skip SLINK header...1st datum should be an ADC header
  // should be n_sgadc * 2 ADC-headers in the buffer
  pdatum = fSpyData + 2;
  adcMode = 0x1   & (*pdatum >> 24);        // 0: latch_all; 1: sparsified
  if( adcMode ){
    //Loop over all ADC Headers
    for(Int_t i = 0; i < fNChip; i++){
      chipID  = 0xf   & (*pdatum >> 26);        // ID of chip (0 or 1)
      adcBlock = (0xfff & (*pdatum >> 12)) - 1; // ADC data block size
      pdatum++;
      switch(fOpMode){
      case ESG_Triple: 
	// Save the 3 sums calculated by the SG-ADC
	// 3 values (sum 0..2) for every ADC channel over pedestal
	for(Int_t j=0; j<adcBlock; j+=3) {
	  index = 16 * chipID + (0xf & (*pdatum >> 26)) + fBaseIndex;//index
	  for(Int_t k=0; k<3; k++){
	    value = 0xffff & *pdatum++;
	    ADCStore(outBuffer,value,index);
	  }
	}
	break;
      case ESG_Sum1:
	// Save single sum calculated by the SADC
	for(Int_t j=0; j<adcBlock; j++){
	  index = 16 * chipID + (0xf & (*pdatum >> 26)) + fBaseIndex;//adc ind
	  value = 0xffff & *pdatum++;
	  ADCStore(outBuffer, value, index);   
	}
	break;
      default:
	// Anything else do nothing
	// It should not happen for the moment!
	return;
      }
    }
    return;
  }
  // Still to code in non-sparcified readout
  else return;
}

//---------------------------------------------------------------------------
void TVME_GeSiCA::ReadIRQorig( void** outBuffer )
{
  // Read GeSiCA (SG-ADC) Data Buffer via Spy Buffer (VMEbus)
  //
  volatile UInt_t* pdatum;
  Int_t adc_mode;                // 0: latch-all, 1: sparsified
  Int_t chip_id;                 // at the moment 1: top / 0: bottom chip
  Int_t sum;                     // total "integral" after ped. subtr.
  Int_t sum_1;                   // 1st "integral"
  Int_t sum_2;                   // 2nd "integral"
  UInt_t adc_blocksize;
  UShort_t index, value;
  UInt_t lvalue;
  //  unsigned Int_t nadc;
  UInt_t n;                  // current channel number
  Int_t k,l;

  UInt_t nword = SpyRead(outBuffer);
  if( !nword ) return;

  // Now decode the contents of the local buffer previously filled
  // from GeSiCA. Skip SLINK header...1st datum should be an ADC header
  // should be n_sgadc * 2 ADC-headers in the buffer
  pdatum = fSpyData + 2;
  //Loop over all ADC Headers
  for(Int_t i = 0; i < fNChip; i++){
    adc_mode = 0x1   & (*pdatum >> 24);      // 0: latch_all; 1: sparsified
    chip_id  = 0xf   & (*pdatum >> 26);      // id of chip (0 or 1 just now)
    adc_blocksize = 0xfff & (*pdatum >> 12); // ADC Block size
    pdatum++;
    
    // ADC in "Latch" mode
    if( adc_mode == 0 ){
      for(n = 0; n<EChanPerBlock; n++){
	index = (16 * chip_id) + n + fBaseIndex; // adc index
	switch(fOpMode){
	case ESG_All:
	  //save the samples + 3 sum outputs
	case ESG_Sample:
	  //save the flash-ADC samples
	case ESG_Sum1:
	case ESG_Sum2:
	  // Store the samples
	  for(k=0; k<fNSampleWd; k++){
	    lvalue = *pdatum++; //3 10-bit words packed into 32-bit datum
	    for( l=0; l<3; l++){
	      value = lvalue & 0x3ff; //10 bits of data
	      lvalue = lvalue >> 10;
	      ADCStore(outBuffer,value,index);
	    }
	  }
	  // Leave out the sums in this case
	  if (fOpMode == ESG_Sample){
	    pdatum +=3;
	    break;
	  }
	case ESG_Triple: //Save the 3 sums calculated by the SG-ADC
	  for(l=0; l<3; l++){ //Store the 3 sums after the sample
	    value = 0xffff & *pdatum++;
	    ADCStore(outBuffer,value,index);
	  }
	  break;
	case ESG_Single: //Save pedestal subtracted signal only
	  if(index!=(0xf & ( *pdatum >> 26)) + fBaseIndex)
	    fprintf(fLogStream, "<E> GESiCA: Decoder at wrong adc index\n");
	  sum_1 = (0xffff & *pdatum++); // get 1st sum
	  sum_2 = (0xffff & *pdatum++); // get 2nd sum
	  pdatum++; // skip 3rd sum
	  sum = sum_2 - Int_t(fPedFactor * sum_1);
	  if( sum <= 0)
	    value = -sum; //Signal at lower voltage than pedestal
	  else
	    value = 0;
	  ADCStore(outBuffer,value,index);
	  break;
	} //switch(sample)
      } //for(n = 0; n<EChanPerBlock; n++)
    }//if(adc_mode == 0)
    
    // SADC in sparce data mode....integrates samples
    // sample specifies the sparce-data mode
    else{
      switch(fOpMode){
      case ESG_All:
	// Save the sample + 3 sum outputs (error)
      case ESG_Sample:
	// Save the flash-ADC sample (error)
	fprintf(fLogStream, "<E> GeSiCA: acqu wants to store samples,\n");
	break; //Jump out here!
      case ESG_Triple: 
	// Save the 3 sums calculated by the SG-ADC
	// 3 values (sum 0..2) for every ADC channel over pedestal
	for(n=0; n<adc_blocksize-1; n+=3) {
	  index = 16 * chip_id + (0xf & (*pdatum >> 26)) + fBaseIndex;// adc index
	  for(l=0; l<3; l++){
	    value = 0xffff & *pdatum++;
	    ADCStore(outBuffer,value,index);
	  }
	}
	break;
      case ESG_Single:
	// Save sum #1 of the 3 (#0,1,2) calculated by SG-ADC
	for( n=0; n<adc_blocksize-1; n+=3) {
	  index = 16 * chip_id + (0xf & (*pdatum >> 26)) + fBaseIndex;//adc ind
	  sum_1 = 0xffff & *pdatum++;
	  sum_2 = 0xffff & *pdatum++;
	  ADCStore( outBuffer, sum_2 - sum_1, index ); // Store ped-subtr. sum
	  *pdatum++;                    // read to stay in phase with the adc.
	}
	break;
      case ESG_Sum2:
	// Save 2 sums calculated by the SADC (not yet tested)
	for(n=0; n<adc_blocksize-1; n+=2){
	  index = 16 * chip_id + (0xf & (*pdatum >> 26)) + fBaseIndex;// adc index
	  sum_1 = (0xffff & *pdatum++); // get 1st sum
	  sum_2 = (0xffff & *pdatum++); // get 1st sum
	  ADCStore(outBuffer,sum_1, index);
	  ADCStore(outBuffer,sum_2, index);
	}
	break;
      case ESG_Sum1:
	// Save single sum calculated by the SADC
	for(n=0; n<adc_blocksize-1; n++){
	  index = 16 * chip_id + (0xf & (*pdatum >> 26)) + fBaseIndex;//adc ind
	  value = 0xffff & *pdatum++;
	  ADCStore(outBuffer, value, index);   
	}
	break;
      }
    } // end of sparce mode
    
  }   // end if i_chips loop
  
  return;
}


//---------------------------------------------------------------------------
void TVME_GeSiCA::PostInit( )
{
  // Initialisation after all config data has been read in.
  //
  if( fIsInit ) return;
  // report on some critical stuff:
  //  printf(" Status of current (ID=0x%x) GeSiCA: 0x%x\n",
  //	 0xffffff & Read(EIBase), Read(EIMStatus));
  //  printf(" Data-Buffer-Status: 0x%x\n", Read(EIDStatus));
  InitReg( GeSiCAReg );
  TVMEmodule::PostInit();
  if( fInitLevel > EExpInit0 ){
    ProgFPGA();           // GeSiCA FPGA file
    ProgSADC();           // SADC FPGA files
    ProgOpMode();         // GeSiCA operational mode 
    ProgSampleSum();      // Sample integration boundaries 
    ProgThresh();         // SADC thresholds
  }
}

//---------------------------------------------------------------------------
void TVME_GeSiCA::SetConfig( Char_t* line, Int_t key )
{
  // Initialise GeSiCA and SG-ADC modules via VMEbus
  // Configuration from file 
  Char_t fpgafile[128];               // name of FPGA program file
  switch( key ){
  case ENADCBoards:
    // Get number SADCs connected to GeSiCA
    if(sscanf(line,"%d", &fNSADC ) != 1 )
      PrintError(line,"<Number ADC Parse Error>",EErrFatal);
    if ((fNSADC <= 0) || (fNSADC > 8))
      PrintError(line,"<Invalid Number ADC>",EErrFatal);
    fNChannel = 32*fNSADC;   // # ADC channels
    fNChip = 2*fNSADC;       // # ADC chips
    fSADCport = new Int_t[fNChannel];
    fSADCchan = new Int_t[fNChannel];
    fSADCthresh = new Int_t[fNChannel];
    break;
  case EGeSiCAFile:
    // GeSiCA FPGA file name
    if(sscanf(line,"%s", fpgafile ) != 1 )
      PrintError(line,"<SetConfig GeSiCA FPGA file parse error>",EErrFatal);
    fFPGAfile = BuildName(fpgafile);
    break;
  case EGeSiCAReg:
    // Write value to GeSiCA register at address offset
    // Not currently supported
    break;
  case EGeSiCAMode:
    // Set readout mode for GeSiCA
    if( sscanf(line,"%d%d%d%d",
	       &fADCmode, &fLatency, &fNSample, &fOpMode )!= 4 )
      PrintError(line,"<SetConfig GeSiCA Mode Parse Error>",EErrFatal);
    fNSampleWd = fNSample/3;// 3 SADC data in each 32-bit word from GeSiCA
    break;
  case ESADCReg: 
    // "Expert Tool": Directly write to SADCs registers
    // Not currently supported
    break;
  case ESADCFile:
    // Read the file for download to the SADCs
   if(sscanf(line,"%s", fpgafile ) != 1 )
      PrintError(line,"<SetConfig SADC FPGA file parse error>",EErrFatal);
    fSADCfile = BuildName(fpgafile);
    break;
  case ESADCSum:
    // Read the three sample regions for integrations:
    if( sscanf(line,"%d%d%d%d%d%d", fSamStart, fSamWidth,
	       fSamStart+1, fSamWidth+1, fSamStart+2, fSamWidth+2 )!= 6 )
      PrintError(line,"<Parse Error sample boundaries>",EErrFatal);
    fPedFactor = fSamWidth[1]/fSamWidth[0];
   break;
  case ESADCThresh:
    // Read SADC threshold
    if( sscanf(line,"%d%d%d", fSADCport+fNthresh, fSADCchan+fNthresh,
	       fSADCthresh+fNthresh) != 3 )
      PrintError(line,"<SetConfig Parse Threshold Error>",EErrFatal);
    fNthresh++;
    break;
  case EGeSiCALog:
    fprintf(fLogStream,"<Separate Log File not supported>\n");
    break;
  default:
    // default try commands of TDAQmodule
    TVME_CATCH::SetConfig(line, key);
    break;
  }
  return;
}

// ----------------------------------------------------------------------------
Char_t* TVME_GeSiCA::LoadRBT(Char_t *filename, Int_t* length)
{
  // Taken via cut and copy from old code.
  // Needs to be revised!
  // this routine reads from a datafile named "datafile",
  // extracts the 1/0s and writes them to an internal
  // buffer. The address of the buffer as well as the
  // length is returned.

  FILE *fp;
  Int_t len=0;
  Char_t line[128];
  Char_t next;
  Int_t nbits=0;
  Char_t *data;
  Int_t i;

  if ( (fp = fopen(filename,"r")) == NULL ) {
    PrintError(filename,"<File not Found>",EErrFatal);
  }
  // 7 is a magic number: the header length is 7 lines ...
  // The header contains the number  of bits to set. This
  // number will be used stopping the readout cycle due to
  // the fact that the file contains no EOF .. */
  for (i=0; i<7; i++) fgets(line, sizeof(line), fp);
  sscanf(line,"%*s %d",&len);
  data = new Char_t[len+1];

  while ( (next=fgetc(fp)) != -1 ) {       // read till EOF
    if ( (next == '0') || (next=='1') ){   // only binary bits
      data[nbits] = (Char_t)atoi(&next);   // char to integer
      if (nbits > len)
	PrintError(filename,"<Inconsistent Data... Check EOF>",EErrFatal);
      nbits++;
    }
  }
  if (nbits != len)
    fprintf(fLogStream,"Inconsistency in datasize!\nShould be: %i, is: %i\n",
	    len,nbits);
  // finished:
  *length = nbits;                         // return length of "rawdata"
  return data;                             // return the buffer
}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::ProgFPGA()
{
  // Send data to FPGA load register
  // Mind that lenght is not a pointer as at read_rbt!
#define PROGRAM 1
#define CCLK 2
#define DATA 4
#define INIT 1
#define DONE 2
  Int_t done=0;
  UShort_t datum=0;
  Char_t* data;
  Int_t length;

  // load conts of FPGA file into memory
  data = LoadRBT(fFPGAfile, &length);
  // Initialize chip
  Write(EIFPGAReg, 7);  // all bits high: PROGRAM=1
  usleep(1);
  Write(EIFPGAReg, CCLK);
  usleep(1);

  // Busywait for FPGA register ready
  // max. 2s. Assume usleep(1) = microsecond.
  for( Int_t i=0; i<=EI2CTimeout; i++ ){
    datum=Read(EIFPGAReg);
    usleep(1);
    if( datum & 0x1 ) break;   
    if( i >= EI2CTimeout )
      PrintError("","< Timeout GeSiCA FPGA register after 2 sec >", EErrFatal);
  }
  // Transmit data
  for (Int_t nb=0; nb<length; nb++) {
    // first step
    datum = data[nb]<<2;
    Write(EIFPGAReg, datum);
    // second step
    datum |= CCLK;
    Write(EIFPGAReg, datum);
    // Checking for done
    if ( (nb > (50*length/51)) ){
      datum = Read(EIFPGAReg);
      if( datum & DONE ) done++;
    }
  }
  if (done <= 2)
    PrintError("","<Insufficient DONE from Chip during FPGA load>",EErrFatal);
  else
    fprintf(fLogStream,"GeSiCA FPGA: %d bits downloaded from file %s\n",
	    length, fFPGAfile);
  delete data;
}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::ProgSADC(){
  // Download FPGA files to SADCs via optical link
  //
  Int_t i, j;
  Int_t dummy = 0 ;
  Int_t status = 0 ;
  UInt_t adc_id[8];
  sleep(1);                        // This sleep is very important
  Write(EIMStatus,0x4);            // enable readout via vme
  usleep(1);
  Int_t nadc;
  // check Module Status:
  nadc = 0;
  status = Read(EIMStatus);
  usleep(1); 
  printf(" SADC detected at Ports -- ");   
  for(i = 0x100, j=0; i<0x10000; i=i<<1, ++j){
    if (status & i) { printf(" %d", j); nadc++; }
  }
  printf("\n");
  // Add check that all attached SADCs are detected, FATAL error if not
  if( nadc != fNSADC )
    PrintError("","<SADC detect on GeSiCA port failed>",EErrFatal);
  if ((status & 0x7) != 0x7)
    PrintError("","<GeSiCA  not Synchronised, Status bits 7\n>",EErrFatal);
  // now: set the id's of the connected adcs and enable then for DAQ
  for(i = 0; i< fNSADC; i++){
    Write(EIHlPort,i);      // select Hot-Link-Port fPHlPort(0x2c)=1   
    //    usleep(2); 
    i2cReset();
    Write(EI2CPortSelect,0xff);      // target id: broadcast -in front of loop?
    //    usleep(2); 
    adc_id[i] = i2cReadB(0);         // read back hard-wired id 
    //    usleep(2); 
    adc_id[i] = adc_id[i] & 0xff ;   // read back hard-wired id 
    i2cWrite(1, adc_id[i] + (i<<8)); // set new id CHANGE! i+1 -> i
    //    usleep(2); 
    dummy = i2cReadB(1);  
    //    usleep(2); 
    dummy = dummy & 0xff;            // cleanup dummy 
    status = Read(EIMStatus);
    //    usleep(2);
    Write(EIMStatus,(status | (0x10000 << i)));
    //    usleep(2);
    printf(" Port %d: id 0x%02x changed to %d.\n ",i, adc_id[i], dummy);
    if ( dummy != i ){  // check...
      printf("ERROR: setting AdcID at port %d failed!\n", i);
    }
  }
  // Program (all) the FPGAs:
  // - Fast, but not without risk; 
  //   it may fail if one of the ADCs has timing problems.
  //   This one will the simply remain unprogrammed. 
  //   Might work around this
  //   by putting it to GeSiCA port 0, i.e. choosing it to be the master.
  // - Implement check as soon as hardware supports it easily
  fprintf(fLogStream, " Programming SADC fpgas (%s)\n", fSADCfile);
  Write(EI2CPortSelect,0xff);    // set broadcast mode
  //  usleep(2);
  Write(EIHlPort,(UInt_t)0x00);  // talk first module, others listening
  //  usleep(2);
  ProgramSgFPGA();
  // make sure i2c is disabled everywhere..
  i2cDisableAll();
  printf("*** SADC FPGA download from GeSiCA %s complete ***\n",
	 this->GetName());
}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::ProgramSgFPGA(){
  UInt_t datum=0;
  Char_t* data;
  Int_t length;
  UShort_t i2c_data = 0;
  UShort_t tmp_i2c_data = 0;

  // Read file contents into memory
  data = LoadRBT(fSADCfile,&length);
  // reset sg_adc fpga
  i2cWriteB(2,0);
  // set program bit
  i2cWriteB(2,4);
  datum = 0;
  while( datum == 0){
    Write(EI2CAddrCtl, 0x294);
    //    usleep(2);
    // wait for end of i2c command
    for( Int_t i=0; i<=EI2CTimeout; i++ ){
      datum = Read(EI2CdwLow);
      usleep(1);  
      if( !(datum & 0x1) )break;
      if( i == EI2CTimeout )
	PrintError("","< Timeout SADC FPGA register after 2 sec >",EErrFatal);
    }
    // check i2c ACK
    for( Int_t i=0; i<3; i++ ){ 
      datum = Read(EI2CStatus); // ----------------------------------- Baya
      usleep(2);                
      if (datum & 0x08);  // ----------------------------------- Baya
      else break;
      if( i == 2 )PrintError("","<SADC I2C timeout on ACK received>",EErrFatal);
    }
    datum = Read(EI2CdrLow); // Init bit xilinx chip   // ---------- Baya
    //    usleep(2);   
    datum = (datum & 0x1);  // ----------------------------------- Baya
  }
  fprintf(fLogStream, " ");
  for (Int_t nb=0; nb<length; nb+=16) {
    i2c_data = 0;
    tmp_i2c_data = 0;
    for(Int_t nbj = 0; nbj<16; nbj++){
      tmp_i2c_data = ((tmp_i2c_data << 1) & 0xfffe) + data[nb+nbj] ;
    }
    i2c_data = ((tmp_i2c_data & 0x00ff) << 8) + ((tmp_i2c_data & 0xff00) >> 8);
    i2cWrite( 3  , i2c_data);
  
    if ((nb & 0xffff) == 0) fprintf(fLogStream, ".");
  }
  // Finished
  fprintf(fLogStream, " SADC download done\n");
  delete data;
}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::ProgOpMode(){
  // Program the GeSiCA operational mode
  // ie read out all samples or integrate samples in 3 windows
  // 3 sums or 1 sum may be read
  i2cEnableAll();

  for(Int_t i=0; i<fNSADC; i++){
    Write(EI2CPortSelect,i);       // set port i
    //    usleep(2);
    Write(EIHlPort,i); 
    //    usleep(2);

    // set all adcs to latch or sparse mode
    i2cWriteChk(0x3, 0, (fADCmode & 0x1)); 
    i2cWriteChk(0x3, 1, (fADCmode & 0x1));
    // set latency
    i2cWriteChk(0x0, 0, fLatency);
    i2cWriteChk(0x0, 1, fLatency);
    // set n_samples @ adcs
    i2cWriteChk(0x1, 0, fNSample);
    i2cWriteChk(0x1, 1, fNSample);

  }
  i2cDisableAll();
  printf(" *** SADC operational mode on GeSiCA %s programmed ***\n",
	 this->GetName());
}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::ProgSampleSum(){
  // Program the sample integration boundaries when operating in
  // sample integrating mode
  i2cEnableAll();

  for(Int_t i=0; i<fNSADC; i++){ // loop over all adc-boards
    Write(EIHlPort,i);           // select port i
    //    usleep(2);
    Write(EI2CPortSelect,i);
    //    usleep(2);
 
    i2cWriteChk(4, 0, fSamStart[0]);  // sample 0 pedestal
    i2cWriteChk(4, 1, fSamStart[0]);
    i2cWriteChk(5, 0, fSamWidth[0]);
    i2cWriteChk(5, 1, fSamWidth[0]);
    i2cWriteChk(6, 0, fSamStart[1]);  // sample 1 signal
    i2cWriteChk(6, 1, fSamStart[1]);
    i2cWriteChk(7, 0, fSamWidth[1]);
    i2cWriteChk(7, 1, fSamWidth[1]);
    i2cWriteChk(8, 0, fSamStart[2]);  // sample 2 tail
    i2cWriteChk(8, 1, fSamStart[2]);
    i2cWriteChk(9, 0, fSamWidth[2]);
    i2cWriteChk(9, 1, fSamWidth[2]);
  }

  i2cDisableAll();
  printf(" *** SADC sample sums on GeSiCA %s programmed ***\n",
	 this->GetName());

}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::ProgThresh(){
  // Program the SADC thresholds
  Int_t chip, reg;
  if( fNthresh != fNChannel )
    PrintError("","<Incompatible number thresholds read in>", EErrFatal);

  i2cEnableAll();
  for( Int_t i=0; i<fNthresh; i++ ){
    Write(EIHlPort, &fSADCport[i]);       // select port
    //    usleep(2);
    Write(EI2CPortSelect, &fSADCport[i]);
    //   usleep(2);

    if( fSADCchan[i] < 16 ){
      chip = 0;
      reg = fSADCchan[i] + 0x10;
    }
    else {
      chip = 1;
      reg = fSADCchan[i] - 16 + 0x10;
    }
    i2cWriteChk(reg, chip, fSADCthresh[i]);
    fprintf(fLogStream, "SADCThreshold: %s %4d %4d %4d\n", 
	    GetName(), fSADCport[i], fSADCchan[i], fSADCthresh[i]);
  }
  i2cDisableAll();
  printf(" *** SADC thresholds on GeSiCA %s programmed ***\n",
	 this->GetName());
}

//-----------------------------------------------------------------------------
void TVME_GeSiCA::i2cWait(){
  // Busywait for i2c transfer to be done. Checks bit 0 of SR for
  // max. 2s. Assume usleep(1) = microsecond.
  Int_t datum = 0;  //----------------------------- Baya
  for( Int_t i=0; i<=EI2CTimeout; i++ ){
    datum = Read(EI2CStatus);  //------------------ Baya
    usleep(1);
    if ( !(datum & 0x1) ) return;
  }
  PrintError("","< I2C Status Bit not returned to 0 within 2 sec >", EErrFatal);
}

//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cCheckACK(){
  // Check that the I2C acknowledge bit(3) is cleared
  Int_t datum = 0;
  datum = Read(EI2CStatus);
  if( datum & 0x8 ) 
    printf("< I2C ACK not received >\n");
}  


//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cSwitchPort(unsigned port){
  // Choose port [0..3] for next i2c operations
  if (port < 8){
    Write(EIHlPort,port);       // set broadcast mode
  }
  else {
    printf(" i2cSwitchPort(%x); port number out of range [0..7]", port);
  }
}

//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cEnableAll()
{
  // enable all I2C ports
  Write(EI2CPortSelect,0xff);       // set broadcast mode
  usleep(2);  
  Write(EI2CAddrCtl, 0xe94);       // set broadcast mode
  i2cWait();
}

//------------------------------------------------------------------------------

void TVME_GeSiCA::i2cEnablePort(unsigned port){
  // set i2c enable for port. (Performs i2c_switch port first)
  //
  if (port < 8){
    i2cSwitchPort(port);
    Write(EI2CAddrCtl, 0xe94);
    i2cWait();
  }
  else {
    printf(" i2cEnablePort(%x); port number out of range [0..7]",port);
  }
}

//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cDisableAll(){
  // Disable all ports 0-7
  Write(EI2CPortSelect,0xff);
  usleep(2);
  Write(EI2CAddrCtl, 0xf94);
  i2cWait();
}

//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cDisablePort(unsigned port){

  // set i2c disable for port. (Performs i2cSwitchPort first)
  if (port < 8){
    i2cSwitchPort(port);
    Write(EI2CAddrCtl, 0xf94);
    i2cWait();
  }
  else {
    printf("i2cDisablePort(%x); port number out of range [0..7]", port);
  }
}

//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cWrite(UShort_t address, UShort_t data)
{
  // Write UShort_t (2 bytes) to I2C
  Write(EI2CdwLow,data);
  usleep(2);
  Write(EI2CAddrCtl, address*0x100 + 0x88);
  usleep(2);
  i2cWait();
  i2cCheckACK();
}

//------------------------------------------------------------------------------
void TVME_GeSiCA::i2cWriteB(UShort_t address, UShort_t data)
{
  // Write single byte to I2C
  Write(EI2CdwLow,data);
  usleep(2);
  Write(EI2CAddrCtl, address*0x100 + 0x84);
  usleep(2);
  i2cWait();
  i2cCheckACK();
}

//------------------------------------------------------------------------------
UShort_t TVME_GeSiCA::i2cRead(UShort_t address){
  // read and return UShort_t 
  UInt_t datum =0;  // ------------------- Baya
  Write(EI2CAddrCtl, address*0x100 + 0x98);
  usleep(2);
  i2cWait();                                 // introduce errorhandling here!
  i2cCheckACK();
  datum = Read(EI2CdrLow); // ------------------ Baya
  return datum;
}

//------------------------------------------------------------------------------
UShort_t TVME_GeSiCA::i2cReadB(UShort_t address)
{
  // read single byte and return in UShort_t
  UInt_t datum =0;  // ------------------- Baya
  Write(EI2CAddrCtl, address*0x100 + 0x94);
  usleep(2);
  i2cWait();                             // introduce errorhandling here
  i2cCheckACK();
  datum = Read(EI2CdrLow); // ------------------ Baya
  return datum;
}

//------------------------------------------------------------------------------

void TVME_GeSiCA::i2cReset(){
  // clear registers and reset controller
  UInt_t datum=0;  //  ------------------------------------------ Baya
  Write(EI2CAddrCtl,(UInt_t)0);
  usleep(2);
  Write(EI2CAddrCtl,0x40);
  usleep(2);
  Write(EI2CAddrCtl,(UInt_t)0);
  usleep(2);// added  -------    --------------------------------Baya
  for( Int_t i=0; i<=EI2CTimeout; i++ ){
    datum = Read(EI2CStatus);  // --------------------------------- Baya
    usleep(2);
    if( (datum & 0x7) == 0x6 ){ //------------- Baya, changed JRMA
      usleep(100000);
      return;
    }
  }
  PrintError("","< I2C Reset not acknowledged within 4 sec >", EErrFatal);
}

//------------------------------------------------------------------------------

void TVME_GeSiCA::i2cWriteChk(Int_t i2cAddr, Int_t adc_side,
				UInt_t configVal )
{
  // Write a value to the I2C and read it back to check
  UInt_t readVal =0; 
  i2cReset();
  i2cWait();
  // Seems to be more reliable to read 1st before writing.
  // ------- I2C READ ------------------------------------
  Write(EI2CdwLow, ((i2cAddr >> 5) & 0xff));
  usleep(2);
  Write(EI2CAddrCtl, (0x4000 + ((adc_side << 13) & 0x2000) +
		      ((i2cAddr << 8) & 0x1f00) + 0x84));
  i2cWait();
  Write(EI2CAddrCtl,(0x4000 + ((adc_side << 13) & 0x2000) +
		     ((i2cAddr << 8) & 0x1f00) + 0x98));
  i2cWait();
  readVal = Read(EI2CdrLow); //------------------------------------ Baya  
  readVal = (readVal  & 0xffff);

  //--- ------- I2C WRITE -----------------------------------
  Write(EI2CdwLow,(((configVal << 8) & 0xff00) + ((i2cAddr >> 5) & 0xff)));
  usleep(2); 
  Write(EI2CdwHigh,((configVal >> 8) & 0xff));
  usleep(2); 
  Write(EI2CAddrCtl,(0x4000 + ((adc_side << 13) & 0x2000) +
		     ((i2cAddr << 8) & 0x1f00) + 0x8c));
  i2cWait();

  // Read back and check
  // ------- I2C READ ------------------------------------
  Write(EI2CdwLow,((i2cAddr >> 5) & 0xff));
  usleep(2);
  Write(EI2CAddrCtl,(0x4000 + ((adc_side << 13) & 0x2000) + 
		     ((i2cAddr << 8) & 0x1f00) + 0x84));
  i2cWait();
  Write(EI2CAddrCtl,(0x4000 + ((adc_side << 13) & 0x2000) + 
		     ((i2cAddr << 8) & 0x1f00) + 0x98));
  i2cWait();
  readVal = Read(EI2CdrLow); //------------------------------------ Baya  
  readVal = (readVal  & 0xffff);
  // Final check
  if( readVal != configVal ){
    printf("I2C Write Error: set value = %d, read-back value = %d\n",
	   configVal, readVal);
    
  }
}


