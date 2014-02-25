*	KELVIN LABORATORY, UNIVERSITY OF GLASGOW
*
*	************ lrs_4434.d ***************
*
*	M Sauer	 	5th Oct 1995
*	M Lang         25th Oct 1995	 
*
*	Parameterfile for the scaler LeCroy 4434
*
**********************************************************************
* Write the Command Register of the 4344. This register has the following
* constellation:
* Bit   Function
* 1-5	first address to be read
* 6	enable load
* 7 	clear all 32 scalers
* 8	readout enable
* 9-13  number of channels to be read out
* 16	test the scaler
*
*1F90
*
********************************************************************
*
*	mode = z for subaddr to be read out 16 bit
*	mode = l for subaddr to be read out 24 bit
*
*   ------------------------------------------------------	
*
*	mode		1st channel		last channel
	z  		0			1f
*
