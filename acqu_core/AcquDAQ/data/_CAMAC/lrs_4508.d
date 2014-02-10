****************************************************************************
*	KELVIN LABORATORY, UNIVERSITY OF GLASGOW
*	*********** lrs_4508.d *****************
*	
*	Sample parameter file for LeCroy 4508 PLU
*	JRM Annand	10th Sep 1993
****************************************************************************
*	1st section for pattern unit read
*	mode = a for subaddr to be read out
*	mode = x for no readout
*	thresholds entered in hexadecimal please
*
*	mode		lower thresh		upper thresh
	a		0			100
	a		0			100
*
*	Now comes the bit to set up the 4508 RAM
*	This 1st section specifies default values which are fed into
*	RAM's 1 (specify -1 for default) and 2 (specify -2 for default)
*	ram(default)	dummy	default data
	-1		0	1
	-2		0	0
*	ram: 1 or 2, corresponding to section 1 or 2 of module
*	addr:0 - ff(hex) the memory location
*	data:0 - ff(hex) the data loaded into location addr
*
*	ram	addr	data
	1	3	2
*
*	dummy stuff for 2
*
	2	1	1
	2	2	2
	2	3	3
	2	4	4
*	end of input
