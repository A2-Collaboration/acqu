****************************************************************************
*	KELVIN LABORATORY, UNIVERSITY OF GLASGOW
*	*********** lrs_2373.d *****************
*	
*	Sample parameter file for LeCroy 2373 MLU
*	JRM Annand	10th Sep 1993
****************************************************************************
*	1st section for pattern unit read
*	mode = a for subaddr to be read out
*	mode = x for no readout
*	thresholds entered in hexadecimal please
*
*	mode		lower thresh		upper thresh
	a		0			ffff
	a		0			ffff
*
*	Now set the 2373 operation mode
*	40(hex)....2373 "strobe" mode (latched outputs)
*	50(hex)....2373 inhibit (for setting up only)
*	60(hex)....2373 "pulse" mode (pulse outputs)
*	dummy	mode
	0	60
*
*	Now comes the bit to set up the 2373 RAM
*	The 1st section specifies default value which is fed into all RAM
*	dummy	default data
	0	0
*
*	Now individually load RAM
*	addr	data
	0	0
	1	1
	2	2
	3	3
	4	4
	5	5
*	end of input
