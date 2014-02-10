********************* scan_setup.d ******************************************************************************
*
* A comment on comment lines: Information comment lines (like this one) begin with a *
* Lines of data which could be used, or have been used but left in as examples are commented out with a #.
*
*
* The first 3 lines contain info about movement,units and conversion factors (to goni counts).
* Shouldn't really need to change this
*									
*Axis_number		1		2		3		4		5		6
*
Movement_direction	Azi-rot		H-rot		Z-trans		V-rot		X-trans		Mol-rot
units			deg		mrad		mm		mrad		mm		deg
lower_limit(units)	-720		-174.0		-2.0		-174.0		0		-10
upper_limit(units)	720.0		174.0		2.0		174.0		150.0		10.0
units2counts		100.0		57.2967		10000.0		57.2967		1000.0		1000.0
*
* The next line gives the zero positions (IN COUNTS) on all axes. For rotation axes 1,2 and 4 these would generally
* be found from scans. For 6, choose 0 as the default. For translations 3, 5, use the "goni" program to put
* the diamond in position, then go into Expert mode and do commands 3tp and 5tp to get the zero positions. Really
* shouldn't need changed from the default settings I've put in.
*
* The zero scans should be used to update 1,2,4 zero numbers.
*
#zero_pos_in_counts	100		-120		0		77		1840		17
#zero_pos_in_counts      -150    	-1384  		0     		839     	1840     	17
#zero_pos_in_counts      250    		-1282   	0    	 	1957      	1840            17 
*zero_pos_in_counts      0	 0	0	0	1840	17
*zero_pos_in_counts      200    -1184   0     1528   1840	17 
zero_pos_in_counts	200	-1303	0	1140	1840	0	
*
*
* Scan information comes in pairs of lines: 
* The first line gives the desired positions on all other axes during the scan (the actual scan axis position
* will be ignored here, since it comes in the second line.   
* The second line gives: scan axis no, start_pos(units), step_size(units), no of steps.
* ALL NUMBERS ARE IN UNITS APPROPRIATE TO THE AXES (SEE ABOVE).
*
* V-rot ZERO SCAN *************************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
#			0		60		0		0		0		0
*
* 	axis_no		start(mrad)	step_size(mrad)		no_of_steps
#	4		-2		1			34	
*******************************************************************************************************************
*
* H-rot ZERO SCAN *************************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
#			0		80		0		15		0		0
*
* 	axis_no		start(mrad)	step_size(mrad)		no_of_steps
#	2		80		2			1	
*******************************************************************************************************************
*
* A-rot ZERO SCAN *************************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
#			0		60		0		0		0		0
*
* 	axis_no		start(deg)	step_size(deg)		no_of_steps
#	1		0		2			180	
*******************************************************************************************************************
*
* PSEUDO_AZI ZERO SCAN *************************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
			0		80		0	80		0		0
*
* 	axis_no		start(deg)	step_size(deg)		no_of_steps
	7		0		2			180	
*******************************************************************************************************************
*
* V-rot ENERGY SCAN ***********************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
*			0		80		0		0		0		0
*
* 	axis_no		start(mrad)	step_size(mrad)		no_of_steps
*	4		-2		1.0			20	
*******************************************************************************************************************
*
* H-rot ENERGY SCAN ***********************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
*			0		0		0		80		0		0
*
* 	axis_no		start(mrad)	step_size(mrad)		no_of_steps
*	2		-2		1.0			20	
*******************************************************************************************************************
*
* X-trans BREM SCAN ***********************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
#			0		60		0		22		0		0
*
* 	axis_no		start(mm)	step_size(mm)		no_of_steps
#	5		-1		0.1			20	
*******************************************************************************************************************
*
* Z-trans BREM SCAN ***********************************************************************************************
* other positions	1(Azi)		2(H-rot)	3(Z-trans)	4(V-rot)	5(X-trans)	6(Moeller)	
* in units of ...	deg		mrad		mm		mrad		mm		deg
#			0		60		0		22		0		0
*
* 	axis_no		start(mm)	step_size(mm)		no_of_steps
#	3		0.0		0.1			12	
*******************************************************************************************************************
