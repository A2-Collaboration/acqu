*	INSTITUTE OF NUCLEAR PHYSICS, UNIVERSITY OF MAINZ
*
*	JRM Annand	 1st May 1996
*
*	4300 FERA module pedestal file
*	Applies ONLY for CAMAC readout of module
*	ie module type "LRS4300B" in CAMAC setup file
*
*	16 lines, A = 0-15, one pedestal value (in hex) per line
*
*	This just a test
	0
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
	ff
**************************************************************************
*	end of pedestals
**************************************************************************
*	This line tells system its an ADC (this is a CAMAC read)
*	window set 0 - 7ff
	b	0	7ff
**************************************************************************
*	Now 1 extra line to specify the mode of FERA readout
*	b800(hex)	Compressed CAMAC readout. Pedestal subtract
*			and overflow suppress
*	3800(hex)	Compressed CAMAC readout. Pedestal subtract
*			but no overflow suppress
*	2000(hex)	UNcompressed CAMAC readout. All FERA channels
*			are transferred up to nchan.
*			Use for pedestal calibration or bebugging
*	other		Any other value or lack of this line results in
*			option 8700 being applied as default
*
*	NB the 1st LRS4300B file in the CAMAC list sets whether
*	compressed or uncompressed data readout is chosen
***************************************************************************
*	The optional 2nd parameter specifies the no. of channels to read
***************************************************************************
*	No compression, 16 channels read..remember HEX!!
	b800	10
*	end of input
