###
Description: CATS readout
File-Name: scratch/CATS 32   	500000	32768   ENet	32768  5019
#File-Name: scratch/CATS 32   	500000	32768 No-Store  32768  5019
#                      Freq
Event-Counters:	       0	0
Start-Run: 0
#
#Control:  Local    0     0
#Control:   ENet  localhost 6020
#Control:   ENetSlave  localhost 6020
Control:    Slave      0 0
#Control:   GUI-Local  0 0
#
#	module	   Name		configFile		BaseAdr		BaseIndex	Nchannel	NBits
Module:	KPH_I686   kphi686_0	kphi686.veto		0		0
Module: KPH_VITEC  vit_0	vitec_0.dat		0000		0		0
Module: CAEN_V965  v965_0	v965_C.noped		24000000   	29032		32		12
Module: CAEN_V1190 v1190_0	v1190_0.dat		29000000   	29192		128		14
# Switch this to VITEC card
IRQ-Ctrl: vit_0
#IRQ-Ctrl: v874_0
Synch-Ctrl: vit_0 400
########### End of File #######################
