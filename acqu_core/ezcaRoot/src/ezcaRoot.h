#ifndef __ROOTEPICS__
#define __ROOTEPICS__


#ifndef __CINT__
//#include <tsDefs.h> this is obsolete since R3.14!
#include <cadef.h>
#include <ezca.h>

int epicsGetStatus(char *pvname,TS_STAMP *timestamp, short *status, short *severity);
int epicsGetWithStatus(char *pvname,char ezcatype, int nelem, void *data_buff, TS_STAMP *timestamp,short *status, short *severity);
#endif

enum EEpicsTypes{kEpicsBYTE, kEpicsSTRING, kEpicsSHORT, kEpicsLONG, kEpicsFLOAT, kEpicsDOUBLE};

//only these will be available through CINT
int epicsGet(char *pvname, char ezcatype,int nelem, void *data_buff);
int epicsGetChar(char *pvname, int nelem, char *data_buff);
int epicsGetString(char *pvname, int nelem, char *data_buff);
int epicsGetShort(char *pvname, int nelem, short *data_buff);
int epicsGetFloat(char *pvname, int nelem, float *data_buff);
int epicsGetDouble(char *pvname, int nelem, double *data_buff);

int epicsPut(char *pvname, char ezcatype,int nelem, void *data_buff);
int epicsPutChar(char *pvname, int nelem, char *data_buff);
int epicsPutString(char *pvname, int nelem, char *data_buff);
int epicsPutShort(char *pvname, int nelem, short *data_buff);
int epicsPutFloat(char *pvname, int nelem, float *data_buff);
int epicsPutDouble(char *pvname, int nelem, double *data_buff);

int epicsGetControlLimits(char *pvname,double *low, double *high);
int epicsGetGraphicLimits(char *pvname,double *low, double *high);
int epicsGetNelem(char *pvname, int *nelem);
int epicsGetPrecision(char *pvname,short *precision);
int epicsGetUnits(char *pvname,char *units); /* units must be at least	EZCA_UNITS_SIZE large */
int epicsPutOldCa(char *pvname, char ezcatype,int nelem, void *data_buff);
int epicsSetMonitor(char *pvname, char ezcatype);

#endif
