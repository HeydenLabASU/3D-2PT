#include "GMXtrrio.h"

int saveOpenRead(FILE **io,char *fn);

int saveOpenReadBin(FILE **io,char *fn);
	
int getFormat(char *coord,int *format);

int openInput(int format,int needVelo,FILE **cin,FILE **vin,XDR *xdrin,char *coord,char *veloxyz);

int closeInput(int format,int needVelo,FILE *cin,FILE *vin,XDR *xdrin);

int checkGROInput(char *buf,t_atom atom,int frame,int element);

int readGRO(FILE *in,t_atom *atoms,int natoms,int needVelo,real *time,int frame,t_box *box);

int readTRR(XDR *in,t_atom *atoms,int natoms,int needVelo,real *time,int frame,t_box *box);

int checkXYZInput(char *buf,t_atom atom,int frame,int element);

int readXYZcrd(FILE *in,t_atom *atoms,int natoms,real *time,int frame);

int readXYZcrdWC(FILE *in,t_atom *atoms,int natoms,t_vecR *WCcrd,int nWC,real *time,int frame);

int readXYZvel(FILE *in,t_atom *atoms,int natoms,int frame);

int readXYZforce(FILE *in,t_atom *atoms,int natoms,int frame);

int readXYZbox(FILE *in,t_box *box,int frame);

int readAMBERcrd(FILE *in,t_atom *atoms,int natoms,real *time,int frame,real dt,t_box *box);

int readAMBERvel(FILE *in,t_atom *atoms,int natoms,int frame);

int getInt(FILE *in,int *i);

int getString(FILE *in,char *buf,int bufLen);

int getReal(FILE *in,real *r);

int getDouble(FILE *in,double *d);

int readDCDheader(FILE *in,int *nFrames,real *time,int nAtoms,real *dt);

int readXST(FILE *in,t_box *box);

int readDCDcrd(FILE *in,t_atom *atoms,int nAtoms,real *time,int frame,real dt,t_box *box);

int readDCDvel(FILE *in,t_atom *atoms,int nAtoms,int frame,real dt);

int prepTRJinput(char *fnCrd,char *fnVel,FILE **cin,FILE **vin,XDR *xdrin,int *format,int *veloFreq,real *time,real *dt,int *nFrames,int needVelo,int nAtoms,int argc);

int readTRJ(int frame,int format,FILE *cin,FILE *vin,XDR *xdrin,real dt,real *time,t_box *box,int nAtoms,int nSites,int needVelo,int veloFreq,int wc,t_atom *atoms,t_vecR *siteCrds);

