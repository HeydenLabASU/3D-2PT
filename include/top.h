int getLineFromTOP(FILE *in,char *buffer,int *lineCnt);

int allocMolContents(t_mol *dest);

int getIntFromTOP(char *fn,FILE *in,int *dest,char *keyword,int *lineCnt);

int getStringFromTOP(char *fn,FILE *in,char *dest,char *keyword,int *lineCnt);

int getAtomLineFromTOP(char *fn,FILE *in,char *name,real *mass,real *charge,real *pol,int *lineCnt);

int getBondLineFromTOP(char *fn,FILE *in,t_bond *bond,int min,int max,int *lineCnt,int bondOffset);

int getAngleLineFromTOP(char *fn,FILE *in,t_angle *angle,int min,int max,int *lineCnt,int bondOffset);

int getDihedralLineFromTOP(char *fn,FILE *in,t_dihed *dihedral,int min,int max,int *lineCnt,int bondOffset);

int getImproperLineFromTOP(char *fn,FILE *in,t_improp *improper,int min,int max,int *lineCnt,int bondOffset);

int cpyMolTop(t_mol src,t_mol *dest,t_atom *atoms,int copy);

int readTOP(char *fn,t_atom **atoms,int *nAtoms,t_vecR **siteCrds,int *nSites,t_mol **mols,int *nMols);

