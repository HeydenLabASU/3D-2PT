typedef struct {
        t_hb	*ref;
	int	nRef;
        int	*hbc;
	int	*hbs;
        int     nTimes;
} t_hbc;

typedef struct {
        real    *time;
        real    *hbc;
	real	*hbs;
        int     nTimes;
        int     cnt;
} t_hbcTotal;

int allocHBC(t_atom *atoms,int nAtoms,t_mol *mols,t_hbc *hbc,int nTimes);

int checkHB(t_atom *atoms,t_box box,real maxDSq,real minAng,int don,int donH,int acc,real *resDSq,real *resAng);

int addHBsToList(t_atom *atoms,t_box box,real maxDSq,real minAng,t_grp grp1, t_grp grp2,t_hb *hb,int *nHB);

int initHBC(t_atom *atoms,int nAtoms,t_box box,real maxDSq,real minAng,t_grp *grps,int grp1Idx,int grp2Idx,t_hbc *hbc);

int initHBCtotal(t_hbcTotal *hbcTotal,int nTimes);

int reinitHBC(t_atom *atoms,int nAtoms,t_box box,real maxDSq,real minAng,t_grp *grps,int grp1Idx,int grp2Idx,t_hbc *hbc,t_hbcTotal *hbcTotal);

int doHBC(t_atom *atoms,t_box box,real maxDSq,real minAng,t_hbc *hbc,int step);

