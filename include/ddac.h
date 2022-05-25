typedef struct {
        t_grp   grp;
        t_vecR  ref;
        real    *ddac;
        int     nTimes;
} t_ddac;

typedef struct {
        real    *time;
        real    *ddac;
        int     nTimes;
        int     cnt;
} t_ddacTotal;

int allocDDAC(t_ddac *ddac,int nTimes);

int getGrpdDM(t_mol *mols,t_grp grp,t_vecR *dDM);

int initDDAC(t_mol *mols,t_grp grp,t_ddac *ddac);

int initDDACtotal(t_ddacTotal *ddacTotal,int nTimes);

int reinitDDAC(t_mol *mols,t_grp grp,t_ddac *ddac,t_ddacTotal *ddacTotal);

int doDDAC(t_mol *mols,t_ddac *ddac,int step);

int ddacSpec(t_ddacTotal *ddacTotal);

int smoothDDacSpec(t_ddacTotal *ddacTotal,real sigma);

