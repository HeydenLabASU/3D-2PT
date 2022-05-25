typedef struct {
        t_grp   grp;
        t_vecR  ref;
        real    *dac;
        int     nTimes;
} t_dac;

typedef struct {
        real    *time;
        real    *dac;
        int     nTimes;
        int     cnt;
} t_dacTotal;

int allocDAC(t_dac *dac,int nTimes);

int getGrpDM(t_mol *mols,t_grp grp,t_vecR *DM);

int initDAC(t_mol *mols,t_grp grp,t_dac *dac);

int initDACtotal(t_dacTotal *dacTotal,int nTimes);

int reinitDAC(t_mol *mols,t_grp grp,t_dac *dac,t_dacTotal *dacTotal);

int doDAC(t_mol *mols,t_dac *dac,int step);

int dacSpec(t_dacTotal *dacTotal);

int smoothDacSpec(t_dacTotal *dacTotal,real sigma);

