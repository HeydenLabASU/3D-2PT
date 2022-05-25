typedef struct {
        t_grp   grp;
        t_vecR  *ref;
        t_vecR  *last;
        real    *msd;
        int     nTimes;
} t_msd;

typedef struct {
        real    *time;
        real    *msd;
        int     nTimes;
        int     cnt;
} t_msdTotal;

int allocMSD(t_msd *msd,int nTimes);

int initMSD(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_msd *msd);

int initMSDtotal(t_msdTotal *msdTotal,int nTimes);

int reinitMSD(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_msd *msd,t_msdTotal *msdTotal);

int doMSD(char *mode,t_atom *atoms,t_mol *mols,t_box box,t_msd *msd,int step);

