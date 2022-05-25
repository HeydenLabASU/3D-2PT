typedef struct {
        t_grp   grp;
        t_vecR  *ref;
        real    *sdac;
        int     nTimes;
} t_sdac;

typedef struct {
        real    *time;
        real    *sdac;
        int     nTimes;
        int     cnt;
} t_sdacTotal;

int allocSDAC(t_sdac *sdac,int nTimes);

int initSDAC(t_mol *mols,t_grp grp,t_sdac *sdac);

int initSDACtotal(t_sdacTotal *sdacTotal,int nTimes);

int reinitSDAC(t_mol *mols,t_grp grp,t_sdac *sdac,t_sdacTotal *sdacTotal);

int doSDAC(t_mol *mols,t_sdac *sdac,int step);

int sdacSpec(t_sdacTotal *sdacTotal);

int smoothSdacSpec(t_sdacTotal *sdacTotal,real sigma);

