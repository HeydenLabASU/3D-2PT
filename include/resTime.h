typedef struct {
        t_grp   grp;
        real    *resTime;
        int     nTimes;
} t_resTime;

typedef struct {
        real    *time;
        real    *resTime;
        int     nTimes;
        int     cnt;
} t_resTimeTotal;

int allocResTime(t_resTime *resTime,int nTimes);

int initResTime(t_grp grp,t_resTime *resTime);

int initResTimeTotal(t_resTimeTotal *resTimeTotal,int nTimes);

int reinitResTime(t_grp grp,t_resTime *resTime,t_resTimeTotal *resTimeTotal);

int doResTime(t_atom *atoms,t_box box,int nAtoms,t_mol *mols,int nMols,t_statGrpDef *statGrpDef,int nStat,t_dynGrpDef dynGrpDef,t_grp *grps,t_resTime *resTime,int step);

