typedef struct {
        t_grp   grp;
        t_vecR  *ref;
        real    *vac;
        int     nTimes;
} t_vac;

typedef struct {
        real    *time;
        real    *vac;
        int     nTimes;
        int     cnt;
} t_vacTotal;

int allocVAC(t_vac *vac,int nTimes);

int initVAC(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_vac *vac);

int initVACtotal(t_vacTotal *vacTotal,int nTimes);

int reinitVAC(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_vac *vac,t_vacTotal *vacTotal);

int doVAC(char *mode,t_atom *atoms,t_mol *mols,t_vac *vac,int step);

int vacDOS(t_vacTotal *vacTotal);

int smoothVacDOS(t_vacTotal *vacTotal,real sigma);

