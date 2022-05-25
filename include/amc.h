typedef struct {
        t_vecR  angMom;
        real  angMom_proj1;
        real  angMom_proj2;
        real  angMom_proj3;
} t_angMom;

typedef struct {
        t_grp   grp;
	t_vecR  *angMom;
	real    *angMom_proj1;
	real    *angMom_proj2;
	real    *angMom_proj3;
        t_vecR  *ref;
	real	*ref_proj1;
	real	*ref_proj2;
	real	*ref_proj3;
        real    *amc;
        real	*amc_proj1;
	real	*amc_proj2;
	real	*amc_proj3;
        int     nTimes;
} t_amc;

typedef struct {
        real    *time;
        real    *amc;
        real	*amc_proj1;
	real	*amc_proj2;
	real	*amc_proj3;
        int     nTimes;
        int     cnt;
} t_amcTotal;

typedef struct {
	real	*time;
	real	*spectrum;
} t_spec;

int initAMC(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_amc *amc);

int initAMCtotal(t_amcTotal *amcTotal,int nTimes);

int initSpec(t_spec *spec,int nTimes);

int reinitAMC(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_amc *amc,t_amcTotal *amcTotal);

int doAMC(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_amc *amc,int step,t_angMom *AngMom);

int amcDOS(t_amcTotal *amcTotal,t_spec *spec,int k);

int smoothAmcDOS(t_amcTotal *amcTotal,t_spec *spec,real sigma);
