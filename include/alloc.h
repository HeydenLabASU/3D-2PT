#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,__FILE__,__LINE__,\
                        (nelem),sizeof(*(ptr)))

#define sfree(ptr) save_free(#ptr,__FILE__,__LINE__,(ptr))

void *save_calloc(char *name,char *file,int line,
                  unsigned nelem,unsigned elsize);

void save_free(char *name,char *file,int line, void *ptr);

void *save_malloc(size_t size);

void *save_realloc(void *ptr,size_t size);

int allocMols(t_mol **mols,int n);

int allocAtoms(t_atom **atoms,int n);

int allocSites(t_site **sites,int n);

int allocBonds(t_bond **bonds,int n);

int allocAngles(t_angle **angles,int n);

int allocDihedrals(t_dihed **dihedrals,int n);

int allocImpropers(t_improp **impropers,int n);

int allocInts(int **ints,int n);

int allocReals(real **reals,int n);

int allocVecRs(t_vecR **vecRs,int n);

int allocStatGrpDefs(t_statGrpDef **statGrpDef,int nStat);

int allocDynGrpDefs(t_dynGrpDef **dynGrpDef,int nDyn);

int allocCombGrpDefs(t_combGrpDef **combGrpDef,int nComb);

