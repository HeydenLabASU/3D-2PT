typedef struct {
	int 	*coordList;
	int	nCoordList;
	real    *order;
	int	*bin;
	double	aver;
	int	cnt;
} t_tetra;

int allocTetra(t_atom *atoms,int nAtoms,t_mol *mols,t_tetra *tetra);

int doTetra(t_atom *atoms,int nAtoms,t_box box,t_grp centerGrp,t_grp coordGrp,t_tetra *tetra);

