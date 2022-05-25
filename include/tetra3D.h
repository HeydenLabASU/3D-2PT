typedef struct {
	int 	*coordList;
	int	nCoordList;
} t_tetra;

int allocTetra(t_atom *atoms,int nAtoms,t_mol *mols,t_tetra *tetra);

int gridPos(t_vecR crd,t_vecR ori,real dX,real dY,real dZ,int *pos,t_box box,int *ng);

int doTetra(t_atom *atoms,int nAtoms,t_box box,t_grp centerGrp,t_grp coordGrp,t_tetra *tetra,t_vecR gridOri,int *ng,real dgrid,double *vox,int *vox2,int *vox3);

