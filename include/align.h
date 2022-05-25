typedef real            rvec[3];

typedef real            matrix[3][3];

static inline void clear_mat(matrix a)
{
  const real nul=0.0;

  a[0][0]=a[0][1]=a[0][2]=nul;
  a[1][0]=a[1][1]=a[1][2]=nul;
  a[2][0]=a[2][1]=a[2][2]=nul;
}

static inline void oprod(const rvec a,const rvec b,rvec c)
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

void jacobi(double **a,int n,double d[],double **v,int *nrot);

void calc_fit_R(int natoms,real *w_rls,rvec *xp,rvec *x,matrix R);

int getAlignRef(t_grp *grp,t_atom *atoms,t_box box,real **w_rls,rvec **refCrd,t_vecR *refCOM,rvec **x);

int alignGroup(t_grp *grp,t_atom *atoms,int nAtoms,t_box *box,real *w_rls,rvec *refCrd,t_vecR refCOM,rvec *x);

int alignGroupVel(t_grp *grp,t_atom *atoms,int nAtoms,t_box *box,real *w_rls,rvec *refCrd,t_vecR refCOM,rvec *x);

int alignGroupNoRot(t_grp *grp,t_atom *atoms,int nAtoms,t_box box,t_vecR refCOM);

int alignGroupNoRotVel(t_grp *grp,t_atom *atoms,int nAtoms,t_box box,t_vecR refCOM);

