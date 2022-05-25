int vecRCpy(t_vecR r1,t_vecR *r2);

int vecRAdd(t_vecR r1,t_vecR r2,t_vecR *r3);

int vecRSub(t_vecR r1,t_vecR r2,t_vecR *r3);

int vecRMulti(t_vecR r1,t_vecR r2,t_vecR *r3);

int vecRProd(t_vecR r1,t_vecR r2,real *p);

int vecRScal(real s,t_vecR r1,t_vecR *r2);

int vecRAddScaled(t_vecR r1,real s,t_vecR r2,t_vecR *r3);

int vecRCross(t_vecR r1,t_vecR r2,t_vecR *r3);

int vecRNormSq(t_vecR r1,real *dSq);

int vecRNorm(t_vecR r1,real *d);

int distRSq(t_vecR r1,t_vecR r2,real *dSq);

int distR(t_vecR r1,t_vecR r2,real *d);

int invMatR(t_vecR c1,t_vecR c2,t_vecR c3,t_vecR *cInv1,t_vecR *cInv2,t_vecR *cInv3);

int colMatDotVecR(t_vecR c1,t_vecR c2,t_vecR c3,t_vecR r,t_vecR *res);

int rowMatDotVecR(t_vecR r1,t_vecR r2,t_vecR r3,t_vecR r,t_vecR *res);

int makePBCtransMat(t_box *box);

int transXYZToPBCbasis(t_box box,t_vecR r,t_vecR *res);

int transPBCToXYZbasis(t_box box,t_vecR r,t_vecR *res);

int wrapTrajectoryAtoms(t_atom *atoms,int nAtoms,t_box box);

int wrapTrajectoryAtomsGrp(t_atom *atoms,t_grp grp,t_box box);

int wrapTrajectoryMol0(t_atom *atoms,t_mol *mols,int nMols,t_box box);

int wrapTrajectoryMol0Grp(t_atom *atoms,t_mol *mols,t_grp grp,t_box box);

int wrapTrajectoryMolCOM(t_atom *atoms,t_mol *mols,int nMols,t_box box);

int wrapTrajectoryMolCOMGrp(t_atom *atoms,t_mol *mols,t_grp grp,t_box box);

int linkPBCfull(t_vecR *link,t_box box);

int linkPBCfullRaw(t_vecR *link,t_box box);

int linkPBCortho(t_vecR *link,t_box box);

int linkPBC(t_vecR r1,t_vecR r2,t_box box,t_vecR *res);

int fixMols(t_atom *atoms,t_mol *mols,int nMols,t_box box);

int distSqPBC(t_vecR r1,t_vecR r2,t_box box,real *dSq);

int distPBC(t_vecR r1,t_vecR r2,t_box box,real *d);

int distDeriv(t_vecR a,t_vecR da,t_vecR b,t_vecR db,t_box box,real *d);

int getAnglePBC(t_vecR a,t_vecR b,t_vecR center,t_box box,real *angle);

int angleDeriv(t_vecR a,t_vecR da,t_vecR b,t_vecR db,t_vecR c,t_vecR dc,t_box box,real *d);

