int allocGrps(t_grpDef grpDef,int nAtoms,int nMols,t_grp **grps);

int checkName(t_nameRule nameRule,char *name,int *res);

int checkNr(int *sel,int nSel,int nr,int *res);

int checkNrRange(int *range,int nr,int *res);

int checkStaticGrp(int idx,t_atom atom,t_mol mol,t_statGrpDef statGrpDef);

int getMolsFromAtoms(t_grp *grp,t_atom *atoms);

int makeStaticGrps(t_atom *atoms,int nAtoms,t_mol *mols,int nMols,t_statGrpDef *statGrpDef,int nStat,t_grp *grps);

int checkDistToCrdAndXYZ(t_vecR trial,t_vecR crd,t_dynGrpDef dynDef,t_box box,int *res);

int checkXYZ(t_vecR trial,t_dynGrpDef dynDef,int *res);

int checkMolAtomDistToCrdAndXYZandDef(t_vecR crd,t_mol mol,t_atom *atoms,t_box box,t_dynGrpDef dynDef,t_statGrpDef preDef,int *res);

int checkMolAtomXYZandDef(t_mol mol,t_atom *atoms,t_dynGrpDef dynDef,t_statGrpDef preDef,int *res);

/*int getGrpCOM(t_atom *atoms,t_grp *grp,t_box box);*/

int distMolsAtomToCrdAndCheckDef(t_vecR crd,t_mol mol,t_atom *atoms,t_box box,t_statGrpDef preDef,real *res);

int distCrdToGrp(t_grp grp,t_vecR trial,t_atom *atoms,t_box box,real *res);

int distMolsAtomToGrpAndCheckDef(t_grp grp,t_mol mol,t_atom *atoms,t_box box,t_statGrpDef preDef,real *res);

int dynAtomsByDistToCrd(t_vecR crd,t_dynGrpDef dynDef,t_atom *atoms,t_box box,t_grp pre,t_grp *res);

int dynAtomsByDistToGrp(t_grp grp,t_dynGrpDef dynDef,t_atom *atoms,t_box box,t_grp pre,t_grp *res);

int dynAtomsByXYZ(t_dynGrpDef dynDef,t_atom *atoms,t_box box,t_grp pre,t_grp *res);

int dynMolsCOMByDistToCrd(t_vecR crd,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsCOMByDistToGrp(t_grp grp,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsCOMByXYZ(t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsAtomByDistToCrd(t_vecR crd,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsAtomByDistToGrp(t_grp grp,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsAtomByXYZ(t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int sortDistances(real *distances,int ndist,int *idxs,int nSel);

int dynAtomsNSelToCrd(t_vecR crd,int nSel,t_atom *atoms,t_box box,t_grp pre,t_grp *res);

int dynAtomsNSelToGrp(t_grp grp,int nSel,t_atom *atoms,t_box box,t_grp pre,t_grp *res);

int dynMolsCOMnSelToCrd(t_vecR crd,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsCOMnSelToGrp(t_grp grp,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsAtomNSelToCrd(t_vecR crd,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

int dynMolsAtomNSelToGrp(t_grp grp,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res);

/* substuted by routine in mol.h */
/*int getMolCOM(t_atom *atoms,t_mol *mol);*/

int makeDynGrps(t_atom *atoms,t_box box,int nAtoms,t_mol *mols,int nMols,t_statGrpDef *statGrpDef,int nStat,t_dynGrpDef *dynGrpDef,int nDyn,t_grp *grps);

int makeCombGrps(t_grp *grps,t_combGrpDef *combDef,int nStat,int nDyn,int nComb,t_atom *atoms,int nAtoms);

int updateGrps(t_atom *atoms,int nAtoms,t_mol *mols,int nMols,t_box box,t_grpDef grpDef,t_grp *grps);

