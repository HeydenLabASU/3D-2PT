#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../include/dataTypes.h"
#include "../../include/geo.h"
#include "../../include/alloc.h"
#include "../../include/fatal.h"
#include "../../include/qsort.h"
#include "../../include/grps.h"

int allocGrps(t_grpDef grpDef,int nAtoms,int nMols,t_grp **grps)
{
	int nGrps,i;

	nGrps=grpDef.nStat+grpDef.nDyn+grpDef.nComb;
	grps[0]=(t_grp*)save_malloc(nGrps*sizeof(t_grp));

	for(i=grpDef.nStat;i<nGrps;i++)
	{
		allocInts(&grps[0][i].atoms,nAtoms);
		allocInts(&grps[0][i].mols,nMols);
	}

	return 0;
}

int checkName(t_nameRule nameRule,char *name,int *res)
{
	int i,j,k;
	int get,getNot,getBeginWC,getNotBeginWC,getEndWC,getNotEndWC;

	get=0;
	getNot=0;
	getBeginWC=0;
	getNotBeginWC=0;
	getEndWC=0;
	getNotEndWC=0;

	k=strlen(name);

	for(i=0;i<nameRule.nGet;i++)
	{
		if(strcmp(nameRule.get[i],name)==0) get=1;
	}
	for(i=0;i<nameRule.nGetNot;i++)
	{
		if(strcmp(nameRule.getNot[i],name)==0) getNot=1;
	}
	for(i=0;i<nameRule.nGetBeginWC;i++)
	{
		j=nameRule.getBeginWClen[i];
		if(k>j && strcmp(nameRule.getBeginWC[i],&name[k-j])==0) getBeginWC=1;
	}
	for(i=0;i<nameRule.nGetNotBeginWC;i++)
	{
		j=nameRule.getNotBeginWClen[i];
		if(k>j && strcmp(nameRule.getNotBeginWC[i],&name[k-j])==0) getNotBeginWC=1;
	}
	for(i=0;i<nameRule.nGetEndWC;i++)
	{
		j=nameRule.getEndWClen[i];
		if(k>j && strncmp(name,nameRule.getEndWC[i],j)==0) getEndWC=1;
	}
	for(i=0;i<nameRule.nGetNotEndWC;i++)
	{
		j=nameRule.getNotEndWClen[i];
		if(k>j && strncmp(name,nameRule.getNotEndWC[i],j)==0) getNotEndWC=1;
	}

	if((get==1 || getBeginWC==1 || getEndWC==1) && (getNot==0 && getNotBeginWC==0 && getNotEndWC==0)) res[0]=1;
	else if(nameRule.nGet==0 && nameRule.nGetBeginWC==0 && nameRule.nGetEndWC==0 && getNot==0 && getNotBeginWC==0 && getNotEndWC==0) res[0]=1;
	else res[0]=0;

	return 0;
}

int checkNr(int *sel,int nSel,int nr,int *res)
{
	int i;

	res[0]=0;
	for(i=0;i<nSel;i++)
	{
		if(nr==sel[i]) res[0]=1;
	}

	return 0;
}

int checkNrRange(int *range,int nr,int *res)
{
	int a,b;

	if(nr>=range[0]) a=1;
	else a=0;
	if(range[1]==-1 || nr<=range[1]) b=1;
	else b=0;

	if(a==1 && b==1) res[0]=1;
	else res[0]=0;

	return 0;
}

int checkStaticGrp(int idx,t_atom atom,t_mol mol,t_statGrpDef statGrpDef)
{
	int chain,resName,atomName,resNr,resNrRange,atomNr,atomNrRange,res;

	if(statGrpDef.nSelChain!=0) checkName(statGrpDef.chainRule,mol.chain,&chain);
        else chain=1;
	if(statGrpDef.nSelResName!=0) checkName(statGrpDef.resRule,atom.resName,&resName);
	else resName=1;
	if(statGrpDef.nSelAtomName!=0) checkName(statGrpDef.atomRule,atom.atomName,&atomName);
	else atomName=1;
	if(statGrpDef.nResNr!=0) checkNr(statGrpDef.selResNr,statGrpDef.nResNr,atom.resNr,&resNr);
	else resNr=1;
	if(statGrpDef.doResNrRange==1) checkNrRange(statGrpDef.resNrRange,atom.resNr,&resNrRange);
	else resNrRange=1;
	if(statGrpDef.nAtomNr!=0) checkNr(statGrpDef.selAtomNr,statGrpDef.nAtomNr,idx,&atomNr);
	else atomNr=1;
	if(statGrpDef.doAtomNrRange==1) checkNrRange(statGrpDef.atomNrRange,idx,&atomNrRange);
	else atomNrRange=1;

	if(chain==1 && resName==1 && atomName==1 && resNr==1 && resNrRange==1 && atomNr==1 && atomNrRange==1) res=1;
	else res=0;

	return res;
}

int getMolsFromAtoms(t_grp *grp,t_atom *atoms)
{
	int i,cur;

	cur=-1;
	grp[0].nMols=0;
	for(i=0;i<grp[0].nAtoms;i++)
	{
		if(atoms[grp[0].atoms[i]].resNr!=cur)
		{
			cur=atoms[grp[0].atoms[i]].resNr;
			grp[0].mols[grp[0].nMols]=cur;
			grp[0].nMols++;
		}
	}

	return 0;
}

int makeStaticGrps(t_atom *atoms,int nAtoms,t_mol *mols,int nMols,t_statGrpDef *statGrpDef,int nStat,t_grp *grps)
{
	int i,j;
	int chain,resName,atomName,resNr,resNrRange,atomNr,atomNrRange;
	t_grp tmp;

	allocInts(&tmp.atoms,nAtoms);
	allocInts(&tmp.mols,nMols);

	for(i=0;i<nStat;i++)
	{
		printf("creating static group #%d ... ",i); fflush(stdout);
		tmp.nAtoms=0;
		for(j=0;j<nAtoms;j++)
		{
			if(checkStaticGrp(j,atoms[j],mols[atoms[j].resNr],statGrpDef[i])==1)
			{
				tmp.atoms[tmp.nAtoms]=j;
				tmp.nAtoms++;
			}
		}
		getMolsFromAtoms(&tmp,atoms);

		strcpy(grps[i].name,statGrpDef[i].name);
		allocInts(&grps[i].atoms,tmp.nAtoms);
		grps[i].nAtoms=tmp.nAtoms;
		for(j=0;j<tmp.nAtoms;j++)
		{
			grps[i].atoms[j]=tmp.atoms[j];
		}

		allocInts(&grps[i].mols,tmp.nMols);
		grps[i].nMols=tmp.nMols;
		for(j=0;j<tmp.nMols;j++)
		{
			grps[i].mols[j]=tmp.mols[j];
		}
		printf("done: %s has %d atoms in %d molecules\n",grps[i].name,grps[i].nAtoms,grps[i].nMols); fflush(stdout);
	}
	free(tmp.atoms);
	free(tmp.mols);

	return 0;
}

int checkDistToCrdAndXYZ(t_vecR trial,t_vecR crd,t_dynGrpDef dynDef,t_box box,int *res)
{
	real dSq;
	int dist,x,y,z;

	dist=0;
	if(dynDef.doX==1) x=0; else x=1;
	if(dynDef.doY==1) y=0; else y=1;
	if(dynDef.doZ==1) z=0; else z=1;
	distSqPBC(crd,trial,box,&dSq);
	if(dSq>=dynDef.distSqRange[0] && dSq<=dynDef.distSqRange[1]) dist=1;
	else dist=0;
	if(dist==1 && dynDef.doX==1)
	{
		if(trial.x>=dynDef.xRange[0] && trial.x<=dynDef.xRange[1]) x=1;
		else x=0;
	}
	if(dist==1 && x==1 && dynDef.doY==1)
        {
                if(trial.y>=dynDef.yRange[0] && trial.y<=dynDef.yRange[1]) y=1;
                else y=0;
        }
	if(dist==1 && x==1 && y==1 && dynDef.doZ==1)
	{
		if(trial.z>=dynDef.zRange[0] && trial.z<=dynDef.zRange[1]) z=1;
		else z=0;
	}
	if(dist==1 && x==1 && y==1 && z==1) res[0]=1;
	else res[0]=0;

	return 0;
}

int checkXYZ(t_vecR trial,t_dynGrpDef dynDef,int *res)
{
	int x,y,z;

        if(dynDef.doX==1) x=0; else x=1;
        if(dynDef.doY==1) y=0; else y=1;
        if(dynDef.doZ==1) z=0; else z=1;
        if(dynDef.doX==1)
        {
                if(trial.x>=dynDef.xRange[0] && trial.x<=dynDef.xRange[1]) x=1;
                else x=0;
        }
        if(x==1 && dynDef.doY==1)
        {
                if(trial.y>=dynDef.yRange[0] && trial.y<=dynDef.yRange[1]) y=1;
                else y=0;
        }
        if(x==1 && y==1 && dynDef.doZ==1)
        {
                if(trial.z>=dynDef.zRange[0] && trial.z<=dynDef.zRange[1]) z=1;
                else z=0;
        }
        if(x==1 && y==1 && z==1) res[0]=1;
        else res[0]=0;

        return 0;
}

int checkMolAtomDistToCrdAndXYZandDef(t_vecR crd,t_mol mol,t_atom *atoms,t_box box,t_dynGrpDef dynDef,t_statGrpDef preDef,int *res)
{
        int i=0;
        int a=0;

        while(a==0 && i<mol.nAtoms)
        {
                if(checkStaticGrp(mol.atoms[i],atoms[mol.atoms[i]],mol,preDef)==1) checkDistToCrdAndXYZ(atoms[mol.atoms[i]].crd,crd,dynDef,box,&a);
                i++;
        }
        res[0]=a;

        return 0;
}

int checkMolAtomXYZandDef(t_mol mol,t_atom *atoms,t_dynGrpDef dynDef,t_statGrpDef preDef,int *res)
{
        int i=0;
        int a=0;

        while(a==0 && i<mol.nAtoms)
        {
                if(checkStaticGrp(mol.atoms[i],atoms[mol.atoms[i]],mol,preDef)==1) checkXYZ(atoms[mol.atoms[i]].crd,dynDef,&a);
                i++;
        }
        res[0]=a;

        return 0;
}

/*int getGrpCOM(t_atom *atoms,t_grp *grp,t_box box)
{
	int i;
	real totMass;
	t_vecR tmp;
	real max=-1.0;
	real dist;

	totMass=0.0;

	grp[0].COM.x=0.0; grp[0].COM.y=0.0; grp[0].COM.z=0.0;
	for(i=0;i<grp[0].nAtoms;i++)
	{
		vecRAddScaled(grp[0].COM,atoms[grp[0].atoms[i]].mass,atoms[grp[0].atoms[i]].crd,&tmp);
		vecRCpy(tmp,&grp[0].COM);
		totMass+=atoms[grp[0].atoms[i]].mass;
	}

	if(totMass!=0.0)
	{
		vecRScal(1.0/totMass,grp[0].COM,&tmp);
		vecRCpy(tmp,&grp[0].COM);
	}

	for(i=0;i<grp[0].nAtoms;i++)
	{
		distRSq(atoms[grp[0].atoms[i]].crd,grp[0].COM,&dist);
		if(dist>max) max=dist;
	}
	grp[0].spread=sqrt(max);

	return 0;
}*/

int distMolsAtomToCrdAndCheckDef(t_vecR crd,t_mol mol,t_atom *atoms,t_box box,t_statGrpDef preDef,real *res)
{
        int i;
        real min,dSq;

        min=999.9;
        for(i=0;i<mol.nAtoms;i++)
        {
                if(checkStaticGrp(mol.atoms[i],atoms[mol.atoms[i]],mol,preDef)==1)
                {
                        distSqPBC(crd,atoms[mol.atoms[0]].crd,box,&dSq);
                        if(dSq<min) min=dSq;
                }
        }

        res[0]=min;

        return 0;
}

int distCrdToGrp(t_grp grp,t_vecR trial,t_atom *atoms,t_box box,real *res)
{
        int i;
        real min,dSq;

        min=999.9;
        for(i=0;i<grp.nAtoms;i++)
        {
                distSqPBC(atoms[grp.atoms[i]].crd,trial,box,&dSq);
                if(dSq<min) min=dSq;
        }

        res[0]=min;

        return 0;
}

int distMolsAtomToGrpAndCheckDef(t_grp grp,t_mol mol,t_atom *atoms,t_box box,t_statGrpDef preDef,real *res)
{
        int i,j;
        real min,dSq;

        min=999.9;
        for(i=0;i<mol.nAtoms;i++)
        {
                if(checkStaticGrp(mol.atoms[i],atoms[mol.atoms[i]],mol,preDef)==1)
                {
                        for(j=0;j<grp.nAtoms;j++)
                        {
                                distSqPBC(atoms[grp.atoms[j]].crd,atoms[mol.atoms[0]].crd,box,&dSq);
                                if(dSq<min) min=dSq;
                        }
                }
        }

        res[0]=min;

        return 0;
}

int dynAtomsByDistToCrd(t_vecR crd,t_dynGrpDef dynDef,t_atom *atoms,t_box box,t_grp pre,t_grp *res)
{
	int i,j;

	res[0].nAtoms=0;
	for(i=0;i<pre.nAtoms;i++)
	{
		checkDistToCrdAndXYZ(atoms[pre.atoms[i]].crd,crd,dynDef,box,&j);
		if(j==1)
		{
			res[0].atoms[res[0].nAtoms]=pre.atoms[i];
			res[0].nAtoms++;
		}
	}
	getMolsFromAtoms(res,atoms);

	getGrpCOM(atoms,res);

	return 0;
}

int dynAtomsByDistToGrp(t_grp grp,t_dynGrpDef dynDef,t_atom *atoms,t_box box,t_grp pre,t_grp *res)
{
        int i,j,k;
	real check;

        res[0].nAtoms=0;
        for(i=0;i<pre.nAtoms;i++)
        {
                j=0;
		distPBC(atoms[pre.atoms[i]].crd,grp.COM,box,&check);
		if(check<grp.spread+dynDef.distRange[1])
		{
			if(dynDef.distRange[0]==0.0)
			{
                		k=0;
                		while(j==0 && k<grp.nAtoms)
                		{
                		        checkDistToCrdAndXYZ(atoms[pre.atoms[i]].crd,atoms[grp.atoms[k]].crd,dynDef,box,&j);
                		        k++;
                		}
			} else {
				checkXYZ(atoms[pre.atoms[i]].crd,dynDef,&k);
				if(k==1)
				{
					distCrdToGrp(grp,atoms[pre.atoms[i]].crd,atoms,box,&check);
					if(check>=dynDef.distSqRange[0] && check<=dynDef.distSqRange[1]) j=1;
				}
			}
		}
                if(j==1)
                {
                        res[0].atoms[res[0].nAtoms]=pre.atoms[i];
                        res[0].nAtoms++;
                }
        }
        getMolsFromAtoms(res,atoms);

        getGrpCOM(atoms,res);

        return 0;
}

int dynAtomsByXYZ(t_dynGrpDef dynDef,t_atom *atoms,t_box box,t_grp pre,t_grp *res)
{
        int i,j;

        res[0].nAtoms=0;
        for(i=0;i<pre.nAtoms;i++)
        {
                checkXYZ(atoms[pre.atoms[i]].crd,dynDef,&j);
                if(j==1)
                {
                        res[0].atoms[res[0].nAtoms]=pre.atoms[i];
                        res[0].nAtoms++;
                }
        }
        getMolsFromAtoms(res,atoms);

        getGrpCOM(atoms,res);

        return 0;
}

int dynMolsCOMByDistToCrd(t_vecR crd,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
	int i,j;
	int idx;

	res[0].nMols=0;
	for(i=0;i<pre.nMols;i++)
	{
		checkDistToCrdAndXYZ(mols[pre.mols[i]].COM,crd,dynDef,box,&j);
                if(j==1)
		{
			res[0].mols[res[0].nMols]=pre.mols[i];
			res[0].nMols++;
		}
	}
	res[0].nAtoms=0;
	for(i=0;i<res[0].nMols;i++)
	{
		for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
		{
			idx=mols[res[0].mols[i]].atoms[j];
			if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
			{
				res[0].atoms[res[0].nAtoms]=idx;
				res[0].nAtoms++;
			}
		}
	}

	getGrpCOM(atoms,res);

	return 0;
}

int dynMolsCOMByDistToGrp(t_grp grp,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
        int i,j,k;
        int idx;
	real check;

        res[0].nMols=0;
        for(i=0;i<pre.nMols;i++)
        {
                j=0;
		distPBC(mols[pre.mols[i]].COM,grp.COM,box,&check);
		if(check<grp.spread+dynDef.distRange[1])
		{
			if(dynDef.distRange[0]==0.0)
			{
				k=0;
                		while(j==0 && k<grp.nAtoms)
                		{
                		        checkDistToCrdAndXYZ(mols[pre.mols[i]].COM,atoms[grp.atoms[k]].crd,dynDef,box,&j);
                		        k++;
                		}
			} else {
				/*checkXYZ(atoms[pre.atoms[i]].crd,dynDef,&k);*/
				checkXYZ(mols[pre.mols[i]].COM,dynDef,&k);
                                if(k==1)
				{
					distCrdToGrp(grp,mols[pre.mols[i]].COM,atoms,box,&check);
					if(check>=dynDef.distSqRange[0] && check<=dynDef.distSqRange[1]) j=1;
				}
			}
		}
                if(j==1)
                {
                        res[0].mols[res[0].nMols]=pre.mols[i];
                        res[0].nMols++;
                }
        }
        res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

        return 0;
}

int dynMolsCOMByXYZ(t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
        int i,j;
        int idx;

        res[0].nMols=0;
        for(i=0;i<pre.nMols;i++)
        {
                checkXYZ(mols[pre.mols[i]].COM,dynDef,&j);
                if(j==1)
                {
                        res[0].mols[res[0].nMols]=pre.mols[i];
                        res[0].nMols++;
                }
        }
        res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

        return 0;
}

int dynMolsAtomByDistToCrd(t_vecR crd,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
	int i,j;
        int idx;

        res[0].nMols=0;
        for(i=0;i<pre.nMols;i++)
        {
		checkMolAtomDistToCrdAndXYZandDef(crd,mols[pre.mols[i]],atoms,box,dynDef,preDef,&j);	
                if(j==1)
                {
                        res[0].mols[res[0].nMols]=pre.mols[i];
                        res[0].nMols++;
                }
        }
	res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

	getGrpCOM(atoms,res);

	return 0;
}

int dynMolsAtomByDistToGrp(t_grp grp,t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
        int i,j,k,l;
        int idx;
	real check;

        res[0].nMols=0;
        for(i=0;i<pre.nMols;i++)
        {
		j=0;
		l=0;
		while(j==0 && l<mols[pre.mols[i]].nAtoms)
		{	
			idx=mols[pre.mols[i]].atoms[l];
			if(j==0 && checkStaticGrp(idx,atoms[idx],mols[pre.mols[i]],preDef)==1)
			{
				distPBC(atoms[idx].crd,grp.COM,box,&check);
				if(check<grp.spread+dynDef.distRange[1])
				{
					if(dynDef.distRange[0]==0.0)
					{
						k=0;
						while(j==0 && k<grp.nAtoms)
						{
							checkDistToCrdAndXYZ(atoms[grp.atoms[k]].crd,atoms[idx].crd,dynDef,box,&j);
							k++;
						}
					} else {
						checkXYZ(atoms[idx].crd,dynDef,&k);
 		                        	if(k==1)
						{
							distCrdToGrp(grp,atoms[idx].crd,atoms,box,&check);
							if(check>=dynDef.distSqRange[0] && check<=dynDef.distSqRange[1]) j=1;
						}
					}
				}
			}
			l++;
		}
                if(j==1)
                {
                        res[0].mols[res[0].nMols]=pre.mols[i];
                        res[0].nMols++;
                }
        }
        res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

	getGrpCOM(atoms,res);

        return 0;
}

int dynMolsAtomByXYZ(t_dynGrpDef dynDef,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
        int i,j;
        int idx;

        res[0].nMols=0;
        for(i=0;i<pre.nMols;i++)
        {
                checkMolAtomXYZandDef(mols[pre.mols[i]],atoms,dynDef,preDef,&j);
                if(j==1)
                {
                        res[0].mols[res[0].nMols]=pre.mols[i];
                        res[0].nMols++;
                }
        }
        res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

        return 0;
}

int sortDistances(real *distances,int ndist,int *idxs,int nSel)
{
        int i,j,maxpos;
        real maxdist;

        for(i=0;i<nSel;i++)
        {
                idxs[i]=i;
        }

        maxdist=distances[idxs[0]];
        maxpos=0;
        for(i=1;i<nSel;i++)
        {
                if(distances[idxs[i]]>maxdist)
                {
                        maxdist=distances[idxs[i]];
                        maxpos=i;
                }
        }

        for(i=nSel;i<ndist;i++)
        {
                if(distances[i]<maxdist)
                {
                        idxs[maxpos]=i;
                        maxdist=distances[idxs[0]];
                        maxpos=0;
                        for(j=1;j<nSel;j++)
                        {
                                if(distances[idxs[j]]>maxdist)
                                {
                                        maxdist=distances[idxs[j]];
                                        maxpos=j;
                                }
                        }
                }
        }
	qsort(idxs,nSel,sizeof(int),int_cmp);

        return 0;
}

int dynAtomsNSelToCrd(t_vecR crd,int nSel,t_atom *atoms,t_box box,t_grp pre,t_grp *res)
{
        int i;
	real *dists;
	int *idxs;

	allocReals(&dists,pre.nAtoms);
	allocInts(&idxs,nSel);

	for(i=0;i<pre.nAtoms;i++)
	{
		distSqPBC(crd,atoms[pre.atoms[i]].crd,box,&dists[i]);
	}
	sortDistances(dists,pre.nAtoms,idxs,nSel);
	res[0].nAtoms=nSel;
	for(i=0;i<nSel;i++)
	{
		res[0].atoms[i]=pre.atoms[idxs[i]];
	}

	free(dists);
	free(idxs);

        getMolsFromAtoms(res,atoms);

        getGrpCOM(atoms,res);

        return 0;
}

int dynAtomsNSelToGrp(t_grp grp,int nSel,t_atom *atoms,t_box box,t_grp pre,t_grp *res)
{
	int i;
        real *dists;
        int *idxs;

        allocReals(&dists,pre.nAtoms);
        allocInts(&idxs,nSel);

        for(i=0;i<pre.nAtoms;i++)
        {
		distCrdToGrp(grp,atoms[pre.atoms[i]].crd,atoms,box,&dists[i]);
        }
        sortDistances(dists,pre.nAtoms,idxs,nSel);
        res[0].nAtoms=nSel;
        for(i=0;i<nSel;i++)
        {
                res[0].atoms[i]=pre.atoms[idxs[i]];
        }

	free(dists);
        free(idxs);

        getMolsFromAtoms(res,atoms);

        getGrpCOM(atoms,res);

        return 0;
}

int dynMolsCOMnSelToCrd(t_vecR crd,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
	int i,j;
	real *dists;
	int *idxs;
	int idx;

	allocReals(&dists,pre.nMols);
	allocInts(&idxs,nSel);

	for(i=0;i<pre.nMols;i++)
	{
		distSqPBC(crd,mols[pre.mols[i]].COM,box,&dists[i]);
	}
	sortDistances(dists,pre.nMols,idxs,nSel);
	res[0].nMols=nSel;
	for(i=0;i<nSel;i++)
        {
                res[0].mols[i]=pre.mols[idxs[i]];
        }

	free(dists);
	free(idxs);

	res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

	return 0;
}

int dynMolsCOMnSelToGrp(t_grp grp,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
	int i,j;
        real *dists;
        int *idxs;
	int idx;

        allocReals(&dists,pre.nMols);
        allocInts(&idxs,nSel);

        for(i=0;i<pre.nMols;i++)
        {
		distCrdToGrp(grp,mols[pre.mols[i]].COM,atoms,box,&dists[i]);
        }
        sortDistances(dists,pre.nMols,idxs,nSel);
        res[0].nMols=nSel;
        for(i=0;i<nSel;i++)
        {
                res[0].mols[i]=pre.mols[idxs[i]];
        }

        free(dists);
        free(idxs);

        res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

        return 0;
}

int dynMolsAtomNSelToCrd(t_vecR crd,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
	int i,j;
	real *dists;
	int *idxs;
	int idx;

        allocReals(&dists,pre.nMols);
	allocInts(&idxs,nSel);

	for(i=0;i<pre.nMols;i++)
	{
		distMolsAtomToCrdAndCheckDef(crd,mols[pre.mols[i]],atoms,box,preDef,&dists[i]);
	}
	sortDistances(dists,pre.nMols,res[0].mols,nSel);
	res[0].nMols=nSel;
        for(i=0;i<nSel;i++)
        {
                res[0].mols[i]=pre.mols[idxs[i]];
        }

        free(dists);
        free(idxs);

	res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

	return 0;
}

int dynMolsAtomNSelToGrp(t_grp grp,int nSel,t_atom *atoms,t_mol *mols,t_box box,t_grp pre,t_statGrpDef preDef,t_grp *res)
{
	int i,j;
        real *dists;
        int *idxs;
	int idx;

        allocReals(&dists,pre.nMols);
        allocInts(&idxs,nSel);

        for(i=0;i<pre.nMols;i++)
        {
		distMolsAtomToGrpAndCheckDef(grp,mols[pre.mols[i]],atoms,box,preDef,&dists[i]);
        }
        sortDistances(dists,pre.nMols,res[0].mols,nSel);
        res[0].nMols=nSel;
        for(i=0;i<nSel;i++)
        {
                res[0].mols[i]=pre.mols[idxs[i]];
        }

        free(dists);
        free(idxs);

        res[0].nAtoms=0;
        for(i=0;i<res[0].nMols;i++)
        {
                for(j=0;j<mols[res[0].mols[i]].nAtoms;j++)
                {
                        idx=mols[res[0].mols[i]].atoms[j];
                        if(checkStaticGrp(idx,atoms[idx],mols[res[0].mols[i]],preDef)==1)
                        {
                                res[0].atoms[res[0].nAtoms]=idx;
                                res[0].nAtoms++;
                        }
                }
        }

        getGrpCOM(atoms,res);

        return 0;
}

/* substituted by functions in mol.c */
/* COM might be needed more frequently than only upon updating the groups */
/* THIS SHOULD BE DONE IN THE MAIN FUNCTION! */
/*int getMolCOM(t_atom *atoms,t_mol *mol)
{
	int i;
	real totMass=0.0;
	t_vecR tmp;

	mol[0].COM.x=0.0; mol[0].COM.y=0.0; mol[0].COM.z=0.0;
	for(i=0;i<mol[0].nAtoms;i++)
	{
		vecRAddScaled(mol[0].COM,atoms[mol[0].atoms[i]].mass,atoms[mol[0].atoms[i]].crd,&tmp);
                vecRCpy(tmp,&mol[0].COM);
                totMass+=atoms[mol[0].atoms[i]].mass;
        }
        if(totMass!=0.0)
        {
                vecRScal(1.0/totMass,mol[0].COM,&tmp);
                vecRCpy(tmp,&mol[0].COM);
        }

        return 0;
}*/

int makeDynGrps(t_atom *atoms,t_box box,int nAtoms,t_mol *mols,int nMols,t_statGrpDef *statGrpDef,int nStat,t_dynGrpDef *dynGrpDef,int nDyn,t_grp *grps)
{
	int i,idx;
	char err[500];

/*	see comment above */
/*	for(i=0;i<nMols;i++)
	{
		getMolCOM(atoms,&mols[i]);
	}*/

	for(i=0;i<nStat;i++)
	{
		getGrpCOM(atoms,&grps[i]);
	}

	for(i=0;i<nDyn;i++)
	{
		strcpy(grps[nStat+i].name,dynGrpDef[i].name);
		if(dynGrpDef[i].doRef==1 && dynGrpDef[i].doCrd==1)
		{
			sprintf(err,"DYNAMIC GROUP %s:\nUSE OF REFERENCE GROUP AND REFERENCE COORD NOT SUPPORTED!\nYOU CAN DO THIS WITH SEPARATE DYNAMIC GROUP DEFINITIONS USING ONE OF THEM AS PRESELECTION GROUP\n -> makeDynGrps\n -> select.c\n",dynGrpDef[i].name);
			fatal(err);
		}
		if(dynGrpDef[i].doDist==1 && dynGrpDef[i].doN==1)
		{
			sprintf(err,"DYNAMIC GROUP %s:\nUSE OF DISTANCE RANGE AND SELECTION OF %d CLOSEST ATOMS/MOLECULES NOT SUPPORTED!\n-> makeDynGrps\n -> select.c\n",dynGrpDef[i].name,dynGrpDef[i].nSel);
			fatal(err);
		}
		if(dynGrpDef[i].doN==1 && (dynGrpDef[i].doX==1 || dynGrpDef[i].doY==1 && dynGrpDef[i].doZ==1))
		{
			sprintf(err,"DYNAMIC GROUP: %s\nUSE OF X/Y/Z RANGE AND SELECTION OF %d CLOSEST ATOMS/MOLECULES NOT SUPPORTED!\n-> makeDynGrps\n -> select.c\n",dynGrpDef[i].name,dynGrpDef[i].nSel);
			fatal(err);
		}
		if(dynGrpDef[i].doDist==0 && dynGrpDef[i].doN==0)
		{
			dynGrpDef[i].doRef=0;
			dynGrpDef[i].doCrd=0;
		}
		if(dynGrpDef[i].doRef==1 && dynGrpDef[i].doDist==1)
		{
			if(strcmp(dynGrpDef[i].refMode,"COM")==0 && strcmp(dynGrpDef[i].selMode,"atoms")==0)
			{
				dynAtomsByDistToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i],atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
			} else if(strcmp(dynGrpDef[i].refMode,"COM")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
			{ 
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMByDistToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMByDistToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
			} else if(strcmp(dynGrpDef[i].refMode,"COM")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
			{
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomByDistToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
                                idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomByDistToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
			} else if(strcmp(dynGrpDef[i].refMode,"closest_atom")==0 && strcmp(dynGrpDef[i].selMode,"atoms")==0)
                        {
				dynAtomsByDistToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i],atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                        } else if(strcmp(dynGrpDef[i].refMode,"closest_atom")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMByDistToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMByDistToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
                        } else if(strcmp(dynGrpDef[i].refMode,"closest_atom")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomByDistToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomByDistToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
                        }
		} else if(dynGrpDef[i].doCrd==1 && dynGrpDef[i].doDist==1)
		{
			if(strcmp(dynGrpDef[i].selMode,"atoms")==0)
			{
				dynAtomsByDistToCrd(dynGrpDef[i].crd,dynGrpDef[i],atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
			} else if(strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
			{
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMByDistToCrd(dynGrpDef[i].crd,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMByDistToCrd(dynGrpDef[i].crd,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
			} else if(strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
			{
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomByDistToCrd(dynGrpDef[i].crd,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomByDistToCrd(dynGrpDef[i].crd,dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
			}
		} else if(dynGrpDef[i].doRef==1 && dynGrpDef[i].doN==1)
		{
			if(strcmp(dynGrpDef[i].refMode,"COM")==0 && strcmp(dynGrpDef[i].selMode,"atoms")==0)
                        {
				dynAtomsNSelToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i].nSel,atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                        } else if(strcmp(dynGrpDef[i].refMode,"COM")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMnSelToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMnSelToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
                        } else if(strcmp(dynGrpDef[i].refMode,"COM")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomNSelToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                                else {
                                idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomNSelToCrd(grps[dynGrpDef[i].refGrpIdx].COM,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
                                }
                        } else if(strcmp(dynGrpDef[i].refMode,"closest_atom")==0 && strcmp(dynGrpDef[i].selMode,"atoms")==0)
                        {
				dynAtomsNSelToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i].nSel,atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                        } else if(strcmp(dynGrpDef[i].refMode,"closest_atom")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMnSelToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                                else {
                                idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMnSelToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
                                }
                        } else if(strcmp(dynGrpDef[i].refMode,"closest_atom")==0 && strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomNSelToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                                else {
                                idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomNSelToGrp(grps[dynGrpDef[i].refGrpIdx],dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
                                }
                        }
		} else if(dynGrpDef[i].doCrd==1 && dynGrpDef[i].doN==1)
		{
			if(strcmp(dynGrpDef[i].selMode,"atoms")==0)
                        {
				dynAtomsNSelToCrd(dynGrpDef[i].crd,dynGrpDef[i].nSel,atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                        } else if(strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMnSelToCrd(dynGrpDef[i].crd,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                                else {
                                idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMnSelToCrd(dynGrpDef[i].crd,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
                                }
                        } else if(strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
                        {
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomNSelToCrd(dynGrpDef[i].crd,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
                                else {
                                idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomNSelToCrd(dynGrpDef[i].crd,dynGrpDef[i].nSel,atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
                                }
                        }
		} else if(dynGrpDef[i].doX==1 || dynGrpDef[i].doY==1 || dynGrpDef[i].doZ==1)
		{
			if(strcmp(dynGrpDef[i].selMode,"atoms")==0)
			{
				dynAtomsByXYZ(dynGrpDef[i],atoms,box,grps[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
			} else if(strcmp(dynGrpDef[i].selMode,"mol_by_COM")==0)
			{
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsCOMByXYZ(dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsCOMByXYZ(dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
			} else if(strcmp(dynGrpDef[i].selMode,"mol_by_atom")==0)
			{
				if(dynGrpDef[i].preGrpIdx<nStat)
				dynMolsAtomByXYZ(dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[dynGrpDef[i].preGrpIdx],&grps[nStat+i]);
				else {
				idx=dynGrpDef[dynGrpDef[i].preGrpIdx-nStat].preGrpIdx;
				dynMolsAtomByXYZ(dynGrpDef[i],atoms,mols,box,grps[dynGrpDef[i].preGrpIdx],statGrpDef[idx],&grps[nStat+i]);
				}
			}
		} else {
			sprintf(err,"NO VALID SELECTION PARAMETERS FOR DYNAMIC GROUP %s!\n -> makeDynGrps\n -> select.c\n",dynGrpDef[i].name);
			fatal(err);
		}
	/*	printf(" -made dynamic selection: %s  %d\n",grps[nStat+i].name,grps[nStat+i].nAtoms);	*/
	}

	return 0;
}

int makeCombGrps(t_grp *grps,t_combGrpDef *combDef,int nStat,int nDyn,int nComb,t_atom *atoms,int nAtoms)
{
	int *tmp;
	int i,j,k,n,cur;

	allocInts(&tmp,2*nAtoms);

	for(i=0;i<nComb;i++)
	{
		n=0;
		for(j=0;j<combDef[i].nGrpNames;j++)
		{	
			for(k=0;k<grps[combDef[i].grpIdx[j]].nAtoms;k++)
			{
				tmp[n]=grps[combDef[i].grpIdx[j]].atoms[k];
				n++;
			}
		}

		qsort(tmp,n,sizeof(int),int_cmp);

		strcpy(grps[nStat+nDyn+i].name,combDef[i].name);
		cur=-1;
		k=0;
		for(j=0;j<n;j++)
		{
			if(tmp[j]!=cur)
			{
				cur=tmp[j];
				grps[nStat+nDyn+i].atoms[k]=cur;
				k++;
			}
		}
		grps[nStat+nDyn+i].nAtoms=k;
		getMolsFromAtoms(&grps[nStat+nDyn+i],atoms);
	}
	free(tmp);

	return 0;
}

int updateGrps(t_atom *atoms,int nAtoms,t_mol *mols,int nMols,t_box box,t_grpDef grpDef,t_grp *grps)
{	
	int i;

	makeDynGrps(atoms,box,nAtoms,mols,nMols,grpDef.statGrpDef,grpDef.nStat,grpDef.dynGrpDef,grpDef.nDyn,grps);

	makeCombGrps(grps,grpDef.combGrpDef,grpDef.nStat,grpDef.nDyn,grpDef.nComb,atoms,nAtoms);

	return 0;
}

