#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/dataTypes.h"
#include "../../include/alloc.h"
#include "../../include/top.h"
#include "../../include/io.h"
#include "../../include/mol.h"
#include "../../include/fatal.h"
#include "../../include/job.h"
#include "../../include/select.h"
#include "../../include/geo.h"
#include "../../include/mol.h"
#include "../../include/align.h"

typedef struct {
	int *idx;
	int n;
} t_list;

typedef struct {
        int     nMol;
        double    averMol;
        double    Uxs,Uss;
} t_voxelProp;

int initVoxelProp(t_voxelProp *p) {
        p[0].nMol=0;
	p[0].averMol=0.0;
        p[0].Uxs=0.0;
        p[0].Uss=0.0;
        return 0;
}

int gridOrigin(int *ng,t_vecR dgrid,t_vecR *ori) {
	ori[0].x=(((real)ng[0])/-2.0)*dgrid.x;
	ori[0].y=(((real)ng[1])/-2.0)*dgrid.y;
	ori[0].z=(((real)ng[2])/-2.0)*dgrid.z;
	return 0;
}

int gridPos(t_vecR crd,t_vecR ori,real dX,real dY,real dZ,int *pos,t_box box,int *ng) {
	t_vecR link;
	int check;
	real tmp1;
	int tmp2;

	/*grid ori is the origin of the grid (one of its corners, not the center)*/
	//linkPBC(ori,crd,box,&link);
	vecRSub(ori,crd,&link);
	
	if(link.x<0.0) link.x+=1e-5;
	if(link.y<0.0) link.y+=1e-5;
	if(link.z<0.0) link.z+=1e-5;
	if(link.x<0.0 || link.y<0.0 || link.z<0.0) {
		check=0; /*outside of grid volume*/
	} else {
		tmp1=link.x/dX;
		tmp2=(int)floor(tmp1);
		if(tmp2>=ng[0]) { tmp1-=1e-5; tmp2=(int)floor(tmp1); }
		pos[0]=tmp2;
		if(pos[0]>=ng[0]) {
			check=0; /*outside of grid volume*/
		} else {
			tmp1=link.y/dY;
			tmp2=(int)floor(tmp1);
			if(tmp2>=ng[1]) { tmp1-=1e-5; tmp2=(int)floor(tmp1); }
			pos[1]=tmp2;
			if(pos[1]>=ng[1]) {
				check=0; /*outside of grid volume*/
			} else {
				tmp1=link.z/dZ;
				tmp2=(int)floor(tmp1);
				if(tmp2>=ng[2]) { tmp1-=1e-5; tmp2=(int)floor(tmp1); }
				pos[2]=tmp2;
				if(pos[2]>=ng[2]) {
					check=0; /*outside of grid volume*/
				} else {
					check=1; /*INSIDE THE GRID VOLUME*/
				}
			}
		}
	}
	return check;
}

int makeLJcombParameters(t_atom *atoms,int nAtoms,int *types,double ***sigSqMat,double ***epsMat) {
        int i,j;
        real *sigList,*epsList;
        int nTypes;

        /*preparing a list of atom types and matrices with the combined LJ potential parameters*/
        /*otherwise we would have to recompute them all the time*/
        sigList=(real*)save_malloc(nAtoms*sizeof(real));
        epsList=(real*)save_malloc(nAtoms*sizeof(real));
        nTypes=0;
        for(i=0;i<nAtoms;i++) {
                if(nTypes==0) {
                        sigList[nTypes]=atoms[i].sigmaLJ;
                        epsList[nTypes]=atoms[i].epsLJ;
                        nTypes++;
                } else {
                        j=0;
                        while((atoms[i].sigmaLJ!=sigList[j] || atoms[i].epsLJ!=epsList[j]) && j<nTypes) {
                                j++;
                        }
                        if(j==nTypes) {
                                sigList[nTypes]=atoms[i].sigmaLJ;
                                epsList[nTypes]=atoms[i].epsLJ;
                                nTypes++;
                        }
                }
        }
        for(i=0;i<nAtoms;i++) {
                j=0;
                while((atoms[i].sigmaLJ!=sigList[j] || atoms[i].epsLJ!=epsList[j]) && j<nTypes) {
                        j++;
                }
                if(j==nTypes) fatal("something went wrong with the atom type recognition\n");
                types[i]=j;
        }
        /*now we finally make the matrix of combined LJ parameters*/
        sigSqMat[0]=(double**)save_malloc(nTypes*sizeof(double));
        epsMat[0]=(double**)save_malloc(nTypes*sizeof(double*));
        for(i=0;i<nTypes;i++) {
                sigSqMat[0][i]=(double*)save_malloc(nTypes*sizeof(double));
                epsMat[0][i]=(double*)save_malloc(nTypes*sizeof(double));
        }
        for(i=0;i<nTypes;i++) {
                for(j=0;j<nTypes;j++) {
                        sigSqMat[0][i][j]=sigList[i]*sigList[j];
                        epsMat[0][i][j]=sqrt(epsList[i]*epsList[j]);
                }
        }
        free(sigList);
        free(epsList);
        return 0;
}

int allocPosGridStuff(t_list **posGridX,t_list **posGridS,t_list **posGridNBListX,t_list **posGridNBListS,
	int ***posGridPos,int dimPosGrid,int dimPosGridCube,int nMols) {
	int i,j,k,l;

        posGridX[0]=(t_list*)save_malloc(dimPosGridCube*sizeof(t_list));
        posGridS[0]=(t_list*)save_malloc(dimPosGridCube*sizeof(t_list));
        posGridNBListX[0]=(t_list*)save_malloc(dimPosGridCube*sizeof(t_list));
        posGridNBListS[0]=(t_list*)save_malloc(dimPosGridCube*sizeof(t_list));
        for(i=0;i<dimPosGridCube;i++) {
                posGridX[0][i].idx=(int*)save_malloc(nMols*sizeof(int));
                posGridS[0][i].idx=(int*)save_malloc(nMols*sizeof(int));
                posGridNBListX[0][i].idx=(int*)save_malloc(nMols*sizeof(int));
                posGridNBListS[0][i].idx=(int*)save_malloc(nMols*sizeof(int));
        }
	/*this array of integer vectors translates the 1D array index into voxel coordinates*/
        posGridPos[0]=(int**)save_malloc(dimPosGridCube*sizeof(int*));
        for(i=0;i<dimPosGridCube;i++) posGridPos[0][i]=(int*)save_malloc(3*sizeof(int));
        for(i=0;i<dimPosGrid;i++) { for(j=0;j<dimPosGrid;j++) { for(k=0;k<dimPosGrid;k++) {
                l=k+dimPosGrid*(j+dimPosGrid*i);
                posGridPos[0][l][0]=i;
                posGridPos[0][l][1]=j;
                posGridPos[0][l][2]=k;
        } } }
	return 0;
}

int prepPosGridParameters(t_list *posGridX,t_list *posGridS,t_list *posGridNBListX,t_list *posGridNBListS,
	t_vecR *dPosGrid,t_vecR *oriPosGrid,int *posGridRange,int dimPosGrid,int dimPosGridCube,t_box box,
	real cutOffList,real cutOffListSq,int **posGridNBtransVec,int *posGridNBtransN) {
	int i,j,k,n;
	t_vecR tmp1,tmp2;
	real tmp3;

	for(i=0;i<dimPosGridCube;i++) {
		posGridX[i].n=0;
		posGridS[i].n=0;
		posGridNBListX[i].n=0;
		posGridNBListS[i].n=0;
	}
	vecRCpy(dPosGrid[0],&tmp1);
	vecRCpy(oriPosGrid[0],&tmp2);
	dPosGrid->x=box.a.x/dimPosGrid; dPosGrid->y=box.b.y/dimPosGrid; dPosGrid->z=box.c.z/dimPosGrid;
	oriPosGrid->x=-0.5*box.a.x; oriPosGrid->y=-0.5*box.b.y; oriPosGrid->z=-0.5*box.c.z;
	posGridRange[0]=ceil(cutOffList/dPosGrid->x);
	posGridRange[1]=ceil(cutOffList/dPosGrid->y);
	posGridRange[2]=ceil(cutOffList/dPosGrid->z);
	if(2*posGridRange[0]+1>=dimPosGrid) posGridRange[0]=(dimPosGrid-1)/2;
	if(2*posGridRange[1]+1>=dimPosGrid) posGridRange[1]=(dimPosGrid-1)/2;
	if(2*posGridRange[2]+1>=dimPosGrid) posGridRange[2]=(dimPosGrid-1)/2;

	if(dPosGrid->x!=tmp1.x || dPosGrid->y!=tmp1.y || dPosGrid->z!=tmp1.z || 
		oriPosGrid->x!=tmp2.x || oriPosGrid->y!=tmp2.y || oriPosGrid->z!=tmp2.z) {
		/*set all vectors to 0 0 0, for the first one, this also the final values*/
		for(i=0;i<dimPosGridCube;i++) {
			for(j=0;j<3;j++) {
				posGridNBtransVec[i][j]=0;
			}
		}
	        n=1;
		/*we also takes the axes as default*/
		for(i=0;i<3;i++) {
			for(j=1;j<=posGridRange[i];j++) {
				posGridNBtransVec[n][i]=-1*j; n++;
				posGridNBtransVec[n][i]=   j; n++;
			}
		}
		/*the planes will be tested*/
		for(i=1;i<=posGridRange[0];i++) {
			for(j=1;j<=posGridRange[1];j++) {
				tmp1.x=i-0.5;
				tmp1.y=j-0.5;
				tmp1.z=0.0;
				vecRMulti(tmp1,dPosGrid[0],&tmp2);
				vecRNormSq(tmp2,&tmp3);
				if(tmp3<cutOffListSq) {
					posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][1]=-1*j; n++;
					posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][1]=   j; n++;
					posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][1]=   j; n++;
					posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][1]=-1*j; n++;
				}
			}
		}
		for(i=1;i<=posGridRange[0];i++) {
	                for(j=1;j<=posGridRange[2];j++) {
	                        tmp1.x=i-0.5;
	                        tmp1.y=0.0;
	                        tmp1.z=j-0.5;
	                        vecRMulti(tmp1,dPosGrid[0],&tmp2);
	                        vecRNormSq(tmp2,&tmp3);
	                        if(tmp3<cutOffListSq) {
	                                posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][2]=-1*j; n++;
	                                posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][2]=   j; n++;
	                                posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][2]=   j; n++;
	                                posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][2]=-1*j; n++;
	                        }
	                }
	        }
		for(i=1;i<=posGridRange[1];i++) {
	                for(j=1;j<=posGridRange[2];j++) {
	                        tmp1.x=0.0;
	                        tmp1.y=i-0.5;
	                        tmp1.z=j-0.5;
	                        vecRMulti(tmp1,dPosGrid[0],&tmp2);
	                        vecRNormSq(tmp2,&tmp3);
	                        if(tmp3<cutOffListSq) {
	                                posGridNBtransVec[n][1]=-1*i; posGridNBtransVec[n][2]=-1*j; n++;
	                                posGridNBtransVec[n][1]=   i; posGridNBtransVec[n][2]=   j; n++;
	                                posGridNBtransVec[n][1]=-1*i; posGridNBtransVec[n][2]=   j; n++;
	                                posGridNBtransVec[n][1]=   i; posGridNBtransVec[n][2]=-1*j; n++;
	                        }
	                }
	        }
		/*and now the rest*/
		for(i=1;i<=posGridRange[0];i++) { 
			tmp1.x=i-0.5;
			for(j=1;j<=posGridRange[1];j++) {
				tmp1.y=j-0.5;
				for(k=1;k<=posGridRange[2];k++) {
					tmp1.z=k-0.5;
					vecRMulti(tmp1,dPosGrid[0],&tmp2);
					vecRNormSq(tmp2,&tmp3);
					if(tmp3<cutOffListSq) {
						posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][1]=-1*j; posGridNBtransVec[n][2]=-1*k; n++;
						posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][1]=-1*j; posGridNBtransVec[n][2]=   k; n++;
						posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][1]=   j; posGridNBtransVec[n][2]=-1*k; n++;
						posGridNBtransVec[n][0]=-1*i; posGridNBtransVec[n][1]=   j; posGridNBtransVec[n][2]=   k; n++;
						posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][1]=-1*j; posGridNBtransVec[n][2]=-1*k; n++;
	                                        posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][1]=-1*j; posGridNBtransVec[n][2]=   k; n++;
	                                        posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][1]=   j; posGridNBtransVec[n][2]=-1*k; n++;
	                                        posGridNBtransVec[n][0]=   i; posGridNBtransVec[n][1]=   j; posGridNBtransVec[n][2]=   k; n++;
					}
				}
			}
		}
		posGridNBtransN[0]=n;
	}
	return 0;
}

int assignGrpMolsToPosGrid(t_grp grp,t_mol *mols,t_vecR oriPosGrid,t_vecR dPosGrid,t_box box,
	int *nPosGrid,t_list *posGrid,int *molPosInPosGrid,int dimPosGrid) {
	int i,j,molIdx;
	int pos[3];

	for(i=0;i<grp.nMols;i++) {
	        molIdx=grp.mols[i];
		if(gridPos(mols[molIdx].COM,oriPosGrid,dPosGrid.x,dPosGrid.y,dPosGrid.z,pos,box,nPosGrid)==1) {
			j=pos[2]+dimPosGrid*(pos[1]+dimPosGrid*pos[0]);
			posGrid[j].idx[posGrid[j].n]=molIdx;
			posGrid[j].n++;
			molPosInPosGrid[molIdx]=j; /*every molecule in the group (molIdx) knows now
	                                             in which voxel of the posGrid it is located in*/
	        } else {
			printf("%d %d %d\n",pos[0],pos[1],pos[2]);
			printf("COM %10.5e %10.5e %10.5e\n",mols[molIdx].COM.x,mols[molIdx].COM.y,mols[molIdx].COM.z);
	                fatal("molecule found outside posGrid: consider using wrap=2\n");
	        }
	}
	return 0;
}

int makeNBLists(t_list *posGridX,t_list *posGridS,t_list *posGridNBListX,t_list *posGridNBListS,
        int **posGridPos,int dimPosGrid,int dimPosGridCube,int *posGridRange,int **posGridNBtransVec,int posGridNBtransN) {
        int i,i2,j,k,l,j2,k2,l2,i3,m;
        #pragma omp parallel for private(j,k,l,i2,j2,k2,l2,i3,m)
        for(i=0;i<dimPosGridCube;i++) {
                j=posGridPos[i][0];
                k=posGridPos[i][1];
                l=posGridPos[i][2];
		for(i2=0;i2<posGridNBtransN;i2++) {
			j2=j-posGridNBtransVec[i2][0];
			k2=k-posGridNBtransVec[i2][1];
			l2=l-posGridNBtransVec[i2][2];
			if(j2<0) j2+=dimPosGrid;
			if(j2>=dimPosGrid) j2-=dimPosGrid;
			if(k2<0) k2+=dimPosGrid;
			if(k2>=dimPosGrid) k2-=dimPosGrid;
			if(l2<0) l2+=dimPosGrid;
			if(l2>=dimPosGrid) l2-=dimPosGrid;
			i3=l2+dimPosGrid*(k2+dimPosGrid*j2);
			for(m=0;m<posGridX[i3].n;m++) {
			        posGridNBListX[i].idx[posGridNBListX[i].n]=posGridX[i3].idx[m];
			        posGridNBListX[i].n++;
			}
			for(m=0;m<posGridS[i3].n;m++) {
			        posGridNBListS[i].idx[posGridNBListS[i].n]=posGridS[i3].idx[m];
			        posGridNBListS[i].n++;
			}
			/*posGridNBListX/S[i] now contains a list of molecules that are possibly within
			  the cutoff distance of molecules in posGridX/S[i]*/
                }
        }
        return 0;
}

int computePot(t_atom *atoms,t_mol *mols,t_box box,int molIdx1,int molIdx2,real cutOffListSq,real cutOffPotSq,
	int LJ,int ES,int PBC,double *U,double **sigSqMat,double **epsMat,int *types,int nAtomsPerResMax) {
	real rSq,r;
	int i,l,m;
	int inRange;
	int atomIdx1,atomIdx2;
	double sigmaSq,eps,sigDivR6;
	/*double sigma;*/
	/*conversion constant convES*(q[e])^2/r[A]) to kJ/mol*/
        double convES=1389.33;
	real **distSqMat;

	distSqMat=(real**)save_malloc(nAtomsPerResMax*sizeof(real*));
        for(i=0;i<nAtomsPerResMax;i++) distSqMat[i]=(real*)save_malloc(nAtomsPerResMax*sizeof(real));

	if(molIdx1!=molIdx2) {
		inRange=0;
		for(l=0;l<mols[molIdx1].nAtoms;l++) {
			atomIdx1=mols[molIdx1].atoms[l];
			for(m=0;m<mols[molIdx2].nAtoms;m++) {
				atomIdx2=mols[molIdx2].atoms[m];
				if(PBC==1) distSqPBC(atoms[atomIdx1].crd,atoms[atomIdx2].crd,box,&rSq);
				else distRSq(atoms[atomIdx1].crd,atoms[atomIdx2].crd,&rSq);
				distSqMat[l][m]=rSq;
				if(rSq<cutOffListSq) {
					if(inRange==0) inRange=1;
					if(rSq<cutOffPotSq) inRange=2;
				}
			}
		}
		if(inRange==2) {
			for(l=0;l<mols[molIdx1].nAtoms;l++) {
				atomIdx1=mols[molIdx1].atoms[l];
				for(m=0;m<mols[molIdx2].nAtoms;m++) {
					atomIdx2=mols[molIdx2].atoms[m];
					rSq=distSqMat[l][m];
					r=sqrt(rSq);
					/*sigma=(atoms[atomIdx1].sigmaLJ+atoms[atomIdx2].sigmaLJ)/2.0;
					sigmaSq=sigma*sigma;*/
					sigmaSq=sigSqMat[types[atomIdx1]][types[atomIdx2]];
					sigDivR6=sigmaSq/rSq;
					sigDivR6=sigDivR6*sigDivR6*sigDivR6;
					/*eps=sqrt(atoms[atomIdx1].epsLJ*atoms[atomIdx2].epsLJ);*/
					eps=epsMat[types[atomIdx1]][types[atomIdx2]];
					if(LJ==1) U[0]+=
					  4.0*eps*(sigDivR6*sigDivR6-sigDivR6);
					//if(r<4.0 && eps!=0.0) printf("sig %10.5e eps %10.5e r %10.5e sigDivR6 %10.5e\n",
					  //sigma,eps,r,sigDivR6);
					if(ES==1) U[0]+=
					  convES*atoms[atomIdx1].charge*atoms[atomIdx2].charge/r;
				}
			}
		}
	}
	for(i=0;i<nAtomsPerResMax;i++) free(distSqMat[i]);
	free(distSqMat);
	return 0;
}

int getAlignRefSpec(char *fnRef,t_grp grp,t_atom *atoms,t_box box,real **w_rls,rvec **refCrd,t_vecR *refCOM,rvec **x)
{
        int i;
        FILE *ref;
        char buffer[200];
        int nRef;
        float xr,yr,zr;
        real totMass=0.0;

        ref=fopen(fnRef,"r");
        fgets(buffer,200,ref);
        fgets(buffer,200,ref);
        sscanf(buffer,"%d",&nRef);
        if(nRef!=grp.nAtoms) {
                sprintf(buffer,"nAtoms in %s (%d) and group %s (%d) don't match!\n",fnRef,nRef,grp.name,grp.nAtoms);
                fatal(buffer);
        }
        allocReals(w_rls,grp.nAtoms);
        refCrd[0]=(rvec*)save_malloc(grp.nAtoms*sizeof(rvec));
        refCOM[0].x=0.0; refCOM[0].y=0.0; refCOM[0].z=0.0;
        for(i=0;i<nRef;i++) {
                fgets(buffer,200,ref);
                sscanf(&buffer[20],"%f %f %f",&xr,&yr,&zr);
                w_rls[0][i]=atoms[grp.atoms[i]].mass;
                refCrd[0][i][0]=10.0*xr;
                refCrd[0][i][1]=10.0*yr;
                refCrd[0][i][2]=10.0*zr;
                refCOM[0].x+=w_rls[0][i]*10.0*xr;
                refCOM[0].y+=w_rls[0][i]*10.0*yr;
                refCOM[0].z+=w_rls[0][i]*10.0*zr;
                totMass+=w_rls[0][i];
        }
        fclose(ref);
        refCOM[0].x/=totMass;
        refCOM[0].y/=totMass;
        refCOM[0].z/=totMass;

        for(i=0;i<nRef;i++)
        {
                refCrd[0][i][0]-=refCOM[0].x;
                refCrd[0][i][1]-=refCOM[0].y;
                refCrd[0][i][2]-=refCOM[0].z;
        }
        x[0]=(rvec*)save_malloc(nRef*sizeof(rvec));

        return 0;
}

int printTitle() {
        printf("THE PURPOSE OF THIS PROGRAM IS TO COMPUTE\n");
	printf("LOCAL POTENTIAL ENERGIES OF SOLVENT MOLECULES\n");
	printf("IN THE ENVIRONMENT OF A SOLUTE\n");
	printf("INTERACTIONS ARE COMPUTED BETWEEN MOLECULES\n");
	printf("IN THE SOLVENT GROUP\n");
	printf("AND THE SOLUTE GROUP\n");
	printf("AND THE NON-SOLUTE GROUP\n");
	printf("THE CORRESPONDING INTERACTION ENERGIES PER\n");
	printf("SOLVENT MOLECULE ARE COMPUTED ON A GRID AND\n");
	printf("REPORTED IN A CUBE FILE\n");
        printf("\n");
        printf("Version 1.1: February 6, 2016\n");
        printf("Author:\n");
        printf(" Dr. Matthias Heyden\n");
        printf(" Max-Planck-Institut fuer Kohlenforschung\n");
        printf(" Kaiser-Wilhelm-Platz 1\n");
        printf(" D-45470 Muelheim an der Ruhr\n");
        printf(" Germany\n");
        printf(" e-mail: heyden@kofo.mpg.de\n");
        printf("\n");
        return 0;
}

int getLineFromCOM(FILE *in,char *buffer,int max)
{
        if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        while(strncmp(buffer,"#",1)==0)
        {
                if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        }
        return 0;
}

int printKeys() {
        printf("fnTop\nfnCrd\nfnVel (if format xyz,crd,dcd)\nfnJob\n");
        printf("nRead\nnSample\nfnRef\n");
        printf("alignGrp\nsoluteGrp\nsolventGrp\nnonSoluteGrp\nrotAlign\nwrap\n");
	printf("ngx\nngy\nngz\ndgridx\ndgridy\ndgridz\ncutOffPot\ncutOffList\nLJ\nES\nPBC\nstemName\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,char *fnCrd,char *fnVel,char *fnJob,
                int *nRead,int *nSample,
                char *fnRef,int *alignGrp,int *soluteGrp,int *solventGrp,
		int *nonSoluteGrp,int *rotAlign,int *wrap,int *ng,t_vecR *dgrid,real *cutOffPotSq,
		real *cutOffList,int *LJ,int *ES,int *PBC,char *fnOut,
                int needVelo) {
        FILE *io;
        char buffer[300];
        int format;
        int i=1;
        float tmp;

        saveOpenRead(&io,fnCOM);

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnTop);
        printf("%2d -> read %20s : %s\n",i,"fnTop",fnTop);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnCrd);
        getFormat(fnCrd,&format);
        printf("%2d -> read %20s : %s  ",i,"fnCrd",fnCrd);i++;
        if(format==1) printf("(xyz)\n");
        else if(format==2) printf("(crd)\n");
        else if(format==3) printf("(gro)\n");
        else if(format==4) printf("(trr)\n");
        else if(format==5) printf("(dcd)\n");
        if(needVelo==1 && (format==1 || format==2 || format==5)) {
                fgets(buffer,300,io);
                sscanf(buffer,"%s",fnVel);
                printf("%2d -> read %20s : %s\n",i,"fnVel",fnVel);i++;
        }

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnJob);
        printf("%2d -> read %20s : %s\n",i,"fnJob",fnJob);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nRead);
        printf("%2d -> read %20s : %d\n",i,"nRead",nRead[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nSample);
        printf("%2d -> read %20s : %d\n",i,"nSample",nSample[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnRef);
        printf("%2d -> read %20s : %s\n",i,"fnRef",fnRef);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",alignGrp);
        printf("%2d -> read %20s : %d\n",i,"alignGrp",alignGrp[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",soluteGrp);
        printf("%2d -> read %20s : %d\n",i,"soluteGrp",soluteGrp[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",solventGrp);
        printf("%2d -> read %20s : %d\n",i,"solventGrp",solventGrp[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nonSoluteGrp);
        printf("%2d -> read %20s : %d\n",i,"nonSoluteGrp",nonSoluteGrp[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",rotAlign);
        printf("%2d -> read %20s : %d\n",i,"rotAlign",rotAlign[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",wrap);
        printf("%2d -> read %20s : %d\n",i,"wrap",wrap[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",&ng[0]);
        printf("%2d -> read %20s : %d\n",i,"ngx",ng[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",&ng[1]);
        printf("%2d -> read %20s : %d\n",i,"ngy",ng[1]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",&ng[2]);
        printf("%2d -> read %20s : %d\n",i,"ngz",ng[2]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
	printf("%2d -> read %20s : %f\n",i,"dGridX",tmp);i++;
	dgrid[0].x=tmp;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"dGridY",tmp);i++;
        dgrid[0].y=tmp;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"dGridZ",tmp);i++;
        dgrid[0].z=tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"cutOffPot",tmp);i++;
        cutOffPotSq[0]=tmp*tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"cutOffList",tmp);i++;
        cutOffList[0]=tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",LJ);
        printf("%2d -> read %20s : %d\n",i,"LJ",LJ[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",ES);
        printf("%2d -> read %20s : %d\n",i,"ES",ES[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",PBC);
        printf("%2d -> read %20s : %d\n",i,"PBC",PBC[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%2d -> read %20s : %s\n",i,"fnOut",fnOut);i++;

        fclose(io);
        return 0;
}

int readTRJbuffer(int *pos,real *time,real *timebuffer,t_box *box,t_box *boxbuffer,int nAtoms,t_atom *atoms,t_atom **atomsbuffer,int nSites,t_vecR *siteCrds,t_vecR **siteCrdsbuffer)
{
        int i,j;
        i=pos[0];
        time[0]=timebuffer[i];
        vecRCpy(boxbuffer[i].a,&box[0].a);
        vecRCpy(boxbuffer[i].b,&box[0].b);
        vecRCpy(boxbuffer[i].c,&box[0].c);
        box[0].ortho=boxbuffer[i].ortho;
        vecRCpy(boxbuffer[i].trans.c1,&box[0].trans.c1);
        vecRCpy(boxbuffer[i].trans.c2,&box[0].trans.c2);
        vecRCpy(boxbuffer[i].trans.c3,&box[0].trans.c3);
        for(j=0;j<nAtoms;j++) {
                vecRCpy(atomsbuffer[i][j].crd,&atoms[j].crd);
                vecRCpy(atomsbuffer[i][j].vel,&atoms[j].vel);
        }
        for(j=0;j<nSites;j++) {
                vecRCpy(siteCrdsbuffer[i][j],&siteCrds[j]);
        }
        pos[0]++;
        return 0;
}

int main(int argc, char *argv[])
{
	FILE  *tin,*cin,*vin,*job,*out;
	XDR xdrin;
	char fnCOM[100],fnTop[100],fnCrd[100],fnVel[100],fnJob[100],fnOut[100],fnRef[100];
	char fnOut2[100];
	t_atom *atoms;
	t_vecR *siteCrds;
	t_mol *mols;
	t_grp *grps;
	int nAtoms,nMols,nSites,format;
	int nFrames=0;
	real dt=1.0;
	real time=0.0;
	t_box box;
	t_grpDef grpDef;
	int needVelo=0;
	int veloFreq=1;
	int wc=0;
	int i=0;
	int j=0;
	int k,l,m,n,s;
	int nRead=5000;
	int nTimesAnalysis=1000;
	int analysisInterval=10;
	int nSample=10;
	int alignGrp;
	int rotAlign;
	int wrap;
	rvec *x;
        real *w_rls;
        rvec *refCrd;
        t_vecR refCOM;
	real *timebuffer;
        t_box *boxbuffer;
        t_atom **atomsbuffer;
        t_vecR **siteCrdsbuffer;
        int bufferpos;

	int soluteGrp,solventGrp,nonSoluteGrp;
	int ng[3];
	t_vecR dgrid;
	real cutOffPotSq;
	real cutOffList,cutOffListSq;
	t_vecR gridOri;
	t_voxelProp *vox;
	t_list *voxMols;
	int ng3,voxIdx;
	int pos[3];
	int cnt;
	real ***volData;
	t_box voxel;
	char volTitle[300];
	char fnVol[300];
	real tmp;
	int molIdx1,molIdx2;
	int LJ,ES,PBC;
	int nAtomsPerResMax;
	double **sigSqMat;
	double **epsMat;
	int *types;
	/*this grid is aligned with the coordinate axes*/
	/*it fits exactly into the (orthorhombic) periodic box*/
	/*its purpose is to speed up the generation of neighbour lists*/
	int dimPosGrid=13;
	int dimPosGridCube;
	t_list *posGridX;
	t_list *posGridS;
	int **posGridPos; /*translation into [i][j][k] of the 1D array position of posGridX/S and posGridNBListX/S*/
	t_vecR dPosGrid;
	t_vecR oriPosGrid;
	int nPosGrid[3];
	t_list *posGridNBListX;
	t_list *posGridNBListS;
	int *molPosInPosGrid;
	int posGridRange[3];
	int **posGridNBtransVec;
	int posGridNBtransN;
	int i2,j2,k2,l2,i3,j3,k3,l3;

	printTitle();

        timebuffer=(real*)save_malloc(100*sizeof(real));
        boxbuffer=(t_box*)save_malloc(100*sizeof(t_box));
        atomsbuffer=(t_atom**)save_malloc(100*sizeof(t_atom*));
        siteCrdsbuffer=(t_vecR**)save_malloc(100*sizeof(t_vecR*));

/*PARTIALLY MODULAR*/
	/*PARSE COMMAND LINE INPUT HERE*/
	if(argc!=2) {
                printKeys();
                fatal("no input file specified\n");
        }
        strcpy(fnCOM,argv[1]);
	getInput(fnCOM,fnTop,fnCrd,fnVel,fnJob,
                &nRead,&nSample,fnRef,&alignGrp,&soluteGrp,&solventGrp,
                &nonSoluteGrp,&rotAlign,&wrap,ng,&dgrid,&cutOffPotSq,
                &cutOffList,&LJ,&ES,&PBC,fnOut,
                needVelo);
	cutOffListSq=cutOffList*cutOffList;
/*PARTIALLY MODULAR*/

	/*READ TOPOLOGY HERE*/
	readTOP(fnTop,&atoms,&nAtoms,&siteCrds,&nSites,&mols,&nMols);
	//for(i=0;i<nAtoms;i++) printf("# %d q %10.5e sigLJ %10.5e epsLJ %10.5e\n",i+1,atoms[i].charge,atoms[i].sigmaLJ,atoms[i].epsLJ);
	for(i=0;i<100;i++) {
                atomsbuffer[i]=(t_atom*)save_malloc(nAtoms*sizeof(t_atom));
                siteCrdsbuffer[i]=(t_vecR*)save_malloc(nSites*sizeof(t_vecR));
        }
	/*mols[nMols] is used to store bonds, angles, dihedrals and impropers between residues*/

	/*READ JOB DATA HERE*/
	readJob(fnJob,&grpDef);

	/*SELECT STATIC GROUPS FROM TOPOLOGY*/
	allocGrps(grpDef,nAtoms,nMols,&grps);
	makeStaticGrps(atoms,nAtoms,mols,nMols,grpDef.statGrpDef,grpDef.nStat,grps);

	/*PREP CRD AND/OR VEL READ HERE*/
	prepTRJinput(fnCrd,fnVel,&cin,&vin,&xdrin,&format,&veloFreq,&time,&dt,&nFrames,needVelo,nAtoms,argc);

/*MODULAR */
	/*we figure out how many atoms the largest residue in the system has*/
	j=-1;
	for(i=0;i<nMols;i++) {
		if(mols[i].nAtoms>j) j=mols[i].nAtoms;
	}
	nAtomsPerResMax=j;

	/*this will be a list of indices indicating where to find the LJ parameters for each pair of atoms in the comination parameter matrix*/
        types=(int*)save_malloc(nAtoms*sizeof(int));
	makeLJcombParameters(atoms,nAtoms,types,&sigSqMat,&epsMat);

	/*here we generate a grid of molecule index lists to speed up neighbor list generation for large systems*/
	dimPosGridCube=dimPosGrid*dimPosGrid*dimPosGrid;
	allocPosGridStuff(&posGridX,&posGridS,&posGridNBListX,&posGridNBListS,&posGridPos,dimPosGrid,dimPosGridCube,nMols);
	nPosGrid[0]=dimPosGrid; nPosGrid[1]=dimPosGrid; nPosGrid[2]=dimPosGrid;
	molPosInPosGrid=(int*)save_malloc(nMols*sizeof(int));
	posGridNBtransVec=(int**)save_malloc(dimPosGridCube*sizeof(int*));
	for(i=0;i<dimPosGridCube;i++) posGridNBtransVec[i]=(int*)save_malloc(3*sizeof(int));
	oriPosGrid.x=0.0; oriPosGrid.y=0.0; oriPosGrid.z=0.0;
	dPosGrid.x=0.0; dPosGrid.y=0.0; dPosGrid.z=0.0;

	if(alignGrp!=-1) getAlignRefSpec(fnRef,grps[alignGrp],atoms,box,&w_rls,&refCrd,&refCOM,&x);
	if(soluteGrp>=grpDef.nStat) fatal("solute group needs to be static\n");
	if(solventGrp>=grpDef.nStat) fatal("solvent group needs to be static\n");
	if(nonSoluteGrp>=grpDef.nStat) fatal("non-solute group needs to be static\n");
	gridOrigin(ng,dgrid,&gridOri);
	ng3=ng[0]*ng[1]*ng[2];
	vox=(t_voxelProp*)save_malloc(ng3*sizeof(t_voxelProp));
	voxMols=(t_list*)save_malloc(ng3*sizeof(t_list));
	for(i=0;i<ng3;i++) {
		initVoxelProp(&vox[i]);
		voxMols[i].idx=(int*)save_malloc(100*sizeof(int));
	}

	cnt=0;
	volData=(real***)save_malloc(ng[0]*sizeof(real**));
	for(i=0;i<ng[0];i++) {
		volData[i]=(real**)save_malloc(ng[1]*sizeof(real*));
		for(j=0;j<ng[1];j++) {
			volData[i][j]=(real*)save_malloc(ng[2]*sizeof(real));
		}
	}
/*MODULAR*/

	/*HAPPY COMPUTING!!!*/
	for(i=0;i<nRead;i++)
	{
		if(i%100==0) {
                        for(j=0;j<100;j++) {
                                readTRJ(i+j,format,cin,vin,&xdrin,dt,&timebuffer[j],&boxbuffer[j],nAtoms,nSites,needVelo,veloFreq,wc,atomsbuffer[j],siteCrdsbuffer[j]);
                        }
                        bufferpos=0;
                }
                readTRJbuffer(&bufferpos,&time,timebuffer,&box,boxbuffer,nAtoms,atoms,atomsbuffer,nSites,siteCrds,siteCrdsbuffer);

		//if(i==0 && alignGrp!=-1) getAlignRef(&grps[alignGrp],atoms,box,&w_rls,&refCrd,&refCOM,&x);
		if(i%nSample==0) {
                	//if(alignGrp!=-1) alignGroup(&grps[alignGrp],atoms,nAtoms,&box,w_rls,refCrd,refCOM,x);
			if(alignGrp!=-1) alignGroupNoRot(&grps[alignGrp],atoms,nAtoms,box,refCOM);
			if(wrap==1 || wrap==2) wrapTrajectoryMol0(atoms,mols,nMols,box);
			getMolCOM(atoms,mols,nMols,needVelo);
			if(wrap==2) wrapTrajectoryMolCOM(atoms,mols,nMols,box);
			/*getMolDM(atoms,mols,nMols,needVelo);*/
			prepPosGridParameters(posGridX,posGridS,posGridNBListX,posGridNBListS,
				&dPosGrid,&oriPosGrid,posGridRange,dimPosGrid,dimPosGridCube,box,
				cutOffList,cutOffListSq,posGridNBtransVec,&posGridNBtransN);
			assignGrpMolsToPosGrid(grps[soluteGrp],mols,oriPosGrid,dPosGrid,box,nPosGrid,posGridX,molPosInPosGrid,dimPosGrid);
			assignGrpMolsToPosGrid(grps[solventGrp],mols,oriPosGrid,dPosGrid,box,nPosGrid,posGridS,molPosInPosGrid,dimPosGrid);
			makeNBLists(posGridX,posGridS,posGridNBListX,posGridNBListS,posGridPos,dimPosGrid,dimPosGridCube,posGridRange,
				posGridNBtransVec,posGridNBtransN);
			//for(j=0;j<dimPosGrid;j++) { for(k=0;k<dimPosGrid;k++) { for(l=0;l<dimPosGrid;l++) {
			//	m=l+dimPosGrid*(k+dimPosGrid*j);
			//	printf("%6d",posGridNBListS[m].n);
			//} printf("\n"); }}
			if(rotAlign==1) {
                                alignGroup(&grps[alignGrp],atoms,nAtoms,&box,w_rls,refCrd,refCOM,x);
                                getMolCOM(atoms,mols,nMols,needVelo);
                        }
		}
		
		if(i%nSample==0)
                {
                        /*updateGrps(atoms,nAtoms,mols,nMols,box,grpDef,grps);*/
/*MODULAR*/
			#pragma omp parallel for
			for(j=0;j<ng3;j++) {
				voxMols[j].n=0;
			}
			for(j=0;j<grps[solventGrp].nMols;j++) {
				molIdx1=grps[solventGrp].mols[j];
				if(gridPos(mols[molIdx1].COM,gridOri,dgrid.x,dgrid.y,dgrid.z,pos,box,ng)==1) {
                                        voxIdx=pos[2]+ng[2]*(pos[1]+ng[1]*pos[0]);
					voxMols[voxIdx].idx[voxMols[voxIdx].n]=molIdx1;
					voxMols[voxIdx].n++;
					vox[voxIdx].nMol++;
				}
			}
			#pragma omp parallel for private(k,molIdx1,l,m,molIdx2)
			for(j=0;j<ng3;j++) {
				for(k=0;k<voxMols[j].n;k++) {
					molIdx1=voxMols[j].idx[k];
					l=molPosInPosGrid[molIdx1];
					for(m=0;m<posGridNBListX[l].n;m++) {
						molIdx2=posGridNBListX[l].idx[m];
						computePot(atoms,mols,box,molIdx1,molIdx2,cutOffListSq,
						cutOffPotSq,LJ,ES,PBC,&vox[j].Uxs,
						sigSqMat,epsMat,types,nAtomsPerResMax);
					}
					for(m=0;m<posGridNBListS[l].n;m++) {
                                                molIdx2=posGridNBListS[l].idx[m];
                                                computePot(atoms,mols,box,molIdx1,molIdx2,cutOffListSq,
                                                cutOffPotSq,LJ,ES,PBC,&vox[j].Uss,
                                                sigSqMat,epsMat,types,nAtomsPerResMax);
                                        }
				}
			}
		        cnt++;
/*MODULAR*/
                }
		printf("step %d of %d%c",i+1,nRead,(char)13); fflush(stdout);
	}
	printf("\n");
	closeInput(format,needVelo,cin,vin,&xdrin);

	for(i=0;i<ng3;i++) {
		vox[i].averMol=((double)vox[i].nMol)/((double)cnt);
		if(vox[i].nMol!=0) {
			vox[i].Uxs/=vox[i].nMol;
			vox[i].Uss/=vox[i].nMol;
		}
	}

	/******************************************************************/
	/* OUTPUT OF CUBE FILES FOR AVERAGE POTENTIAL ENERGY DATA         */
	/******************************************************************/
	voxel.a.x=dgrid.x; voxel.a.y=0.0; voxel.a.z=0.0;
	voxel.b.x=0.0; voxel.b.y=dgrid.y; voxel.b.z=0.0;
	voxel.c.x=0.0; voxel.c.y=0.0; voxel.c.z=dgrid.z;

	sprintf(volTitle,"CUBE file of average solvent molecule number density (nm^-3)");
        sprintf(fnVol,"%s_pot3D-numDens.cube",fnOut);
	l=0;
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[l].averMol/(dgrid.x*dgrid.y*dgrid.z*0.001);
		l++;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average solute-solvent interaction energy per solvent molecule (kJ/mol)");
        sprintf(fnVol,"%s_pot3D-Uxs.cube",fnOut);
	l=0;
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[l].Uxs;
		l++;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average solvens-solvent interaction energy per solvent molecule (kJ/mol)");
        sprintf(fnVol,"%s_pot3D-Uss.cube",fnOut);
	l=0;
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[l].Uss;
		l++;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	return 0;
}

