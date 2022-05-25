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
        t_vecR  COM;
        t_vecR  axes[3];
	t_vecR  bDir[2];
        t_vecR  vCOM;
        t_vecR  vel[3];
        t_vecR  angMom;
        real	angMomProj[3];
} t_molProp;

typedef struct {
        int	nMol;
	real	averMol;
        real    *COMmsd;
        real    *rotCorr[3];
        real    *vCOMcorr;
        real    *angMomCorr;
        real    *angMomProjCorr[3];
	real	*atVelCorr[3];
	t_vecR	toProt,polV;
	real	polToProt,polNorm,OHtoProt,molPlaneToProt;
} t_voxelProp;

int initVoxelProp(t_voxelProp *p,int nDelay) {
        int i,j;
        p[0].nMol=0;
        p[0].COMmsd=(real*)save_malloc(nDelay*sizeof(real));
        p[0].vCOMcorr=(real*)save_malloc(nDelay*sizeof(real));
        p[0].angMomCorr=(real*)save_malloc(nDelay*sizeof(real));
	for(i=0;i<nDelay;i++) {
		p[0].COMmsd[i]=0.0;
		p[0].vCOMcorr[i]=0.0;
		p[0].angMomCorr[i]=0.0;
	}
        for(i=0;i<3;i++) {
                p[0].rotCorr[i]=(real*)save_malloc(nDelay*sizeof(real));
                p[0].angMomProjCorr[i]=(real*)save_malloc(nDelay*sizeof(real));
		p[0].atVelCorr[i]=(real*)save_malloc(nDelay*sizeof(real));
		for(j=0;j<nDelay;j++) {
			p[0].rotCorr[i][j]=0.0;
			p[0].angMomProjCorr[i][j]=0.0;
			p[0].atVelCorr[i][j]=0.0;
		}
        }
	p[0].polV.x=0.0; p[0].polV.y=0.0; p[0].polV.z=0.0;
	p[0].polNorm=0.0;
        p[0].polToProt=0.0;
        p[0].OHtoProt=0.0;
        p[0].molPlaneToProt=0.0;
        return 0;
}

int getDirToProt(t_voxelProp ***p,rvec *refCrd,int nProt,int *ng,t_vecR gridOri,real dgrid) {
	int i,j,k,l;
	t_vecR pt;
	real min,cur,tmp;
	int idx;

	for(i=0;i<ng[0];i++) {
		for(j=0;j<ng[1];j++) {
			for(k=0;k<ng[2];k++) {
				pt.x=gridOri.x+(i+0.5)*dgrid;
				pt.y=gridOri.y+(j+0.5)*dgrid;
				pt.z=gridOri.z+(k+0.5)*dgrid;
				idx=0;
				min=999999.9;
				for(l=0;l<nProt;l++) {
					tmp=pt.x-refCrd[l][0];
					cur=tmp*tmp;
					tmp=pt.y-refCrd[l][1];
					cur+=tmp*tmp;
					tmp=pt.z-refCrd[l][2];
					cur+=tmp*tmp;
					if(cur<min) {
						min=cur;
						idx=l;
					}
				}
				p[i][j][k].toProt.x=refCrd[l][0]-pt.x;
				p[i][j][k].toProt.y=refCrd[l][1]-pt.y;
				p[i][j][k].toProt.z=refCrd[l][2]-pt.z;
				vecRNorm(p[i][j][k].toProt,&tmp);
				p[i][j][k].toProt.x/=tmp;
				p[i][j][k].toProt.y/=tmp;
				p[i][j][k].toProt.z/=tmp;
			}
		}
	}
	return 0;
}


int gridOrigin(int *ng,real dgrid,t_vecR *ori) {
	ori[0].x=(((real)ng[0])/-2.0)*dgrid;
	ori[0].y=(((real)ng[1])/-2.0)*dgrid;
	ori[0].z=(((real)ng[2])/-2.0)*dgrid;
	return 0;
}

int gridPos(t_vecR crd,t_vecR ori,real dgrid,int *pos,t_box box,int *ng) {
	t_vecR link;
	int check;

	/*grid ori is the origin of the grid (one of tis corners, not the center)*/
	//linkPBC(ori,crd,box,&link);
	vecRSub(ori,crd,&link);
	if(link.x<0.0 || link.y<0.0 || link.z<0.0) {
		check=0; /*outside of grid volume*/
	} else {
		pos[0]=(int)floor(link.x/dgrid);
		if(pos[0]>=ng[0]) {
			check=0; /*outside of grid volume*/
		} else {
			pos[1]=(int)floor(link.y/dgrid);
			if(pos[1]>=ng[1]) {
				check=0; /*outside of grid volume*/
			} else {
				pos[2]=(int)floor(link.z/dgrid);
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

int checkWater(t_atom *atoms,t_mol mol) {
	int i,n;
	int nO=0;
	int nH=0;
	for(i=0;i<mol.nAtoms;i++) {
		n=mol.atoms[i];
		if(strncmp(atoms[n].atomName,"O",1)==0) nO++;
		if(strncmp(atoms[n].atomName,"H",1)==0) nH++;
		if(nO==0 && nH!=0) { 
			for(i=0;i<mol.nAtoms;i++) {
				n=mol.atoms[i];
				printf("%s\n",atoms[n].atomName);
			}
			fatal("expect O H H atom order in water molecules\n");
		}
	}
	if(nO!=1 || nH!=2) fatal("found non water molecule\n");
	return 0;
}

int checkWaterGroup(t_atom *atoms,t_mol *mols,t_grp grp,int *OHH) {
	int i,n,j;
	for(i=0;i<grp.nMols;i++) {
		checkWater(atoms,mols[grp.mols[i]]);
	}
	if(grp.nMols==0) fatal("no molecules in analysis group\n");
	j=0;
	for(i=0;i<mols[grp.mols[0]].nAtoms;i++) {
		n=mols[grp.mols[0]].atoms[i];
		if(strncmp(atoms[n].atomName,"O",1)==0) OHH[0]=i;
		if(strncmp(atoms[n].atomName,"H",1)==0) {
			if(j==0) {
				OHH[1]=i;
				j=1;
			} else {
				OHH[2]=i;
			}
		}
	}
	return 0;
}

int getMolProp(t_atom *atoms,t_mol mol,t_molProp *p,int *OHH) {
	int i,n;
	real totM;
	t_vecR b1,b2,tmp1,tmp3,tmp4;
	real tmp2;

	p->COM.x=0.0; p->vCOM.x=0.0;
	p->COM.y=0.0; p->vCOM.y=0.0;
	p->COM.z=0.0; p->vCOM.z=0.0;
	totM=mol.mass;
	for(i=0;i<3;i++) {
		n=mol.atoms[OHH[i]];
		p->COM.x+=atoms[n].mass*atoms[n].crd.x;
		p->COM.y+=atoms[n].mass*atoms[n].crd.y;
		p->COM.z+=atoms[n].mass*atoms[n].crd.z;
		p->vCOM.x+=atoms[n].mass*atoms[n].vel.x;
                p->vCOM.y+=atoms[n].mass*atoms[n].vel.y;
                p->vCOM.z+=atoms[n].mass*atoms[n].vel.z;
		vecRCpy(atoms[n].vel,&p[0].vel[i]);
	}
	p->COM.x/=totM; p->vCOM.x/=totM;
	p->COM.y/=totM; p->vCOM.y/=totM;
	p->COM.z/=totM; p->vCOM.z/=totM;

	vecRSub(atoms[mol.atoms[OHH[0]]].crd,atoms[mol.atoms[OHH[1]]].crd,&b1);
	vecRSub(atoms[mol.atoms[OHH[0]]].crd,atoms[mol.atoms[OHH[2]]].crd,&b2);
	vecRNorm(b1,&tmp2);
	p->bDir[0].x=b1.x/tmp2;
	p->bDir[0].y=b1.y/tmp2;
	p->bDir[0].z=b1.z/tmp2;
	vecRNorm(b2,&tmp2);
        p->bDir[1].x=b2.x/tmp2;
        p->bDir[1].y=b2.y/tmp2;
        p->bDir[1].z=b2.z/tmp2;

	vecRAdd(b1,b2,&tmp1);
	vecRNorm(tmp1,&tmp2);
	p->axes[2].x=-1.0*tmp1.x/tmp2;
	p->axes[2].y=-1.0*tmp1.y/tmp2;
	p->axes[2].z=-1.0*tmp1.z/tmp2;
	vecRCross(b1,b2,&tmp1);
	vecRNorm(tmp1,&tmp2);
        p->axes[1].x=tmp1.x/tmp2;
        p->axes[1].y=tmp1.y/tmp2;
        p->axes[1].z=tmp1.z/tmp2;
	vecRCross(p->axes[1],p->axes[2],&tmp1);
	vecRNorm(tmp1,&tmp2);
        p->axes[0].x=tmp1.x/tmp2;
        p->axes[0].y=tmp1.y/tmp2;
        p->axes[0].z=tmp1.z/tmp2;

	p->angMom.x=0.0;
	p->angMom.y=0.0;
	p->angMom.z=0.0;
	for(i=0;i<3;i++) {
		n=mol.atoms[OHH[i]];
		vecRSub(p->COM,atoms[n].crd,&tmp1);
		vecRSub(p->vCOM,atoms[n].vel,&tmp3);
		vecRScal(atoms[n].mass,tmp3,&tmp4);
		vecRCross(tmp1,tmp4,&tmp3);
		p->angMom.x+=tmp3.x;
		p->angMom.y+=tmp3.y;
		p->angMom.z+=tmp3.z;
	}
	for(i=0;i<3;i++) {
		vecRProd(p->axes[i],p->angMom,&p[0].angMomProj[i]);
	}
	return 0;
}

int addCorr(t_molProp *molProp,t_voxelProp *vox,int idxCnt,int nDelay,t_box box) {
	int i,j,k;
	real tmp1,tmp2;

	vox[0].nMol++;
	#pragma omp parallel for private(j,k,tmp1)
	for(i=0;i<nDelay;i++) {
		j=i+idxCnt;
		if(j>=nDelay) j-=nDelay;
		//distSqPBC(molProp[idxCnt].COM,molProp[j].COM,box,&tmp1);
		distRSq(molProp[idxCnt].COM,molProp[j].COM,&tmp1);
		vox[0].COMmsd[i]+=tmp1;
		vecRProd(molProp[idxCnt].vCOM,molProp[j].vCOM,&tmp1);
		vox[0].vCOMcorr[i]+=tmp1;
		vecRProd(molProp[idxCnt].angMom,molProp[j].angMom,&tmp1);
                vox[0].angMomCorr[i]+=tmp1;
		for(k=0;k<3;k++) {
			vecRProd(molProp[idxCnt].axes[k],molProp[j].axes[k],&tmp1);
			vox[0].rotCorr[k][i]+=tmp1;
			vecRProd(molProp[idxCnt].vel[k],molProp[j].vel[k],&tmp1);
			vox[0].atVelCorr[k][i]+=tmp1;
			tmp1=molProp[idxCnt].angMomProj[k]*molProp[j].angMomProj[k];
			vox[0].angMomProjCorr[k][i]+=tmp1;
		}
	}

	/*here we compute dipole moment and eventually the average dipole moment) */
	/*vector of water (length normalized to 1) molecules in this voxel, to */
	/*see if there is an average polarization. axes[2] of each molecule is the */
	/*negative of the dipole moment vector (assuming a rigid non-polarizable */
	/*water molecule, therefore the factor -1.0*/
	vox[0].polV.x+=-1.0*molProp[idxCnt].axes[2].x;
	vox[0].polV.y+=-1.0*molProp[idxCnt].axes[2].y;
	vox[0].polV.z+=-1.0*molProp[idxCnt].axes[2].z;

	/*here we compute the projection of the (normalized) dipole moment on a */
	/*vector pointing to the protein surface*/
	/*the '-=' is used again, because axes[2] is the negative dipole moment direction*/
	vecRProd(molProp[idxCnt].axes[2],vox[0].toProt,&tmp1);
	vox[0].polToProt-=tmp1;

	vecRProd(molProp[idxCnt].bDir[0],vox[0].toProt,&tmp1);
	vecRProd(molProp[idxCnt].bDir[1],vox[0].toProt,&tmp2);
	/*select here the OH bond vector that points more to the protein*/
	if(tmp1>tmp2) {
		vox[0].OHtoProt+=tmp1;
	} else {
		vox[0].OHtoProt+=tmp2;
	}
	/*here chose the normal vector of the water molecular plane that points to */
	/*the protein surface*/
	vecRProd(molProp[idxCnt].axes[0],vox[0].toProt,&tmp1);
	if(tmp1>0) {
		vox[0].molPlaneToProt+=tmp1;
	} else {
		vox[0].molPlaneToProt-=tmp1;
	}
	return 0;
}

/*array is an fftwf_complex array with 2*n-1 elements			*/
/*  => to be allocated in main						*/
/*ft is the Fourier transform plan for 2*n-1 data points		*/
/*  => prepare plan in main						*/
/*data is the real input data, to be replaced by result			*/
/*n is the number of elements in data					*/
/*norm is what the FT is going to be divided by, e.g. sqrt(2*n-1)	*/
int symFT(fftwf_plan ft,fftwf_complex *array,real *data,int n,real norm) {
	int i,j;
	for(i=0;i<n;i++) {
		array[i][0]=data[i];
		array[i][1]=0.0;
		if(i!=0) {
			j=2*n-1-i;
			array[j][0]=data[i];
			array[j][1]=0.0;
		}
	}
	fftwf_execute(ft);
	for(i=0;i<n;i++) {
		data[i]=array[i][0]/norm;
	}
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
        printf("THE PURPOSE OF THIS PROGRAM IS TO ANALYZE WATER DYNAMICS\n");
	printf("1. COM MSD\n2. ROT. TCF (+FT) FOR MOL. AXES\n3. ATOMIC VELOCITY TCF\n");
	printf("4. ATOMIC VDOS\n5. COM VELOCITY TCF\n6. COM VDOS\n7. ANG. MOM. TCF\n");
	printf("8. ANG. MOM. VDOS\n9. ANG. MOM. PROJ. ON MOL. AXES TCF\n");
	printf("10.ANG. MOM. PROJ. ON MOL. AXES VDOS\n");
	printf("ON A GRID IN THE 3-DIMENSIONAL ENVIRONMENT OF A BIOMOLECULE\n");
	printf("WITHIN A SPECIFIED TIME WINDOW (nDelay) AND TO WRITE OUT\n");
	printf("THE RESULTS IN A TABULATED DATA FORMAT AS WELL AS A CUBE FILE\n");
        printf("\n");
        printf("Version 2.0: July, 2019\n");
        printf("Author:\n");
        printf(" Dr. Matthias Heyden\n");
        printf(" Arizona State University\n");
        printf(" School of Molecular Sciences\n");
        printf(" Tempe, AZ 85287\n");
        printf(" USA\n");
        printf(" e-mail: mheyden1@asu.edu\n");
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
        printf("alignGrp\nanalyzeGrp\nwrap\nnDelay\nngx\nngy\nngz\ndgrid\nrad\nstemName\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,char *fnCrd,char *fnVel,char *fnJob,
                int *nRead,int *nSample,
                char *fnRef,int *alignGrp,int *analyzeGrp,int *wrap,
		int *nDelay,int *ng,real *dgrid,real *radSq,char *fnOut,
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
        sscanf(buffer,"%d",analyzeGrp);
        printf("%2d -> read %20s : %d\n",i,"analyzeGrp",analyzeGrp[0]);i++;

        getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",wrap);
        printf("%2d -> read %20s : %d\n",i,"wrap",wrap[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nDelay);
        printf("%2d -> read %20s : %d\n",i,"nDelay",nDelay[0]);i++;

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
	printf("%2d -> read %20s : %f\n",i,"dGrid",tmp);i++;
	dgrid[0]=tmp;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",&tmp);
        printf("%2d -> read %20s : %f\n",i,"rad",tmp);i++;
        radSq[0]=tmp*tmp;

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
	int needVelo=1;
	int veloFreq=1;
	int wc=0;
	int i=0;
	int j=0;
	int k,l,m,n,s;
	int nRead=5000;
	int nTimesAnalysis=1000;
	int analysisInterval=10;
	int nSample=10;
	int alignGrp,analyzeGrp;
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

	int nDelay;
	real tstep,fstep;
	real delayTime;
	real *times,*freqs;
	int ng[3];
	real dgrid,radSq;
	t_vecR gridOri;
	t_voxelProp ***vox;
	t_molProp **molProp;
	int OHH[3];
	int pos[3];
	int idxCnt;
	int begin;
	int cnt;
	real ***volData;
	t_box voxel;
	char volTitle[300];
	char fnVol[300];
	fftwf_plan ft;
	fftwf_complex *array;
	real FTnorm;
	real tmp;

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
                &nRead,&nSample,fnRef,
                &alignGrp,&analyzeGrp,&wrap,
                &nDelay,ng,&dgrid,&radSq,fnOut,
                needVelo);
/*PARTIALLY MODULAR*/

	/*READ TOPOLOGY HERE*/
	readTOP(fnTop,&atoms,&nAtoms,&siteCrds,&nSites,&mols,&nMols);
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
	getAlignRefSpec(fnRef,grps[alignGrp],atoms,box,&w_rls,&refCrd,&refCOM,&x);
	if(nDelay%nSample!=0) fatal("'nDelay' needs to be multiple of 'nSample'\n");
	nDelay=nDelay/nSample;
	if(analyzeGrp>=grpDef.nStat) fatal("analysis group needs to be static\n");
	checkWaterGroup(atoms,mols,grps[analyzeGrp],OHH);
	molProp=(t_molProp**)save_malloc(grps[analyzeGrp].nMols*sizeof(t_molProp*));
	for(i=0;i<grps[analyzeGrp].nMols;i++) {
		molProp[i]=(t_molProp*)save_malloc(nDelay*sizeof(t_molProp));
	}
	gridOrigin(ng,dgrid,&gridOri);
	vox=(t_voxelProp***)save_malloc(ng[0]*sizeof(t_voxelProp**));
	for(i=0;i<ng[0];i++) {
		vox[i]=(t_voxelProp**)save_malloc(ng[1]*sizeof(t_voxelProp*));
		for(j=0;j<ng[1];j++) {
			vox[i][j]=(t_voxelProp*)save_malloc(ng[2]*sizeof(t_voxelProp));
			for(k=0;k<ng[2];k++) {
				initVoxelProp(&vox[i][j][k],nDelay);
			}
		}
	}
	getDirToProt(vox,refCrd,grps[alignGrp].nAtoms,ng,gridOri,dgrid);
	idxCnt=0;
	begin=0;
	cnt=0;
	volData=(real***)save_malloc(ng[0]*sizeof(real**));
	for(i=0;i<ng[0];i++) {
		volData[i]=(real**)save_malloc(ng[1]*sizeof(real*));
		for(j=0;j<ng[1];j++) {
			volData[i][j]=(real*)save_malloc(ng[2]*sizeof(real));
		}
	}
	array=(fftwf_complex*)fftwf_malloc((2*nDelay-1)*sizeof(fftwf_complex));
	ft=fftwf_plan_dft_1d(2*nDelay-1,array,array,FFTW_FORWARD,FFTW_ESTIMATE);
	FTnorm=sqrt(((real)(2*nDelay-1)));
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
		if(i==1) {
			tstep=nSample*(timebuffer[1]-timebuffer[0]);
			delayTime=nDelay*tstep;
			fstep=1.0/((2*nDelay-1)*tstep)*33.3564;
			times=(real*)save_malloc(nDelay*sizeof(real));
			freqs=(real*)save_malloc(nDelay*sizeof(real));
			for(j=0;j<nDelay;j++) {
				times[j]=j*tstep;
				freqs[j]=j*fstep;
			}
		}

		//if(i==0 && alignGrp!=-1) getAlignRef(&grps[alignGrp],atoms,box,&w_rls,&refCrd,&refCOM,&x);
		if(i%nSample==0) {
                	/*if(alignGrp!=-1) alignGroup(&grps[alignGrp],atoms,nAtoms,&box,w_rls,refCrd,refCOM,x);*/
			if(alignGrp!=-1) alignGroupNoRotVel(&grps[alignGrp],atoms,nAtoms,box,refCOM);
			if(wrap==1 || wrap ==2) wrapTrajectoryMol0(atoms,mols,nMols,box);
			getMolCOM(atoms,mols,nMols,needVelo);
			if(wrap==2) wrapTrajectoryMolCOM(atoms,mols,nMols,box);
			/*getMolDM(atoms,mols,nMols,needVelo);*/
		}
		
		if(i%nSample==0)
                {
                        /*updateGrps(atoms,nAtoms,mols,nMols,box,grpDef,grps);*/
/*MODULAR*/
			#pragma omp parallel for
			for(j=0;j<grps[analyzeGrp].nMols;j++) {
				getMolProp(atoms,mols[grps[analyzeGrp].mols[j]],&molProp[j][idxCnt],OHH);
			}
			idxCnt++;
			if(idxCnt==nDelay) { idxCnt=0; begin=1; }
			if(begin==1) {
				for(j=0;j<grps[analyzeGrp].nMols;j++) {
					/*this checks if the molecule is within the grid volume*/
					/*'pos' will contain the voxel index coordinates*/
					/*idxCnt points to the first element of the 'rotating' list of past nDelay data points*/
					if(gridPos(molProp[j][idxCnt].COM,gridOri,dgrid,pos,box,ng)==1) {
						addCorr(molProp[j],&vox[pos[0]][pos[1]][pos[2]],idxCnt,nDelay,box);
					}
				}
				cnt++;
			}
/*MODULAR*/
                }
		printf("step %d of %d%c",i+1,nRead,(char)13); fflush(stdout);
	}
	printf("\n");
	closeInput(format,needVelo,cin,vin,&xdrin);

	for(i=0;i<ng[0];i++) {
		for(j=0;j<ng[1];j++) {
			for(k=0;k<ng[2];k++) {
				vox[i][j][k].averMol=((double)vox[i][j][k].nMol)/((double)cnt);
				if(vox[i][j][k].nMol!=0) {
					#pragma omp parallel for private(m)
					for(l=0;l<nDelay;l++) {
						vox[i][j][k].COMmsd[l]/=vox[i][j][k].nMol;
						vox[i][j][k].vCOMcorr[l]/=vox[i][j][k].nMol;
						vox[i][j][k].angMomCorr[l]/=vox[i][j][k].nMol;
						for(m=0;m<3;m++) {
							vox[i][j][k].rotCorr[m][l]/=vox[i][j][k].nMol;
							vox[i][j][k].angMomProjCorr[m][l]/=vox[i][j][k].nMol;
							vox[i][j][k].atVelCorr[m][l]/=vox[i][j][k].nMol;
						}
					}
					vox[i][j][k].polV.x/=vox[i][j][k].nMol;
					vox[i][j][k].polV.y/=vox[i][j][k].nMol;
					vox[i][j][k].polV.z/=vox[i][j][k].nMol;
					vecRNorm(vox[i][j][k].polV,&vox[i][j][k].polNorm);
					vox[i][j][k].polToProt/=vox[i][j][k].nMol;
					vox[i][j][k].OHtoProt/=vox[i][j][k].nMol;
					vox[i][j][k].molPlaneToProt/=vox[i][j][k].nMol;
				}
			}
		}
	}

	/******************************************************************/
	/* OUTPUT OF CUBE FILES FOR SOME PART OF THE DATA                 */
	/******************************************************************/
	voxel.a.x=dgrid; voxel.a.y=0.0; voxel.a.z=0.0;
	voxel.b.x=0.0; voxel.b.y=dgrid; voxel.b.z=0.0;
	voxel.c.x=0.0; voxel.c.y=0.0; voxel.c.z=dgrid;

	sprintf(volTitle,"CUBE file of average water number density (nm^-3)",delayTime);
        sprintf(fnVol,"%s_water3D-numDens.cube",fnOut,delayTime);
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[i][j][k].averMol/(dgrid*dgrid*dgrid*0.001);
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average polarization (length of averaged normalized dipole moment vector) (0 to 1)",delayTime);
        sprintf(fnVol,"%s_water3D-pol.cube",fnOut,delayTime);
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[i][j][k].polNorm;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average polarization towards protein (-1 to +1)",delayTime);
        sprintf(fnVol,"%s_water3D-pol2Prot.cube",fnOut,delayTime);
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[i][j][k].polToProt;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average OH-orientation (-1 to +1)",delayTime);
	/*-1 not possible, actual minimum value reached when both hydrogen are always */
	/*equally far away from protein => cos(180 degree - 0.5 HOH angle in water) */
        sprintf(fnVol,"%s_water3D-OH2Prot.cube",fnOut,delayTime);
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[i][j][k].OHtoProt;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average mol plane orientation (0 to +1)",delayTime);
        sprintf(fnVol,"%s_water3D-plane2Prot.cube",fnOut,delayTime);
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[i][j][k].molPlaneToProt;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of water COM MSD's after %.3fps in A^2",delayTime);
	sprintf(fnVol,"%s_water3D-COMmsd-%.3fps.cube",fnOut,delayTime);
	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
		volData[i][j][k]=vox[i][j][k].COMmsd[nDelay-1];
	}}}
	printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	for(m=0;m<3;m++) {
		sprintf(volTitle,"CUBE file of water rot. corr. (normalized) around axis %d after %.3fps",m+1,delayTime);
        	sprintf(fnVol,"%s_water3D-rot%d-%.3fps.cube",fnOut,m+1,delayTime);
        	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
			if(vox[i][j][k].rotCorr[m][0]!=0.0) {
				volData[i][j][k]=vox[i][j][k].rotCorr[m][nDelay-1]/vox[i][j][k].rotCorr[m][0];
			} else volData[i][j][k]=0.0;
        	}}}
        	printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);
	}

        /******************************************************************/
        /* OUTPUT OF AVERAGE WATER PROPERTY DATA                          */
        /******************************************************************/

	/* num density (A^-3)*/
	sprintf(fnVol,"%s_water3D-numDens.dat",fnOut,tstep,delayTime);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"%9.5f %9.5f %9.5f %10.5e\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,
			vox[i][j][k].averMol/(dgrid*dgrid*dgrid*0.001));
        }}}
        fclose(out);

	/* polarization (averaged normalized vector and its length, 0 to 1) */
	/* averaged poarization towards protein, OH-bond direction preference */
	/* and molecular plane orientation preference*/
        sprintf(fnVol,"%s_water3D-pol-and-orientation.dat",fnOut,tstep,delayTime);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"%9.5f %9.5f %9.5f %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",
                        gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,
                        vox[i][j][k].polV.x,vox[i][j][k].polV.y,vox[i][j][k].polV.z,vox[i][j][k].polNorm,
			vox[i][j][k].polToProt,vox[i][j][k].OHtoProt,vox[i][j][k].molPlaneToProt);
        }}}
        fclose(out);

        /******************************************************************/
        /* OUTPUT OF DATA FILES FOR ALL TIME DOMAIN DATA                  */
        /******************************************************************/

/*	sprintf(fnVol,"%s_water3D-COMmsd_dt-%.3fps_maxt-%.3fps.dat",fnOut,tstep,delayTime);
	out=fopen(fnVol,"w");
	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
		fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
		for(l=0;l<nDelay;l++) {
			fprintf(out," %10.5e",vox[i][j][k].COMmsd[l]);
		}
		fprintf(out,"\n");
	}}}
	fclose(out);

	sprintf(fnVol,"%s_water3D-COMvel-TCF_dt-%.3fps_maxt-%.3fps.dat",fnOut,tstep,delayTime);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",vox[i][j][k].vCOMcorr[l]);
                }
                fprintf(out,"\n");
        }}}
        fclose(out);

	sprintf(fnVol,"%s_water3D-AngMom-TCF_dt-%.3fps_maxt-%.3fps.dat",fnOut,tstep,delayTime);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",vox[i][j][k].angMomCorr[l]);
                }
                fprintf(out,"\n");
        }}}
        fclose(out);

	for(m=0;m<3;m++) {
		sprintf(fnVol,"%s_water3D-Rot%d-TCF_dt-%.3fps_maxt-%.3fps.dat",fnOut,m+1,tstep,delayTime);
        	out=fopen(fnVol,"w");
        	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
        	        fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
				gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
        	        for(l=0;l<nDelay;l++) {
        	                fprintf(out," %10.5e",vox[i][j][k].rotCorr[m][l]);
        	        }
        	        fprintf(out,"\n");
        	}}}
        	fclose(out);

		sprintf(fnVol,"%s_water3D-AngMomProj%d-TCF_dt-%.3fps_maxt-%.3fps.dat",fnOut,m+1,tstep,delayTime);
                out=fopen(fnVol,"w");
                for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                        fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
				gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
                        for(l=0;l<nDelay;l++) {
                                fprintf(out," %10.5e",vox[i][j][k].angMomProjCorr[m][l]);
                        }
                        fprintf(out,"\n");
                }}}
                fclose(out);
	}

	sprintf(fnVol,"%s_water3D-OxyVel-TCF_dt-%.3fps_maxt-%.3fps.dat",fnOut,tstep,delayTime);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",vox[i][j][k].atVelCorr[0][l]);
                }
                fprintf(out,"\n");
        }}}
        fclose(out);

	sprintf(fnVol,"%s_water3D-HydVel-TCF_dt-%.3fps_maxt-%.3fps.dat",fnOut,tstep,delayTime);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f dt= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,tstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",0.5*(vox[i][j][k].atVelCorr[1][l]+vox[i][j][k].atVelCorr[2][l]));
                }
                fprintf(out,"\n");
        }}}
        fclose(out);
*/
        /******************************************************************/
        /* FOURIER TRANSFORM OF TIME DOMAIN DATA                          */
        /******************************************************************/

	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
		symFT(ft,array,vox[i][j][k].vCOMcorr,nDelay,FTnorm);
	}}}

	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                symFT(ft,array,vox[i][j][k].angMomCorr,nDelay,FTnorm);
        }}}

	for(m=0;m<3;m++) {
		for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
			symFT(ft,array,vox[i][j][k].rotCorr[m],nDelay,FTnorm);
		}}}
		for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                        symFT(ft,array,vox[i][j][k].angMomProjCorr[m],nDelay,FTnorm);
                }}}
		for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                        symFT(ft,array,vox[i][j][k].atVelCorr[m],nDelay,FTnorm);
                }}}
	}

	/******************************************************************/
        /* OUTPUT OF DATA FILES FOR ALL FREQUENCY DOMAIN DATA             */
        /******************************************************************/

	sprintf(fnVol,"%s_water3D-COMvel-VDOS_df-%.3fwn_maxf-%.3fwn.dat",fnOut,fstep,nDelay*fstep);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f df= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,fstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",vox[i][j][k].vCOMcorr[l]);
                }
                fprintf(out,"\n");
        }}}
        fclose(out);

/*	sprintf(fnVol,"%s_water3D-AngMom-VDOS_df-%.3fwn_maxf-%.3fwn.dat",fnOut,fstep,nDelay*fstep);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f df= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,fstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",vox[i][j][k].angMomCorr[l]);
                }
                fprintf(out,"\n");
        }}}
        fclose(out);
*/
	for(m=0;m<3;m++) {
/*		sprintf(fnVol,"%s_water3D-Rot%d-TCF-FT_df-%.3fwn_maxf-%.3fwn.dat",fnOut,m+1,fstep,nDelay*fstep);
        	out=fopen(fnVol,"w");
        	for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
        	        fprintf(out,"# %9.5f %9.5f %9.5f df= %9.5f n= %d\n",
				gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,fstep,vox[i][j][k].nMol);
        	        for(l=0;l<nDelay;l++) {
        	                fprintf(out," %10.5e",vox[i][j][k].rotCorr[m][l]);
        	        }
        	        fprintf(out,"\n");
        	}}}
        	fclose(out);
*/
		sprintf(fnVol,"%s_water3D-AngMomProj%d-VDOS_df-%.3fwn_maxf-%.3fwn.dat",fnOut,m+1,fstep,nDelay*fstep);
                out=fopen(fnVol,"w");
                for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                        fprintf(out,"# %9.5f %9.5f %9.5f df= %9.5f n= %d\n",
				gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,fstep,vox[i][j][k].nMol);
                        for(l=0;l<nDelay;l++) {
                                fprintf(out," %10.5e",vox[i][j][k].angMomProjCorr[m][l]);
                        }
                        fprintf(out,"\n");
                }}}
                fclose(out);
	}

        sprintf(fnVol,"%s_water3D-OxyVel-VDOS_df-%.3fwn_maxf-%.3fwn.dat",fnOut,fstep,nDelay*fstep);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f df= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,fstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",vox[i][j][k].atVelCorr[0][l]);
                }
                fprintf(out,"\n");
        }}}
        fclose(out);

        sprintf(fnVol,"%s_water3D-HydVel-VDOS_df-%.3fwn_maxf-%.3fwn.dat",fnOut,fstep,nDelay*fstep);
        out=fopen(fnVol,"w");
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                fprintf(out,"# %9.5f %9.5f %9.5f df= %9.5f n= %d\n",
			gridOri.x+(i+0.5)*dgrid,gridOri.y+(j+0.5)*dgrid,gridOri.z+(k+0.5)*dgrid,fstep,vox[i][j][k].nMol);
                for(l=0;l<nDelay;l++) {
                        fprintf(out," %10.5e",0.5*(vox[i][j][k].atVelCorr[1][l]+vox[i][j][k].atVelCorr[2][l]));
                }
                fprintf(out,"\n");
        }}}
        fclose(out);

	return 0;
}

