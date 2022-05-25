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
        int	nAtom;
        double	averAtom;
	double	mass;
} t_voxelProp;

int initVoxelProp(t_voxelProp *p) {
        p[0].nAtom=0;
	p[0].averAtom=0.0;
        p[0].mass=0.0;
        return 0;
}

int gridOrigin(int *ng,real dgrid,t_vecR *ori) {
	ori[0].x=(((real)ng[0])/-2.0)*dgrid;
	ori[0].y=(((real)ng[1])/-2.0)*dgrid;
	ori[0].z=(((real)ng[2])/-2.0)*dgrid;
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
	printf("ATOM NUMBER DENSITIES AND MASS DENSITIES IN THE\n");
	printf("ENVIRONMENT OF A SOLUTE FOR A SELECTED STATIC\n");
	printf("GROUP ON A GRID REPORTED IN A CUBE FILE\n");
        printf("\n");
        printf("Version 1.1: July 11, 2019\n");
        printf("Author:\n");
        printf(" Dr. Matthias Heyden\n");
        printf(" Arizona State University\n");
        printf(" School of Molecular Sciences\n");
        printf(" Tempe, AZ, USA\n");
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
        printf("alignGrp\ndensGrp\nrotAlign\nwrap\n");
	printf("ngx\nngy\nngz\ndgrid\nstemName\n");
        return 0;
}

int getInput(char *fnCOM,char *fnTop,char *fnCrd,char *fnVel,char *fnJob,
                int *nRead,int *nSample,
                char *fnRef,int *alignGrp,int *densGrp,
		int *rotAlign,int *wrap,int *ng,real *dgrid,
		char *fnOut,
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
        sscanf(buffer,"%d",densGrp);
        printf("%2d -> read %20s : %d\n",i,"densGrp",densGrp[0]);i++;

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
	printf("%2d -> read %20s : %f\n",i,"dGrid",tmp);i++;
	dgrid[0]=tmp;

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

	int densGrp;
	int ng[3];
	real dgrid;
	t_vecR gridOri;
	t_voxelProp *vox;
	int ng3,voxIdx;
	int pos[3];
	int cnt;
	real ***volData;
	t_box voxel;
	char volTitle[300];
	char fnVol[300];
	real tmp;
	int atomIdx;

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
                &nRead,&nSample,fnRef,&alignGrp,&densGrp,
                &rotAlign,&wrap,ng,&dgrid,
                fnOut,
                needVelo);
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
	if(alignGrp!=-1) getAlignRefSpec(fnRef,grps[alignGrp],atoms,box,&w_rls,&refCrd,&refCOM,&x);
	if(densGrp>=grpDef.nStat) fatal("dens group needs to be static\n");
	gridOrigin(ng,dgrid,&gridOri);
	ng3=ng[0]*ng[1]*ng[2];
	vox=(t_voxelProp*)save_malloc(ng3*sizeof(t_voxelProp));
	for(i=0;i<ng3;i++) {
		initVoxelProp(&vox[i]);
	}

	cnt=0;
	volData=(real***)save_malloc(ng[0]*sizeof(real**));
	for(i=0;i<ng[1];i++) {
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
			if(rotAlign==1) {
                                alignGroup(&grps[alignGrp],atoms,nAtoms,&box,w_rls,refCrd,refCOM,x);
                                getMolCOM(atoms,mols,nMols,needVelo);
                        }
		}
		
		if(i%nSample==0)
                {
                        /*updateGrps(atoms,nAtoms,mols,nMols,box,grpDef,grps);*/
/*MODULAR*/
			for(j=0;j<grps[densGrp].nAtoms;j++) {
				atomIdx=grps[densGrp].atoms[j];
				if(gridPos(atoms[atomIdx].crd,gridOri,dgrid,dgrid,dgrid,pos,box,ng)==1) {
					voxIdx=pos[2]+ng[1]*(pos[1]+ng[0]*pos[0]);
					vox[voxIdx].nAtom++;
					vox[voxIdx].mass+=atoms[atomIdx].mass;
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
		vox[i].averAtom=((double)vox[i].nAtom)/((double)cnt);
		if(vox[i].nAtom!=0) {
			vox[i].mass/=((double)cnt);
		}
	}

	/******************************************************************/
	/* OUTPUT OF CUBE FILES FOR DENSITY DATA                          */
	/******************************************************************/
	voxel.a.x=dgrid; voxel.a.y=0.0; voxel.a.z=0.0;
	voxel.b.x=0.0; voxel.b.y=dgrid; voxel.b.z=0.0;
	voxel.c.x=0.0; voxel.c.y=0.0; voxel.c.z=dgrid;

	sprintf(volTitle,"CUBE file of average atom number density (nm^-3)");
        sprintf(fnVol,"%s_dens3D-numDensAtom.cube",fnOut);
	l=0;
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[l].averAtom/(dgrid*dgrid*dgrid*0.001);
		l++;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	sprintf(volTitle,"CUBE file of average mass density (g/mL)");
        sprintf(fnVol,"%s_dens3D-numDensMass.cube",fnOut);
        l=0;
        for(i=0;i<ng[0];i++) { for(j=0;j<ng[1];j++) { for(k=0;k<ng[2];k++) {
                volData[i][j][k]=vox[l].mass/(dgrid*dgrid*dgrid)*1.66058;
                l++;
        }}}
        printVolumeData(fnVol,volTitle,gridOri,voxel,ng[0],ng[1],ng[2],volData);

	return 0;
}

