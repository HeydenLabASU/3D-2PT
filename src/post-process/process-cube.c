#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int readCube(char *fn,float *grid,int gN) {
	FILE *io;
	char buffer[300];
	int i;
	io=fopen(fn,"r");
	for(i=0;i<7;i++) fgets(buffer,300,io);
	for(i=0;i<gN;i++) fscanf(io,"%f",&grid[i]);
	fclose(io);
	return 0;
}

float distSq(float *a,float *b) {
	int i;
	float tmp;
	float res=0.0;

	for(i=0;i<3;i++) {
		tmp=a[i]-b[i];
		res+=tmp*tmp;
	}
	return res;
}

int printVolumeData(char *fn,char *title,float *ori,float *vox,int nx,int ny,int nz,float *data) {
        FILE *cube;
	int a,b,c;
        char name[100];
	double oriOut[3];

	/*center of voxel[0][0][0], instead of corner (used internally)*/
        oriOut[0]=((double)ori[0])+((double)vox[0])/2.0;
        oriOut[1]=((double)ori[1])+((double)vox[1])/2.0;
        oriOut[2]=((double)ori[2])+((double)vox[2])/2.0;

        cube=fopen(fn,"w");
        fprintf(cube,"%s\n",title);
        fprintf(cube,"created by M. Heyden\n");
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",1,oriOut[0]/0.5292,oriOut[1]/0.5292,oriOut[2]/0.5292);
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",nx,vox[0]/0.5292,0.0,0.0);
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",ny,0.0,vox[1]/0.5292,0.0);
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",nz,0.0,0.0,vox[2]/0.5292);
        fprintf(cube,"%d%12.6f%12.6f%12.6f%12.6f\n",8,0.0,0.0,0.0,0.0);
	
        for (a=0;a<nx;a++) {
                for (b=0;b<ny;b++) {
                        for (c=0;c<nz;c++) {
                                fprintf(cube," %12.5e",(double)data[a*ny*nz+b*nz+c]);
                                if (c % 6 == 5) fprintf(cube,"\n");
                        }
                        fprintf(cube,"\n");
                }
        }
        fclose(cube);
        sprintf(name,"%-50s",fn);
        printf("   wrote file: %s\n",name);
}

int main(int argc,char *argv[]) {
	FILE *io;
	char buffer[300];
	char fileName[100];
	char **atName;
	char **resName;
	float **atCrd;
	int nAtoms;
	int i,j,k,m;
	int nx,ny,nz;
	float gOri[3];
	float gVox[3];
	float vol;
	float dum1,dum2;
	float **gCrd;
	float *dens1,*dens2,*Uxs,*Uss,*dUss,*dH,*totS,*mTdStot,*mTdSxs,*dG,*dGdesolv,*minDist;
	float temp,minDens;
	float bulkUss=0.0;
	float bulkStot=0.0;
	float bulkDens1=0.0;
	float bulkDens2=0.0;
	float bulkDist,bulkDistSq,dSq;
	int *bulkSel;
	int *intSel;
	double integral;
	int nBulkSel=0;
	float norm;
	int gN;
	float dSqMin;

	printf("reference gro file (all atoms):  ");
	scanf("%s",fileName);
	printf("%s\n",fileName);
	
	io=fopen(fileName,"r");
	fgets(buffer,300,io);
	fgets(buffer,300,io);
	sscanf(buffer,"%d",&nAtoms);
	
	atName=(char**)malloc(nAtoms*sizeof(char*));
	resName=(char**)malloc(nAtoms*sizeof(char*));
	atCrd=(float**)malloc(nAtoms*sizeof(float*));
	for(i=0;i<nAtoms;i++) {
		atName[i]=(char*)malloc(10*sizeof(char));
		resName[i]=(char*)malloc(10*sizeof(char));
		atCrd[i]=(float*)malloc(3*sizeof(float));
	}
	for(i=0;i<nAtoms;i++) {
		fgets(buffer,300,io);
		sscanf(&buffer[5],"%s %s",resName[i],atName[i]);
		sscanf(&buffer[20],"%f %f %f",&atCrd[i][0],&atCrd[i][1],&atCrd[i][2]);
		atCrd[i][0]*=10.0;	
		atCrd[i][1]*=10.0;
		atCrd[i][2]*=10.0;
		//printf("%f %f %f\n",atCrd[i][0],atCrd[i][1],atCrd[i][2]);
	}
	printf("read %d atoms\n",nAtoms);
	fclose(io);

	printf("sample cube file (any property): ");
	scanf("%s",fileName);
	printf("%s\n",fileName);
	io=fopen(fileName,"r");
	fgets(buffer,300,io);
	fgets(buffer,300,io);
	fgets(buffer,300,io);	
	sscanf(buffer,"%d %f %f %f",&i,&gOri[0],&gOri[1],&gOri[2]);
	gOri[0]*=0.5292;
	gOri[1]*=0.5292;
	gOri[2]*=0.5292;
	fgets(buffer,300,io);
        sscanf(buffer,"%d %f %f %f",&nx,&gVox[0],&dum1,&dum2);
	gVox[0]*=0.5292;
	fgets(buffer,300,io);
        sscanf(buffer,"%d %f %f %f",&ny,&dum1,&gVox[1],&dum2);
	gVox[1]*=0.5292;
	fgets(buffer,300,io);
        sscanf(buffer,"%d %f %f %f",&nz,&dum1,&dum2,&gVox[2]);
	gVox[2]*=0.5292;
	fgets(buffer,300,io);
	fclose(io);

	/*convert to corner of voxel[0][0][0] for internal use*/
	gOri[0]-=0.5*gVox[0];
	gOri[1]-=0.5*gVox[1];
	gOri[2]-=0.5*gVox[2];

	printf("grid origin (corner): %f %f %f\n",gOri[0],gOri[1],gOri[2]);
	printf("voxel dimensions:     %f %f %f\n",gVox[0],gVox[1],gVox[2]);

	gN=nx*ny*nz;
	dens1=(float*)malloc(gN*sizeof(float));
	dens2=(float*)malloc(gN*sizeof(float));
	Uxs=(float*)malloc(gN*sizeof(float));
	Uss=(float*)malloc(gN*sizeof(float));
	dUss=(float*)malloc(gN*sizeof(float));
	dH=(float*)malloc(gN*sizeof(float));
	totS=(float*)malloc(gN*sizeof(float));
	mTdStot=(float*)malloc(gN*sizeof(float));
	mTdSxs=(float*)malloc(gN*sizeof(float));
	dG=(float*)malloc(gN*sizeof(float));
	dGdesolv=(float*)malloc(gN*sizeof(float));
	minDist=(float*)malloc(gN*sizeof(float));
	bulkSel=(int*)malloc(gN*sizeof(int));
	intSel=(int*)malloc(gN*sizeof(int));
	gCrd=(float**)malloc(gN*sizeof(float*));
	for(i=0;i<gN;i++) {
		gCrd[i]=(float*)malloc(3*sizeof(float));
		bulkSel[i]=1;
		intSel[i]=0;
	}
	m=0;
	for(i=0;i<nx;i++) {
		for(j=0;j<ny;j++) {
			for(k=0;k<nz;k++) {
				gCrd[m][0]=gOri[0]+(i+0.5)*gVox[0];
				gCrd[m][1]=gOri[1]+(j+0.5)*gVox[1];
				gCrd[m][2]=gOri[2]+(k+0.5)*gVox[2];
				//printf("%f %f %f\n",gCrd[m][0],gCrd[m][1],gCrd[m][2]);
				m++;
			}
		}
	}
	printf("bulk distance to closest atom:   ");
        scanf("%f",&bulkDist);
	printf(" %f\n",bulkDist);
	bulkDistSq=bulkDist*bulkDist;
	for(i=0;i<gN;i++) {
		dSqMin=1000000.0;
		for(j=0;j<nAtoms;j++) {
			dSq=distSq(gCrd[i],atCrd[j]);
			//printf("%f  ",dSq);
			if(dSq<=bulkDistSq) bulkSel[i]=0;
			if(dSq<dSqMin) dSqMin=dSq;
		}
		minDist[i]=sqrt(dSqMin);
	}
	for(i=0;i<gN;i++) {
		if(bulkSel[i]==1) nBulkSel++;
	}
	printf("assigned %d voxels as bulk-like\n",nBulkSel);

	printf("pot3D density cube file:         ");
	scanf("%s",fileName);
	printf("%s\n",fileName);
	readCube(fileName,dens1,gN);
	for(i=0;i<gN;i++) {
		bulkDens1+=bulkSel[i]*dens1[i];
	}
	bulkDens1/=nBulkSel;
	printf("bulk number density (pot3D):   %10.5f nm^-3\n",bulkDens1);

	printf("pot3D Uxs cube file:             ");
        scanf("%s",fileName);
	printf("%s\n",fileName);
        readCube(fileName,Uxs,gN);

	printf("pot3D Uss cube file:             ");
        scanf("%s",fileName);
	printf("%s\n",fileName);
        readCube(fileName,Uss,gN);
	norm=0.0;
	for(i=0;i<gN;i++) {
                bulkUss+=bulkSel[i]*dens1[i]*Uss[i];
		norm+=bulkSel[i]*dens1[i];
        }
        bulkUss/=norm;
        printf("bulk Uss: %10.5f kJ/mol\n",bulkUss);

	printf("water3D density cube file:       ");
        scanf("%s",fileName);
	printf("%s\n",fileName);
        readCube(fileName,dens2,gN);
	for(i=0;i<gN;i++) {
                bulkDens2+=bulkSel[i]*dens2[i];
        }
        bulkDens2/=nBulkSel;
        printf("bulk number density (water3D): %10.5f nm^-3\n",bulkDens2);

	printf("water3D totS cube file:          ");
        scanf("%s",fileName);
	printf("%s\n",fileName);
        readCube(fileName,totS,gN);
	norm=0.0;
	for(i=0;i<gN;i++) {
                bulkStot+=bulkSel[i]*dens2[i]*totS[i];
		norm+=bulkSel[i]*dens2[i];
        }
        bulkStot/=norm;
        printf("bulk water entropy (water3D):  %10.5f J/(K mol)\n",bulkStot);

	printf("simulation temperature (Kelvin): ");
	scanf("%f",&temp);
	
	vol=gVox[0]*gVox[1]*gVox[2]*0.001;
	for(i=0;i<gN;i++) {
		dUss[i]=0.5*(Uss[i]-bulkUss);
		dH[i]=Uxs[i]+dUss[i];
		mTdStot[i]=-1.0*temp*(totS[i]-bulkStot)*0.001;
		mTdSxs[i]=mTdStot[i]+dUss[i];
		dG[i]=Uxs[i]+mTdSxs[i];
		dGdesolv[i]=dG[i]*-0.5*(dens1[i]+dens2[i])*vol;
	}

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-minDist.cube",nx,ny,nz);
        printVolumeData(fileName,"CUBE file of minimum distance to solute (A)",gOri,gVox,nx,ny,nz,minDist);

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-dUss.cube",nx,ny,nz);
	printVolumeData(fileName,"CUBE file of 0.5*(Uss-UssBulk) per solvent molecule (kJ/mol)",gOri,gVox,nx,ny,nz,dUss);

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-dH.cube",nx,ny,nz);
        printVolumeData(fileName,"CUBE file of Uxs+0.5*(Uss-UssBulk) per solvent molecule (kJ/mol)",gOri,gVox,nx,ny,nz,dH);

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-mTdStot.cube",nx,ny,nz);
	printVolumeData(fileName,"CUBE file of -TdS (total) per water molecule (kJ/mol)",gOri,gVox,nx,ny,nz,mTdStot);

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-mTdSxs.cube",nx,ny,nz);
	printVolumeData(fileName,"CUBE file of -TdSxs (solute-solvent) per water molecule (kJ/mol)",gOri,gVox,nx,ny,nz,mTdSxs);

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-dG.cube",nx,ny,nz);
	printVolumeData(fileName,"CUBE file of dG per water molecule (kJ/mol)",gOri,gVox,nx,ny,nz,dG);

	sprintf(fileName,"grid-%dx%dx%d_3D-2PT-dG-desolv.cube",nx,ny,nz);
	printVolumeData(fileName,"CUBE file of desolvation deltaG per voxel (kJ/mol)",gOri,gVox,nx,ny,nz,dGdesolv);

	printf("min fraction of bulk density:    ");
	scanf("%f",&minDens);
	for(i=0;i<gN;i++) {
		if((dens1[i]+dens2[i])*0.5>=minDens*0.5*(bulkDens1+bulkDens2)) {
			intSel[i]=1;
		}
	}

	sprintf(fileName,"results_bulk-distance-to-atom-%.1fA_minDensRatio-%.2f.dat",bulkDist,minDens);
	io=fopen(fileName,"w");
	fprintf(io,"bulk_dens (nm^-3)\n%14.5e\n",0.5*(bulkDens1+bulkDens2));
	fprintf(io,"bulk Uss (kJ/mol)\n%14.5e\n",bulkUss);
	fprintf(io,"bulk_S_total (J/K/mol)\n%14.5e\n",bulkStot);
	integral=0.0;
	for(i=0;i<gN;i++) {
		integral+=(dens1[i]*Uxs[i]+dens2[i]*mTdSxs[i])*vol;
	}
	fprintf(io,"dG_solv (kJ/mol)\n%14.5e\n",integral);
	integral=0.0;
        for(i=0;i<gN;i++) {
                integral+=(dens1[i]*(Uxs[i]+dUss[i]))*vol;
        }
        fprintf(io,"dH_solv (kJ/mol)\n%14.5e\n",integral);
	integral=0.0;
        for(i=0;i<gN;i++) {
                integral+=(dens2[i]*mTdStot[i])*vol;
        }
        fprintf(io,"-TdS_solv (kJ/mol)\n%14.5e\n",integral);
	integral=0.0;
        for(i=0;i<gN;i++) {
                integral+=(dens1[i]*Uxs[i])*vol;
        }
        fprintf(io,"dHxs_solv (kJ/mol)\n%14.5e\n",integral);
	integral=0.0;
        for(i=0;i<gN;i++) {
                integral+=(dens2[i]*mTdSxs[i])*vol;
        }
        fprintf(io,"-TdSxs_solv (kJ/mol)\n%14.5e\n",integral);
	fclose(io);

	return 0;
}
