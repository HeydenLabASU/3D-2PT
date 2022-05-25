#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../../include/dataTypes.h"
#include "../../include/fatal.h"
#include "../../include/alloc.h"
#include "../../include/io.h"

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
        printf("nFiles\nlist of filenames\nnPts/data points per line (nDelay/nSample)\n");
	printf("ngx ngy ngz\noutput filename\n");
        return 0;
}

int getInput(char *fnCOM,int *nFiles,char ***fn,int *nPts,int *ng,char *fnOut) {
	FILE *io;
        char buffer[300];
	int i=1;
	int j;

        saveOpenRead(&io,fnCOM);
	
	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nFiles);
        printf("%3d -> read %20s : %d\n",i,"nFiles",nFiles[0]);i++;

	fn[0]=(char**)save_malloc(nFiles[0]*sizeof(char*));
	for(j=0;j<nFiles[0];j++) {
		fn[0][j]=(char*)save_malloc(300*sizeof(char));
		getLineFromCOM(io,buffer,300);
		sscanf(buffer,"%s",fn[0][j]);
        	printf("%3d -> read %15s #%3d : %s\n",i,"input file",j+1,fn[0][j]);i++;
	}

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",nPts);
        printf("%3d -> read %20s : %d\n",i,"data points per line",nPts[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d %d %d",&ng[0],&ng[1],&ng[2]);
        printf("%3d -> read %20s : %d %d %d\n",i,"ngx ngy ngz",ng[0],ng[1],ng[2]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%3d -> read %20s : %s\n",i,"output file",fnOut);i++;

	return 0;
}

int main(int argc, char *argv[]) {
	FILE *io;
	char buffer[500];
	int nFiles;
	char **fn;
	int nPts;
	int ng[3];
	int ngT;
	char fnOut[100];
	char fnCOM[100];
	float vec[3];
	t_vecR *crds;
	t_vecR *check;
	char dtdf1[10],dtdf2[10],dum[10];
	float delta1,delta2;
	int weight;
	int *nMol;
	float dataPt;
	double **dataPts;
	int i,j,k;

	if(argc!=2) {
                printKeys();
                fatal("no input file specified\n");
        }
        strcpy(fnCOM,argv[1]);
	getInput(fnCOM,&nFiles,&fn,&nPts,ng,fnOut);
	ngT=ng[0]*ng[1]*ng[2];

	crds=(t_vecR*)save_malloc(ngT*sizeof(t_vecR));
	nMol=(int*)save_malloc(ngT*sizeof(int));
	dataPts=(double**)save_malloc(ngT*sizeof(double*));
	for(i=0;i<ngT;i++) {
		nMol[i]=0;
		dataPts[i]=(double*)save_malloc(nPts*sizeof(double));
		for(j=0;j<nPts;j++) dataPts[i][j]=0.0;
	}

	for(i=0;i<nFiles;i++) {
		printf("opening file: %s\n",fn[i]);fflush(stdout);
		saveOpenRead(&io,fn[i]);
		for(j=0;j<ngT;j++) {
			fgets(buffer,500,io);
			if(strncmp(buffer,"#",1)!=0) {
				fatal("wrong file format\n");
			}
			if(sscanf(&buffer[1],"%f %f %f %s %f %s %d",
				&vec[0],&vec[1],&vec[2],dtdf1,&delta1,dum,&weight)!=7) fatal("wrong file format");
			if(i==0) {
				crds[j].x=vec[0];crds[j].y=vec[1];crds[j].z=vec[2];
			} else if(crds[j].x!=vec[0] || crds[j].y!=vec[1] || crds[j].z!=vec[2]) {
				fatal("grid point coordinates don't match in input files\n");
			}
			if(i==0 && j==0) {
				strcpy(dtdf2,dtdf1);
				delta2=delta1;
			} else if(strcmp(dtdf2,dtdf1)!=0 || delta1!=delta1) {
				fatal("time/frequency steps don't match in input files\n");
			}
			nMol[j]+=weight;
			for(k=0;k<nPts;k++) {
				fscanf(io,"%f",&dataPt);
				dataPts[j][k]+=((double)weight*dataPt);
			}
			fgets(buffer,500,io);
		}
		fclose(io);
	}
	for(i=0;i<ngT;i++) {
		if(nMol[i]!=0) {
			for(j=0;j<nPts;j++) {
				dataPts[i][j]/=((double)nMol[i]);
			}
		}
	}
	io=fopen(fnOut,"w");
	for(i=0;i<ngT;i++) {
		fprintf(io,"# %9.5f %9.5f %9.5f %s %9.5f n= %d\n",
                        crds[i].x,crds[i].y,crds[i].z,dtdf1,delta1,nMol[i]);
                for(j=0;j<nPts;j++) {
                        fprintf(io," %10.5e",dataPts[i][j]);
                }
                fprintf(io,"\n");
	}
	fclose(io);

	return 0;
}

