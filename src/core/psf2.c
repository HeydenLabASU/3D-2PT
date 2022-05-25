#include <stdio.h>
#include <string.h>
#include "../../include/fatal.h"
#include "../../include/dataTypes.h"
#include "../../include/alloc.h"
#include "../../include/io.h"

typedef struct {
	int	atomNr;
	char	chain[6];
	int	resNr;
	char	resName[6];
	char	atomName[6];
	char	atomType[6];
	float	charge;
	float	mass;
} t_atomLine;

int goToSection(FILE *in,char *keyword, int *n)
{
	char buffer[300];

	fgets(buffer,300,in);
	while(strstr(buffer,keyword)==NULL) fgets(buffer,300,in);

	sscanf(buffer,"%d",n);

	return 0;
}

int readAtoms(FILE *in,t_atomLine *atomInfo,int nAtoms)
{
	int i;
	char buffer[300];

	for(i=0;i<nAtoms;i++)
	{
		fgets(buffer,300,in);
		sscanf(buffer,"%d %s %d %s %s %s %f %f",&atomInfo[i].atomNr,atomInfo[i].chain,&atomInfo[i].resNr,atomInfo[i].resName,atomInfo[i].atomName,atomInfo[i].atomType,&atomInfo[i].charge,&atomInfo[i].mass);
	}

	return 0;
}

int readBonds(FILE *in,t_bond *bonds,int n)
{
	int i;
	char buffer[300];

	for(i=0;i<n;i++) fscanf(in,"%d %d",&bonds[i][0],&bonds[i][1]);
	fgets(buffer,300,in);

	return 0;
}

int readAngles(FILE *in,t_angle *angles,int n)
{
        int i;
        char buffer[300];

        for(i=0;i<n;i++) fscanf(in,"%d %d %d",&angles[i][0],&angles[i][1],&angles[i][2]);
        fgets(buffer,300,in);

        return 0;
}

int readDihedrals(FILE *in,t_dihed *diheds,int n)
{
        int i;
        char buffer[300];

        for(i=0;i<n;i++) fscanf(in,"%d %d %d %d",&diheds[i][0],&diheds[i][1],&diheds[i][2],&diheds[i][3]);
        fgets(buffer,300,in);

        return 0;
}

int readImpropers(FILE *in,t_improp *improps,int n)
{
        int i;
        char buffer[300];

        for(i=0;i<n;i++) fscanf(in,"%d %d %d %d",&improps[i][0],&improps[i][1],&improps[i][2],&improps[i][3]);
        fgets(buffer,300,in);

        return 0;
}

int main(int argc,char *argv[])
{
	FILE *io;
	char err[300];
	t_atomLine *atoms;
	t_bond *bonds;
	t_angle *angles;
	t_dihed *dihedrals;
	t_improp *impropers;
	int *bTag,*aTag,*dTag,*iTag;
	int i,j,k,l,b,e,nAtoms,nBonds,nAngles,nDihed,nImprop;
	int nMols;
	char curChain[6];
	char **at;
	int nAt=0;
	real *atm;
	int *resNr;

	if(argc!=3) fatal("WRONG COMMAND LINE ARGUMENTS:\nNEEDED: input.psf   output.top\n -> main\n -> psf.c\n");
	if(strcmp(argv[1],argv[2])==0) fatal("COWARDLY REFUSING TO OVERWRITE INPUT FILE\n -> main\n -> psf.c\n");
	

	saveOpenRead(&io,argv[1]);

	goToSection(io,"!NATOM",&nAtoms);
	atoms=(t_atomLine*)save_malloc(nAtoms*sizeof(t_atomLine));
	readAtoms(io,atoms,nAtoms);

	goToSection(io,"!NBOND",&nBonds);
        bonds=(t_bond*)save_malloc(nBonds*sizeof(t_bond));
	bTag=(int*)save_malloc(nBonds*sizeof(int));
	for(i=0;i<nBonds;i++) bTag[i]=0;
        readBonds(io,bonds,nBonds);
	
	goToSection(io,"!NTHETA",&nAngles);
        angles=(t_angle*)save_malloc(nAngles*sizeof(t_angle));
	aTag=(int*)save_malloc(nAngles*sizeof(int));
        for(i=0;i<nAngles;i++) aTag[i]=0;
        readAngles(io,angles,nAngles);

	goToSection(io,"!NPHI",&nDihed);
        dihedrals=(t_dihed*)save_malloc(nDihed*sizeof(t_dihed));
	dTag=(int*)save_malloc(nDihed*sizeof(int));
        for(i=0;i<nDihed;i++) dTag[i]=0;
        readDihedrals(io,dihedrals,nDihed);

	goToSection(io,"!NIMPHI",&nImprop);
        impropers=(t_improp*)save_malloc(nImprop*sizeof(t_improp));
	iTag=(int*)save_malloc(nImprop*sizeof(int));
        for(i=0;i<nImprop;i++) iTag[i]=0;
        readImpropers(io,impropers,nImprop);

	/*printf("nBonds:\t\t%d\nnAngles:\t%d\nnDihedrals:\t%d\nnImpropers:\t%d\n",nBonds,nAngles,nDihed,nImprop);*/

	fclose(io);

	nMols=0;
	j=-1;
	strcpy(curChain,"X123X");
	for(i=0;i<nAtoms;i++)
	{
		if(atoms[i].resNr!=j || strcmp(atoms[i].chain,curChain)!=0)
		{
			j=atoms[i].resNr;
			strcpy(curChain,atoms[i].chain);
			nMols++;
		}
	}

	at=(char**)save_malloc(nAtoms*sizeof(char*));
	for(i=0;i<nAtoms;i++) at[i]=(char*)save_malloc(6*sizeof(char));
	atm=(real*)save_malloc(nAtoms*sizeof(real));
	resNr=(int*)save_malloc(nAtoms*sizeof(int));
	j=0;
	for(i=0;i<nAtoms;i++)
        {
		k=1;
		l=0;
		while(l<j && k!=0)
		{
			k=strcmp(atoms[i].atomType,at[l]);
			l++;
		}
		if(k!=0)
		{
			strcpy(at[j],atoms[i].atomType);
			atm[j]=atoms[i].mass;
			j++;
		}
	}
	nAt=j;

	j=atoms[0].resNr;
        strcpy(curChain,atoms[0].chain);
	b=1;
        for(i=0;i<=nAtoms;i++)
        {
                if(atoms[i].resNr!=j || strcmp(curChain,atoms[i].chain)!=0)
                {
			j=atoms[i].resNr;
			b++;
		}
		resNr[i]=b;
	}

	io=fopen(argv[2],"w");
	fprintf(io,";File %s\n;converted from: %s\n;by converPSF2.exe\n",argv[2],argv[1]);
	fprintf(io,"\n[ defaults ]\nLJ       Geometric\n\n");
	fprintf(io,"[ atomtypes ]\n;name        mass      charge   ptype   c6      c12\n");
	for(i=0;i<nAt;i++)
	{
		fprintf(io,"%5s%12.5f       0.0     A       0.0     0.0\n",at[i],atm[i]);
	}
	fprintf(io,"\n\n");

/*	fprintf(io,"[ nonbond_params ]\n");
	fprintf(io,"\n");*/

	fprintf(io,"[ moleculetype ]\n");
	fprintf(io,"; molname       nrexcl\nSYS       3\n\n");

	fprintf(io,"[ atoms ]\n");
	fprintf(io,"; id    at type res nr  residu name     at name         cg nr   charge\n");
	for(i=0;i<nAtoms;i++)
	{
		fprintf(io,"%-8d%-8s%-8d%-16s%-16s%-8d%-8.5f\n",i+1,atoms[i].atomType,resNr[i],atoms[i].resName,atoms[i].atomName,i+1,atoms[i].charge);
	}
	fprintf(io,"\n\n");

	fprintf(io,"[ bonds ]\n");
	for(i=0;i<nBonds;i++)
	{
		fprintf(io,"%-8d%-8d%-4d%-12.5f%-12.5f\n",bonds[i][0],bonds[i][1],2,0.100,10000.0);
	}
	fprintf(io,"\n\n");

	fprintf(io,"[ angles ]\n");
        for(i=0;i<nAngles;i++)
        {
                fprintf(io,"%-8d%-8d%-8d%-4d%-12.5f%-12.5f\n",angles[i][0],angles[i][1],angles[i][2],2,100.0,100.0);
        }
        fprintf(io,"\n\n");

	fprintf(io,"[ dihedrals ]\n");
        for(i=0;i<nDihed;i++)
        {
                fprintf(io,"%-8d%-8d%-8d%-8d%-4d%-12.5f%-12.5f%-12.5f%-12.5f%-12.5f%-12.5f\n",dihedrals[i][0],dihedrals[i][1],dihedrals[i][2],dihedrals[i][3],1,0.0,0.0,0.0,0.0,0.0,0.0);
        }
        fprintf(io,"\n\n");

	fprintf(io,"[ dihedrals ]\n");
        for(i=0;i<nImprop;i++)
        {
                fprintf(io,"%-8d%-8d%-8d%-8d%-4d%-12.5f%-12.5f%-12.5f%-12.5f\n",impropers[i][0],impropers[i][1],impropers[i][2],impropers[i][3],2,0.0,0.0,0.0,0.0);
        }
        fprintf(io,"\n\n");

	fprintf(io,"[ system ]\n; Name\n;Protein in water\n\n");

	fprintf(io,"[ molecules ]\n; Compound      #mols\nSYS       1\n");

	fclose(io);

	return 0;
}

