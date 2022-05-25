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
	int i,j,k,b,e,nAtoms,nBonds,nAngles,nDihed,nImprop;
	int nMols;
	char curChain[6];

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

	io=fopen(argv[2],"w");
	fprintf(io,"#File %s\n#converted from: %s\n#by converPSF.exe\n",argv[2],argv[1]);
	fprintf(io,"nAtoms:\t\t%d\nnSites:\t\t%d\nnMols:\t\t%d\nnMolTypes:\t%d\n",nAtoms,0,nMols,nMols);

        j=atoms[0].resNr;
	strcpy(curChain,atoms[0].chain);
	b=0;
        for(i=0;i<=nAtoms;i++)
        {
                if(atoms[i].resNr!=j || strcmp(curChain,atoms[i].chain)!=0)
                {
			fprintf(io,"[molecule]\t1\nmolName\t\t%s\nchain\t\t%s\n",atoms[b].resName,atoms[b].chain);
			e=i-1;
			fprintf(io,"[atoms]\t\t%d\n",e-b+1);
			for(j=b;j<i;j++) fprintf(io,"%-7s %8.4f   %8.4f   %8.4f\n",atoms[j].atomName,atoms[j].mass,atoms[j].charge,0.0);
			fprintf(io,"[sites]\t\t0\n");

			k=0;
			for(j=0;j<nBonds;j++) {
                                if(bonds[j][0]>=atoms[b].atomNr && bonds[j][0]<=atoms[e].atomNr
                                && bonds[j][1]>=atoms[b].atomNr && bonds[j][1]<=atoms[e].atomNr)
                                {
					k++;
                                        bTag[j]=1;
                                }
                        }
			fprintf(io,"[bonds]\t\t%d\n",k);
			for(j=0;j<nBonds;j++) {
				if(bonds[j][0]>=atoms[b].atomNr && bonds[j][0]<=atoms[e].atomNr
				&& bonds[j][1]>=atoms[b].atomNr && bonds[j][1]<=atoms[e].atomNr)
				{
					fprintf(io,"%d\t%d\n",bonds[j][0]-b,bonds[j][1]-b);
				}
			}

			k=0;
                        for(j=0;j<nAngles;j++) {
                                if(angles[j][0]>=atoms[b].atomNr && angles[j][0]<=atoms[e].atomNr
                                && angles[j][1]>=atoms[b].atomNr && angles[j][1]<=atoms[e].atomNr
                                && angles[j][2]>=atoms[b].atomNr && angles[j][2]<=atoms[e].atomNr)
                                {
					k++;
                                        aTag[j]=1;
                                }
                        }
			fprintf(io,"[angles]\t%d\n",k);
                        for(j=0;j<nAngles;j++) {
                                if(angles[j][0]>=atoms[b].atomNr && angles[j][0]<=atoms[e].atomNr
				&& angles[j][1]>=atoms[b].atomNr && angles[j][1]<=atoms[e].atomNr
				&& angles[j][2]>=atoms[b].atomNr && angles[j][2]<=atoms[e].atomNr)
                                {
                                        fprintf(io,"%d\t%d\t%d\n",angles[j][0]-b,angles[j][1]-b,angles[j][2]-b);
                                }
			}

			k=0;
			for(j=0;j<nDihed;j++) {
                                if(dihedrals[j][0]>=atoms[b].atomNr && dihedrals[j][0]<=atoms[e].atomNr
                                && dihedrals[j][1]>=atoms[b].atomNr && dihedrals[j][1]<=atoms[e].atomNr
                                && dihedrals[j][2]>=atoms[b].atomNr && dihedrals[j][2]<=atoms[e].atomNr
                                && dihedrals[j][3]>=atoms[b].atomNr && dihedrals[j][3]<=atoms[e].atomNr)
                                {
					k++;
                                        dTag[j]=1;
                                }
                        }
			fprintf(io,"[dihedrals]\t%d\n",k);
                        for(j=0;j<nDihed;j++) {
                                if(dihedrals[j][0]>=atoms[b].atomNr && dihedrals[j][0]<=atoms[e].atomNr
				&& dihedrals[j][1]>=atoms[b].atomNr && dihedrals[j][1]<=atoms[e].atomNr
				&& dihedrals[j][2]>=atoms[b].atomNr && dihedrals[j][2]<=atoms[e].atomNr
				&& dihedrals[j][3]>=atoms[b].atomNr && dihedrals[j][3]<=atoms[e].atomNr)
                                {
                                        fprintf(io,"%d\t%d\t%d\t%d\n",dihedrals[j][0]-b,dihedrals[j][1]-b,dihedrals[j][2]-b,dihedrals[j][3]-b);
                                }
                        }

			k=0;
			for(j=0;j<nImprop;j++) {
                                if(impropers[j][0]>=atoms[b].atomNr && impropers[j][0]<=atoms[e].atomNr
                                && impropers[j][1]>=atoms[b].atomNr && impropers[j][1]<=atoms[e].atomNr
                                && impropers[j][2]>=atoms[b].atomNr && impropers[j][2]<=atoms[e].atomNr
                                && impropers[j][3]>=atoms[b].atomNr && impropers[j][3]<=atoms[e].atomNr)
                                {
					k++;
                                        iTag[j]=1;
                                }
                        }
			fprintf(io,"[impropers]\t%d\n",k);
                        for(j=0;j<nImprop;j++) {
                                if(impropers[j][0]>=atoms[b].atomNr && impropers[j][0]<=atoms[e].atomNr                                                    
                                && impropers[j][1]>=atoms[b].atomNr && impropers[j][1]<=atoms[e].atomNr
                                && impropers[j][2]>=atoms[b].atomNr && impropers[j][2]<=atoms[e].atomNr
                                && impropers[j][3]>=atoms[b].atomNr && impropers[j][3]<=atoms[e].atomNr)
                                {
                                        fprintf(io,"%d\t%d\t%d\t%d\n",impropers[j][0]-b,impropers[j][1]-b,impropers[j][2]-b,impropers[j][3]-b);
                                        iTag[j]=1;
                                }
                        }
                        j=atoms[i].resNr;
			strcpy(curChain,atoms[i].chain);
			b=i;
		}
	}
	j=0;
	for(i=0;i<nBonds;i++) j+=bTag[i];
	fprintf(io,"[inter_bonds]\t%d\n",nBonds-j);
	for(i=0;i<nBonds;i++)
	{
		if(bTag[i]==0)  fprintf(io,"%d\t%d\t%s\t%s\n",bonds[i][0],bonds[i][1],atoms[bonds[i][0]].resName,atoms[bonds[i][1]].resName);
	}

	j=0;
        for(i=0;i<nAngles;i++) j+=aTag[i];
        fprintf(io,"[inter_angles]\t%d\n",nAngles-j);
        for(i=0;i<nAngles;i++) 
        {
                if(aTag[i]==0)  fprintf(io,"%d\t%d\t%d\n",angles[i][0],angles[i][1],angles[i][2]);
        }
	
	j=0;
	for(i=0;i<nDihed;i++) j+=dTag[i];
	fprintf(io,"[inter_dihedrals]\t%d\n",nDihed-j);
	for(i=0;i<nDihed;i++)
	{
		if(dTag[i]==0) fprintf(io,"%d\t%d\t%d\t%d\n",dihedrals[i][0],dihedrals[i][1],dihedrals[i][2],dihedrals[i][3]);
	}

	j=0;
        for(i=0;i<nImprop;i++) j+=iTag[i];
        fprintf(io,"[inter_impropers]\t%d\n",nImprop-j);
        for(i=0;i<nImprop;i++)
        {
                if(iTag[i]==0) fprintf(io,"%d\t%d\t%d\t%d\n",impropers[i][0],impropers[i][1],impropers[i][2],impropers[i][3]);
        }

	fprintf(io,"\n");

	fclose(io);

	return 0;
}

