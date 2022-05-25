#include <stdio.h>
#include <string.h>
#include "../../include/fatal.h"
#include "../../include/dataTypes.h"
#include "../../include/io.h"
#include "../../include/alloc.h"

int getLineFromTOP(FILE *in,char *buffer,int *lineCnt)
{
	fgets(buffer,300,in);
	lineCnt[0]++;
	while(strncmp(buffer,"#",1)==0) 
	{
		fgets(buffer,300,in);
		lineCnt[0]++;
	}
	
	return 0;
}

int allocMolContents(t_mol *dest)
{
	allocInts(&dest[0].atoms,dest[0].nAtoms);
        allocSites(&dest[0].sites,dest[0].nSites);
        allocBonds(&dest[0].bonds,dest[0].nBonds);
        allocAngles(&dest[0].angles,dest[0].nAngles);
        allocDihedrals(&dest[0].dihed,dest[0].nDihed);
        allocImpropers(&dest[0].improp,dest[0].nImprop);

	return 0;
}

int getIntFromTOP(char *fn,FILE *in,int *dest,char *keyword,int *lineCnt)
{
	char buffer[300];
	int readCnt;
	char object[100];
	char err[300];

	getLineFromTOP(in,buffer,lineCnt);
	readCnt=sscanf(buffer,"%s %d",object,dest);
	if(strcmp(object,keyword)!=0 || readCnt!=2) {
		sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getIntFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
		fatal(err);
	}
	return 0;
}

int getStringFromTOP(char *fn,FILE *in,char *dest,char *keyword,int *lineCnt)
{
        char buffer[300];
        int readCnt;
        char object[100];
        char err[300];

        getLineFromTOP(in,buffer,lineCnt);
        readCnt=sscanf(buffer,"%s %s",object,dest);
        if(strcmp(object,keyword)!=0 || readCnt!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getStringFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
        return 0;
}

int getAtomLineFromTOP(char *fn,FILE *in,char *name,real *mass,real *charge,real *pol,int *lineCnt)
{
	char buffer[300];
	char err[400];
	float dumMass,dumCharge,dumPol;
	
	getLineFromTOP(in,buffer,lineCnt);
        if(sscanf(buffer,"%s %f %f %f",name,&dumMass,&dumCharge,&dumPol)!=4) {
        	sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getAtomLineFromTOP\n -> top.c",fn,lineCnt[0],buffer);
        	fatal(err);
        }
	mass[0]=(real)dumMass;
	charge[0]=(real)dumCharge;
	pol[0]=(real)dumPol;
	return 0;
}

int getBondLineFromTOP(char *fn,FILE *in,t_bond *bond,int min,int max,int *lineCnt,int bondOffset)
{
	char buffer[300];
	char err[400];

	getLineFromTOP(in,buffer,lineCnt);
	if(sscanf(buffer,"%d %d",&bond[0][0],&bond[0][1])!=2) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getBondLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	bond[0][0]--;
	bond[0][1]--;
	if(bond[0][0]<min || bond[0][1]<min || bond[0][0]>max || bond[0][1]>max) {
		sprintf(err,"BOND OUT OF BOUNDS In FILE: %s\nLINE # %d: %s\n -> getBondLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
		fatal(err);
	}
	bond[0][0]+=bondOffset;
	bond[0][1]+=bondOffset;
	return 0;
}

int getAngleLineFromTOP(char *fn,FILE *in,t_angle *angle,int min,int max,int *lineCnt,int bondOffset)
{
        char buffer[300];
	char err[400];

        getLineFromTOP(in,buffer,lineCnt);
        if(sscanf(buffer,"%d %d %d",&angle[0][0],&angle[0][1],&angle[0][2])!=3) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getAngleLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	angle[0][0]--;
	angle[0][1]--;
	angle[0][2]--;
        if(angle[0][0]<min || angle[0][1]<min || angle[0][2]<min || angle[0][0]>max || angle[0][1]>max || angle[0][2]>max) {
                sprintf(err,"ANGLE OUT OF BOUNDS In FILE: %s\nLINE # %d: %s\n -> getAngleLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	angle[0][0]+=bondOffset;
	angle[0][1]+=bondOffset;
	angle[0][2]+=bondOffset;
        return 0;
}

int getDihedralLineFromTOP(char *fn,FILE *in,t_dihed *dihedral,int min,int max,int *lineCnt,int bondOffset)
{
        char buffer[300];
	char err[400];

        getLineFromTOP(in,buffer,lineCnt);
        if(sscanf(buffer,"%d %d %d %d",&dihedral[0][0],&dihedral[0][1],&dihedral[0][2],&dihedral[0][3])!=4) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getDihedralLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	dihedral[0][0]--;
	dihedral[0][1]--;
	dihedral[0][2]--;
	dihedral[0][3]--;
        if(dihedral[0][0]<min || dihedral[0][1]<min || dihedral[0][2]<min || dihedral[0][3]<min || dihedral[0][0]>max || dihedral[0][1]>max || dihedral[0][2]>max || dihedral[0][3]>max) {
                sprintf(err,"DIHEDRAL OUT OF BOUNDS In FILE: %s\nLINE # %d: %s\n -> getDihedralLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	dihedral[0][0]+=bondOffset;
	dihedral[0][1]+=bondOffset;
	dihedral[0][2]+=bondOffset;
	dihedral[0][3]+=bondOffset;
        return 0;
}

int getImproperLineFromTOP(char *fn,FILE *in,t_improp *improper,int min,int max,int *lineCnt,int bondOffset)
{
        char buffer[300];
	char err[400];

        getLineFromTOP(in,buffer,lineCnt);
        if(sscanf(buffer,"%d %d %d %d",&improper[0][0],&improper[0][1],&improper[0][2],&improper[0][3])!=4) {
                sprintf(err,"WRONG FORMAT IN FILE: %s\nLINE #%d: %s\n -> getImproperLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	improper[0][0]--;
	improper[0][1]--;
	improper[0][2]--;
	improper[0][3]--;
        if(improper[0][0]<min || improper[0][1]<min || improper[0][2]<min || improper[0][3]<min || improper[0][0]>max || improper[0][1]>max || improper[0][2]>max || improper[0][3]>max) {
                sprintf(err,"IMPROPER OUT OF BOUNDS In FILE: %s\nLINE # %d: %s\n -> getImproperLineFromTOP\n -> top.c\n",fn,lineCnt[0],buffer);
                fatal(err);
        }
	improper[0][0]+=bondOffset;
	improper[0][1]+=bondOffset;
	improper[0][2]+=bondOffset;
	improper[0][3]+=bondOffset;
        return 0;
}

int cpyMolTop(t_mol src,t_mol *dest,t_atom *atoms,int copy)
{
	int i;

	dest[0].nAtoms=src.nAtoms;
	dest[0].nSites=src.nSites;
	dest[0].nBonds=src.nBonds;
	dest[0].nAngles=src.nAngles;
	dest[0].nDihed=src.nDihed;
	dest[0].nImprop=src.nImprop;
	dest[0].mass=src.mass;

	allocMolContents(dest);

	strcpy(dest[0].molName,src.molName);
	strcpy(dest[0].chain,src.chain);

	for(i=0;i<src.nAtoms;i++)
	{
		dest[0].atoms[i]=src.atoms[i]+copy*src.nAtoms;
		atoms[dest[0].atoms[i]].resNr=atoms[src.atoms[i]].resNr+copy;
		strcpy(atoms[dest[0].atoms[i]].resName,src.molName);
		strcpy(atoms[dest[0].atoms[i]].atomName,atoms[src.atoms[i]].atomName);
		atoms[dest[0].atoms[i]].atomNr=atoms[src.atoms[i]].atomNr+copy*src.nAtoms;
		atoms[dest[0].atoms[i]].mass=atoms[src.atoms[i]].mass;
		atoms[dest[0].atoms[i]].charge=atoms[src.atoms[i]].charge;
		atoms[dest[0].atoms[i]].pol=atoms[src.atoms[i]].pol;
	}
	for(i=0;i<src.nSites;i++)
	{
		strcpy(dest[0].sites[i].siteName,src.sites[i].siteName);
		dest[0].sites[i].mass=src.sites[i].mass;
		dest[0].sites[i].charge=src.sites[i].charge;
		dest[0].sites[i].pol=src.sites[i].pol;
	}
	for(i=0;i<src.nBonds;i++)
	{
		dest[0].bonds[i][0]=src.bonds[i][0]+copy*src.nAtoms;
		dest[0].bonds[i][1]=src.bonds[i][1]+copy*src.nAtoms;
	}
	for(i=0;i<src.nAngles;i++)
        {
                dest[0].angles[i][0]=src.angles[i][0]+copy*src.nAtoms;
                dest[0].angles[i][1]=src.angles[i][1]+copy*src.nAtoms;
		dest[0].angles[i][2]=src.angles[i][2]+copy*src.nAtoms;
        }
	for(i=0;i<src.nDihed;i++)
        {
                dest[0].dihed[i][0]=src.dihed[i][0]+copy*src.nAtoms;
                dest[0].dihed[i][1]=src.dihed[i][1]+copy*src.nAtoms;
                dest[0].dihed[i][2]=src.dihed[i][2]+copy*src.nAtoms;
		dest[0].dihed[i][3]=src.dihed[i][3]+copy*src.nAtoms;
        }
	for(i=0;i<src.nImprop;i++)
        {
                dest[0].improp[i][0]=src.improp[i][0]+copy*src.nAtoms;
                dest[0].improp[i][1]=src.improp[i][1]+copy*src.nAtoms;
                dest[0].improp[i][2]=src.improp[i][2]+copy*src.nAtoms;
                dest[0].improp[i][3]=src.improp[i][3]+copy*src.nAtoms;
        }

	return 0;
}	

int readTOP(char *fn,t_atom **atoms,int *nAtoms,t_vecR **siteCrds,int *nSites,t_mol **mols,int *nMols)
{
	FILE *in;
	char buffer[300],err[400],object[300];
	int readCnt,atomCnt,molCnt,nMolTypes;
	int i,j,k;
	int copy,ncopy;
	int min,max;
	int lineCnt=0;
	int bondOffset;
	
	saveOpenRead(&in,fn);

	getIntFromTOP(fn,in,nAtoms,"nAtoms:",&lineCnt);
        allocAtoms(atoms,nAtoms[0]);

	getIntFromTOP(fn,in,nSites,"nSites:",&lineCnt);
        allocVecRs(siteCrds,nSites[0]);

	getIntFromTOP(fn,in,nMols,"nMols:",&lineCnt);
	allocMols(mols,nMols[0]+1);

	getIntFromTOP(fn,in,&nMolTypes,"nMolTypes:",&lineCnt);

	molCnt=0;
	atomCnt=0;
	for(i=0;i<nMolTypes;i++) {
		getIntFromTOP(fn,in,&ncopy,"[molecule]",&lineCnt);
		getStringFromTOP(fn,in,mols[0][molCnt].molName,"molName",&lineCnt);
		getStringFromTOP(fn,in,mols[0][molCnt].chain,"chain",&lineCnt);
		getIntFromTOP(fn,in,&mols[0][molCnt].nAtoms,"[atoms]",&lineCnt);
		allocInts(&mols[0][molCnt].atoms,mols[0][molCnt].nAtoms);
		mols[0][molCnt].mass=0.0;
		for(j=0;j<mols[0][molCnt].nAtoms;j++) {
			mols[0][molCnt].atoms[j]=atomCnt;
			atoms[0][atomCnt].resNr=molCnt;
			strcpy(atoms[0][atomCnt].resName,mols[0][molCnt].molName);
			getAtomLineFromTOP(fn,in,atoms[0][atomCnt].atomName,&atoms[0][atomCnt].mass,&atoms[0][atomCnt].charge,&atoms[0][atomCnt].pol,&lineCnt);
			atoms[0][atomCnt].atomNr=atomCnt;
			mols[0][molCnt].mass+=atoms[0][atomCnt].mass;
			atomCnt++;
		}
		min=0;
		max=mols[0][molCnt].nAtoms-1;

		bondOffset=atomCnt-mols[0][molCnt].nAtoms;

		getIntFromTOP(fn,in,&mols[0][molCnt].nSites,"[sites]",&lineCnt);
                allocSites(&mols[0][molCnt].sites,mols[0][molCnt].nSites);
		for(j=0;j<mols[0][molCnt].nSites;j++)
		{
			getAtomLineFromTOP(fn,in,mols[0][molCnt].sites[j].siteName,&mols[0][molCnt].sites[j].mass,&mols[0][molCnt].sites[j].charge,&mols[0][molCnt].sites[j].pol,&lineCnt);
			mols[0][molCnt].mass+=mols[0][molCnt].sites[j].mass;
		}
		
		getIntFromTOP(fn,in,&mols[0][molCnt].nBonds,"[bonds]",&lineCnt);
		allocBonds(&mols[0][molCnt].bonds,mols[0][molCnt].nBonds);
		for(j=0;j<mols[0][molCnt].nBonds;j++) getBondLineFromTOP(fn,in,&mols[0][molCnt].bonds[j],min,max,&lineCnt,bondOffset);

		getIntFromTOP(fn,in,&mols[0][molCnt].nAngles,"[angles]",&lineCnt);
		allocAngles(&mols[0][molCnt].angles,mols[0][molCnt].nAngles);
                for(j=0;j<mols[0][molCnt].nAngles;j++) getAngleLineFromTOP(fn,in,&mols[0][molCnt].angles[j],min,max,&lineCnt,bondOffset);
		
		getIntFromTOP(fn,in,&mols[0][molCnt].nDihed,"[dihedrals]",&lineCnt);
		allocDihedrals(&mols[0][molCnt].dihed,mols[0][molCnt].nDihed);
		for(j=0;j<mols[0][molCnt].nDihed;j++) getDihedralLineFromTOP(fn,in,&mols[0][molCnt].dihed[j],min,max,&lineCnt,bondOffset);	

		getIntFromTOP(fn,in,&mols[0][molCnt].nImprop,"[impropers]",&lineCnt);
		allocImpropers(&mols[0][molCnt].improp,mols[0][molCnt].nImprop);
                for(j=0;j<mols[0][molCnt].nImprop;j++) getImproperLineFromTOP(fn,in,&mols[0][molCnt].improp[j],min,max,&lineCnt,bondOffset);

		k=molCnt;
		for(copy=1;copy<ncopy;copy++) { 
			cpyMolTop(mols[0][k],&mols[0][k+copy],atoms[0],copy);
			atomCnt+=mols[0][k].nAtoms;
			molCnt++;
		}
		molCnt++;
	}

	j=0;
	for(i=0;i<nMols[0];i++) j+=mols[0][i].nAtoms;
	if(j!=nAtoms[0]) fatal("total atom number inconsistent!\n");
	
	j=0;
        for(i=0;i<nMols[0];i++) j+=mols[0][i].nSites;
        if(j!=nSites[0]) fatal("total site number inconsistent!\n");

	strcpy(mols[0][nMols[0]].molName,"NONE");
	strcpy(mols[0][nMols[0]].chain,"NONE");
	mols[0][nMols[0]].nAtoms=0;
	mols[0][nMols[0]].nSites=0;
	getIntFromTOP(fn,in,&mols[0][nMols[0]].nBonds,"[inter_bonds]",&lineCnt);
	allocBonds(&mols[0][nMols[0]].bonds,mols[0][nMols[0]].nBonds);
	for(j=0;j<mols[0][nMols[0]].nBonds;j++) getBondLineFromTOP(fn,in,&mols[0][nMols[0]].bonds[j],0,nAtoms[0]-1,&lineCnt,0);

	getIntFromTOP(fn,in,&mols[0][nMols[0]].nAngles,"[inter_angles]",&lineCnt);
	allocAngles(&mols[0][nMols[0]].angles,mols[0][nMols[0]].nAngles);
	for(j=0;j<mols[0][nMols[0]].nAngles;j++) getAngleLineFromTOP(fn,in,&mols[0][nMols[0]].angles[j],0,nAtoms[0]-1,&lineCnt,0);

	getIntFromTOP(fn,in,&mols[0][nMols[0]].nDihed,"[inter_dihedrals]",&lineCnt);
	allocDihedrals(&mols[0][nMols[0]].dihed,mols[0][nMols[0]].nDihed);
	for(j=0;j<mols[0][nMols[0]].nDihed;j++) getDihedralLineFromTOP(fn,in,&mols[0][nMols[0]].dihed[j],0,nAtoms[0]-1,&lineCnt,0);

	getIntFromTOP(fn,in,&mols[0][nMols[0]].nImprop,"[inter_impropers]",&lineCnt);
	allocImpropers(&mols[0][nMols[0]].improp,mols[0][nMols[0]].nImprop);
	for(j=0;j<mols[0][nMols[0]].nImprop;j++) getImproperLineFromTOP(fn,in,&mols[0][nMols[0]].improp[j],0,nAtoms[0]-1,&lineCnt,0);

	fclose(in);

	return 0;
}


	
	
	
