#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/fatal.h"
#include "../include/geo.h"

int assignWannierByDist(t_vecR *siteCrds,int nSites,t_atom *atoms,t_mol *mols,int nMols,t_box box,real cut)
{
	int i,j,k,l;
	real distSq;
	int *tags;
	char err[200];
	t_vecR tmp;

	allocInts(&tags,nSites);
	for(i=0;i<nSites;i++) tags[i]=0;

	for(i=0;i<nMols;i++)
	{
		k=0;
		for(j=0;j<mols[i].nAtoms;j++)
		{
			if(atoms[mols[i].atoms[j]].mass>3.0)
			{
				for(l=0;l<nSites;l++)
				{
					if(tags[l]==0)
					{
						distSqPBC(atoms[mols[i].atoms[j]].crd,siteCrds[l],box,&distSq);
						if(distSq<=cut)
						{
							linkPBC(atoms[mols[i].atoms[j]].crd,siteCrds[l],box,&tmp);
							vecRAdd(atoms[mols[i].atoms[j]].crd,tmp,&mols[i].sites[k].crd);
							tags[l]=1;
							k++;
						}
					}
				}
			}
		}
		if(k!=mols[i].nSites)
		{
			sprintf(err,"NOT ALL WANNIER CENTERS WERE FOUND FOR MOLECULE %d WITHIN DISTANCE OF %f: found %d of %d\n -> assignWannierByDist\n -> assignWC.c\n",i,sqrt(cut),k,mols[i].nSites);
			fatal(err);
		}
	}
	for(i=0;i<nSites;i++)
	{
		if(tags[i]!=1)
		{
			sprintf(err,"NOT ALL WANNIER CENTERS WERE ASSIGNED TO A MOLECULE\n -> assignWannierByDist\n -> assignWC.c\n");
			fatal(err);
		}
	}
	free(tags);
	return 0;
}
