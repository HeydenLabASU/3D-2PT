#include <stdio.h>
#include <stdlib.h>
#include "../../include/dataTypes.h"
#include "../../include/fatal.h"

void *save_calloc(char *name,char *file,int line,
                  unsigned nelem,unsigned elsize)
{
  void *p;
  char err[300];

  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL)
      {
         sprintf(err,"calloc for %s (nelem=%d, elsize=%d, file %s line %d)\n",name,nelem,elsize,file,line);
         fatal(err);
      }
    }
  return p;
}

void save_free(char *name,char *file,int line, void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}

void *save_malloc(size_t size)
{
        void *ptr;
        char buffer[200];

        if((ptr=malloc(size))==NULL)
        {
                sprintf(buffer,"MEMORY ALLOCATION FAILED!!!\ntried to allocate %lu bytes\n -> save_malloc\n -> alloc.c\n",size);
                fatal(buffer);
        }

        return ptr;
}

void *save_realloc(void *ptr,size_t size)
{
        char buffer[200];

        if((ptr=realloc(ptr,size))==NULL && size!=0)
        {
                sprintf(buffer,"MEMORY REALLOCATION FAILED!!!\ntried to reallocate %lu bytes\n -> save_realloc\n -> alloc.c\n",size);
                fatal(buffer);
        }

        return ptr;
}

int allocMols(t_mol **mols,int n)
{
	mols[0]=(t_mol*)save_malloc(n*sizeof(t_mol));
	return 0;
}

int allocAtoms(t_atom **atoms,int n)
{
        atoms[0]=(t_atom*)save_malloc(n*sizeof(t_atom));
        return 0;
}

int allocSites(t_site **sites,int n)
{
        sites[0]=(t_site*)save_malloc(n*sizeof(t_site));
        return 0;
}

int allocBonds(t_bond **bonds,int n)
{
        bonds[0]=(t_bond*)save_malloc(n*sizeof(t_bond));
        return 0;
}

int allocAngles(t_angle **angles,int n)
{
        angles[0]=(t_angle*)save_malloc(n*sizeof(t_angle));
        return 0;
}

int allocDihedrals(t_dihed **dihedrals,int n)
{
        dihedrals[0]=(t_dihed*)save_malloc(n*sizeof(t_dihed));
        return 0;
}

int allocImpropers(t_improp **impropers,int n)
{
        impropers[0]=(t_improp*)save_malloc(n*sizeof(t_improp));
        return 0;
}

int allocInts(int **ints,int n)
{
        ints[0]=(int*)save_malloc(n*sizeof(int));
        return 0;
}

int allocReals(real **reals,int n)
{
        reals[0]=(real*)save_malloc(n*sizeof(real));
        return 0;
}

int allocVecRs(t_vecR **vecRs,int n)
{
	vecRs[0]=(t_vecR*)save_malloc(n*sizeof(t_vecR));
	return 0;
}

int allocStatGrpDefs(t_statGrpDef **statGrpDef,int nStat)
{
	statGrpDef[0]=(t_statGrpDef*)save_malloc(nStat*sizeof(t_statGrpDef));
	return 0;
}

int allocDynGrpDefs(t_dynGrpDef **dynGrpDef,int nDyn)
{
	dynGrpDef[0]=(t_dynGrpDef*)save_malloc(nDyn*sizeof(t_dynGrpDef));
	return 0;
}

int allocCombGrpDefs(t_combGrpDef **combGrpDef,int nComb)
{
	combGrpDef[0]=(t_combGrpDef*)save_malloc(nComb*sizeof(t_combGrpDef));
	return 0;
}

