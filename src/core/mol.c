#include "../../include/dataTypes.h"
#include "../../include/geo.h"
#include "../../include/fatal.h"

int getMolCOM(t_atom *atoms,t_mol *mols,int nMols,int needVelo)
{
	int i,j;
	t_vecR tmp;

	for(i=0;i<nMols;i++)
	{
		mols[i].COM.x=0.0; mols[i].COM.y=0.0; mols[i].COM.z=0.0;
		for(j=0;j<mols[i].nAtoms;j++)
		{
			vecRAddScaled(mols[i].COM,atoms[mols[i].atoms[j]].mass,atoms[mols[i].atoms[j]].crd,&tmp);
                	vecRCpy(tmp,&mols[i].COM);
		}

		for(j=0;j<mols[i].nSites;j++)
                {
                        vecRAddScaled(mols[i].COM,mols[i].sites[j].mass,mols[i].sites[j].crd,&tmp);
                        vecRCpy(tmp,&mols[i].COM);
                }

		if(mols[i].mass!=0.0)
        	{
        	        vecRScal(1.0/mols[i].mass,mols[i].COM,&tmp);
        	        vecRCpy(tmp,&mols[i].COM);
        	}

		mols[i].dCOM.x=0.0; mols[i].dCOM.y=0.0; mols[i].dCOM.z=0.0;
		if(needVelo==1)
		{
			for(j=0;j<mols[i].nAtoms;j++)
 	        	{
                        	vecRAddScaled(mols[i].dCOM,atoms[mols[i].atoms[j]].mass,atoms[mols[i].atoms[j]].vel,&tmp);
                        	vecRCpy(tmp,&mols[i].dCOM);
			}

			for(j=0;j<mols[i].nSites;j++)
			{
				if(mols[i].sites[j].mass!=0.0) fatal("NO VELOCITIES FOR SITES: UNABLE TO COMPUTE dCOM!\n -> getMolCOM\n -> mol.c\n");
			}

			if(mols[i].mass!=0.0)
                	{
                        	vecRScal(1.0/mols[i].mass,mols[i].dCOM,&tmp);
                        	vecRCpy(tmp,&mols[i].dCOM);
                	}
                }
	}
	return 0;
}

int getMolDM(t_atom *atoms,t_mol *mols,int nMols,int needVelo)
{
        int i,j;
        t_vecR tmp;

        for(i=0;i<nMols;i++)
        {
                mols[i].DM.x=0.0; mols[i].DM.y=0.0; mols[i].DM.z=0.0;
                for(j=0;j<mols[i].nAtoms;j++)
                {
                        vecRAddScaled(mols[i].DM,atoms[mols[i].atoms[j]].charge,atoms[mols[i].atoms[j]].crd,&tmp);
                        vecRCpy(tmp,&mols[i].DM);
                }

                for(j=0;j<mols[i].nSites;j++)
                {
                        vecRAddScaled(mols[i].DM,mols[i].sites[j].charge,mols[i].sites[j].crd,&tmp);
                        vecRCpy(tmp,&mols[i].DM);
                }

                mols[i].dDM.x=0.0; mols[i].dDM.y=0.0; mols[i].dDM.z=0.0;
                if(needVelo==1)
                {
                        for(j=0;j<mols[i].nAtoms;j++)
                        {
                                vecRAddScaled(mols[i].dDM,atoms[mols[i].atoms[j]].charge,atoms[mols[i].atoms[j]].vel,&tmp);
                                vecRCpy(tmp,&mols[i].dDM);
                        }

			for(j=0;j<mols[i].nSites;j++)
			{
				if(mols[i].sites[j].charge!=0.0) printf("NO VELOCITIES FOR SITES: UNABLE TO COMPUTE dDM!\n -> getMolDM\n -> mol.c\n");
			}
                }
        }
        return 0;
}

