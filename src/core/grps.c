#include <math.h>
#include "../../include/dataTypes.h"
#include "../../include/geo.h"
#include "../../include/fatal.h"
#include "../../include/grps.h"

int getGrpCOM(t_atom *atoms,t_grp *grp)
{
        int i;
        real totMass;
        t_vecR tmp;
        real max=-1.0;
        real dist;

        totMass=0.0;

        grp[0].COM.x=0.0; grp[0].COM.y=0.0; grp[0].COM.z=0.0;
        for(i=0;i<grp[0].nAtoms;i++)
        {
                vecRAddScaled(grp[0].COM,atoms[grp[0].atoms[i]].mass,atoms[grp[0].atoms[i]].crd,&tmp);
                vecRCpy(tmp,&grp[0].COM);
                totMass+=atoms[grp[0].atoms[i]].mass;
        }

        if(totMass!=0.0)
        {
                vecRScal(1.0/totMass,grp[0].COM,&tmp);
                vecRCpy(tmp,&grp[0].COM);
        }

        for(i=0;i<grp[0].nAtoms;i++)
        {
                distRSq(atoms[grp[0].atoms[i]].crd,grp[0].COM,&dist);
                if(dist>max) max=dist;
        }
        grp[0].spread=sqrt(max);

        return 0;
}

int getGrpCOMvel(t_atom *atoms,t_grp *grp)
{
	int i,j;
	real m,totMass;

	totMass=0.0;

	grp[0].COMvel.x=0.0; grp[0].COMvel.y=0.0; grp[0].COMvel.z=0.0;
	for(i=0;i<grp[0].nAtoms;i++) {
                j=grp[0].atoms[i];
                m=atoms[j].mass;
                totMass+=m;
                grp[0].COMvel.x+=m*atoms[j].vel.x;
                grp[0].COMvel.y+=m*atoms[j].vel.y;
                grp[0].COMvel.z+=m*atoms[j].vel.z;
        }
	if(totMass!=0.0) {
        	grp[0].COMvel.x/=totMass;
        	grp[0].COMvel.y/=totMass;
        	grp[0].COMvel.z/=totMass;
	}
	return 0;
}

int getGrpInertiaTensor(t_atom *atoms,t_grp *grp)
{
/* expects updated COM of group */
	int i,j;
	real m;
	t_vecR COM,dR;

	grp[0].T.c1.x=0.0;
	grp[0].T.c1.y=0.0;
	grp[0].T.c1.z=0.0;
	grp[0].T.c2.x=0.0;
        grp[0].T.c2.y=0.0;
        grp[0].T.c2.z=0.0;
	grp[0].T.c3.x=0.0;
        grp[0].T.c3.y=0.0;
        grp[0].T.c3.z=0.0;

	for(i=0;i<grp[0].nAtoms;i++) {
		j=grp[0].atoms[i];
		vecRSub(grp[0].COM,atoms[j].crd,&dR);
		m=atoms[j].mass;
		grp[0].T.c1.x+=m*(dR.y*dR.y+dR.z*dR.z);
		grp[0].T.c1.y+=m*(-1.0*dR.x*dR.y);
		grp[0].T.c1.z+=m*(-1.0*dR.x*dR.z);
		grp[0].T.c2.x+=m*(-1.0*dR.y*dR.x);
		grp[0].T.c2.y+=m*(dR.x*dR.x+dR.z*dR.z);
		grp[0].T.c2.z+=m*(-1.0*dR.y*dR.z);
		grp[0].T.c3.x+=m*(-1.0*dR.z*dR.x);
		grp[0].T.c3.y+=m*(-1.0*dR.z*dR.y);
		grp[0].T.c3.z+=m*(dR.x*dR.x+dR.y*dR.y);
	}
	return 0;
}

int getGrpAngularMomentum(t_atom *atoms,t_grp *grp)
{
/* expects updated COM of group */
/* expects updated COMvel of group */
	t_vecR dR,dV,tmp,tmp2;
	real totMass=0.0;
	real m;
	int i,j;

	grp[0].L.x=0.0;
	grp[0].L.y=0.0;
	grp[0].L.z=0.0;

	for(i=0;i<grp[0].nAtoms;i++) {
                j=grp[0].atoms[i];
                m=atoms[j].mass;
		vecRSub(grp[0].COM,atoms[j].crd,&dR);
		vecRSub(grp[0].COMvel,atoms[j].vel,&dV);
		vecRScal(m,dV,&tmp);
		vecRCross(dR,tmp,&tmp2);
		grp[0].L.x+=tmp2.x;
		grp[0].L.y+=tmp2.y;
		grp[0].L.z+=tmp2.z;
	}
	return 0;
}

int removeGrpCOMveloAndAngularMomentum(t_atom *atoms,t_grp *grp)
{
	t_matR invT;
	t_vecR omega;
	real omegaNorm;
	t_vecR omegaDir;
	int i,j;
	t_vecR dR,dRperp,dV,subvelo,tmp;
	real proj;
	
	getGrpCOM(atoms,grp);
	getGrpCOMvel(atoms,grp);
	getGrpInertiaTensor(atoms,grp);
	invMatR(grp[0].T.c1,grp[0].T.c2,grp[0].T.c3,&invT.c1,&invT.c2,&invT.c3);
	getGrpAngularMomentum(atoms,grp);
	colMatDotVecR(invT.c1,invT.c2,invT.c3,grp[0].L,&omega);
	vecRNorm(omega,&omegaNorm);
	vecRScal(1/omegaNorm,omega,&omegaDir);

	for(i=0;i<grp[0].nAtoms;i++) {
		j=grp[0].atoms[i];
		vecRSub(grp[0].COM,atoms[j].crd,&dR);
		vecRProd(dR,omegaDir,&proj);
		proj*=-1.0;
		vecRAddScaled(dR,proj,omegaDir,&dRperp);
		vecRCross(omega,dRperp,&subvelo);
		vecRSub(atoms[j].vel,grp[0].COMvel,&tmp);
		vecRSub(tmp,subvelo,&atoms[j].vel);
	}
	return 0;
}

		
