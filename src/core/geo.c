#include <stdio.h>
#include <math.h>
#include "../../include/dataTypes.h"
#include "../../include/fatal.h"

int vecRCpy(t_vecR r1,t_vecR *r2)
{
        r2->x=r1.x;
        r2->y=r1.y;
        r2->z=r1.z;
        return 0;
}

int vecRAdd(t_vecR r1,t_vecR r2,t_vecR *r3)
{
        r3->x=r1.x+r2.x;
        r3->y=r1.y+r2.y;
        r3->z=r1.z+r2.z;
        return 0;
}

int vecRSub(t_vecR r1,t_vecR r2,t_vecR *r3)
{
        r3->x=r2.x-r1.x;
        r3->y=r2.y-r1.y;
        r3->z=r2.z-r1.z;
        return 0;
}

int vecRMulti(t_vecR r1,t_vecR r2,t_vecR *r3)
{
        r3->x=r1.x*r2.x;
        r3->y=r1.y*r2.y;
        r3->z=r1.z*r2.z;
        return 0;
}

int vecRProd(t_vecR r1,t_vecR r2,real *p)
{
        p[0]=r1.x*r2.x+r1.y*r2.y+r1.z*r2.z;
        return 0;
}

int vecRScal(real s,t_vecR r1,t_vecR *r2)
{
        r2->x=s*r1.x;
        r2->y=s*r1.y;
        r2->z=s*r1.z;
        return 0;
}

int vecRAddScaled(t_vecR r1,real s,t_vecR r2,t_vecR *r3)
{
        r3[0].x=r1.x+s*r2.x;
        r3[0].y=r1.y+s*r2.y;
        r3[0].z=r1.z+s*r2.z;
        return 0;
}

int vecRCross(t_vecR r1,t_vecR r2,t_vecR *r3)
{
        r3[0].x=r1.y*r2.z-r1.z*r2.y;
        r3[0].y=r1.z*r2.x-r1.x*r2.z;
        r3[0].z=r1.x*r2.y-r1.y*r2.x;
        return 0;
}

int vecRNormSq(t_vecR r1,real *dSq)
{
        dSq[0]=r1.x*r1.x+r1.y*r1.y+r1.z*r1.z;
        return 0;
}

int vecRNorm(t_vecR r1,real *d)
{
        d[0]=(real)sqrt((double)(r1.x*r1.x+r1.y*r1.y+r1.z*r1.z));
        return 0;
}

int distRSq(t_vecR r1,t_vecR r2,real *dSq)
{
        t_vecR link;

        vecRSub(r1,r2,&link);
        vecRNormSq(link,dSq);
        return 0;
}

int distR(t_vecR r1,t_vecR r2,real *d)
{
        t_vecR link;

        vecRSub(r1,r2,&link);
        vecRNorm(link,d);
        return 0;
}

/* if input vectors describe columns => output vectors will be the colums of the inverse matrix */
/* if input vectors describe rows    => output vectors will be the rows of the inverse matrix */
int invMatR(t_vecR c1,t_vecR c2,t_vecR c3,t_vecR *cInv1,t_vecR *cInv2,t_vecR *cInv3)
{
	real n;

	n=0.0;
	n+=c1.x*c2.y*c3.z;
	n+=c1.y*c2.z*c3.x;
	n+=c1.z*c2.x*c3.y;
	n-=c1.z*c2.y*c3.x;
	n-=c1.x*c2.z*c3.y;
	n-=c1.y*c2.x*c3.z;

	cInv1[0].x=(c2.y*c3.z-c2.z*c3.y)/n;
	cInv1[0].y=(c1.z*c3.y-c1.y*c3.z)/n;
	cInv1[0].z=(c1.y*c2.z-c1.z*c2.y)/n;
	cInv2[0].x=(c2.z*c3.x-c2.x*c3.z)/n;
	cInv2[0].y=(c1.x*c3.z-c1.z*c3.x)/n;
	cInv2[0].z=(c1.z*c2.x-c1.x*c2.z)/n;
	cInv3[0].x=(c2.x*c3.y-c2.y*c3.x)/n;
	cInv3[0].y=(c1.y*c3.x-c1.x*c3.y)/n;
	cInv3[0].z=(c1.x*c2.y-c1.y*c2.x)/n;

	return 0;
}

int colMatDotVecR(t_vecR c1,t_vecR c2,t_vecR c3,t_vecR r,t_vecR *res)
{
	res[0].x=c1.x*r.x+c2.x*r.y+c3.x*r.z;
	res[0].y=c1.y*r.x+c2.y*r.y+c3.y*r.z;
	res[0].z=c1.z*r.x+c2.z*r.y+c3.z*r.z;
	return 0;
}

int rowMatDotVecR(t_vecR r1,t_vecR r2,t_vecR r3,t_vecR r,t_vecR *res)
{
	res[0].x=r1.x*r.x+r1.y*r.y+r1.z*r.z;
	res[0].y=r2.x*r.x+r2.y*r.y+r2.z*r.z;
	res[0].z=r3.x*r.x+r3.y*r.y+r3.z*r.z;
	return 0;
}

int makePBCtransMat(t_box *box)
{
	invMatR(box[0].a,box[0].b,box[0].c,&box[0].trans.c1,&box[0].trans.c2,&box[0].trans.c3);
	return 0;
}

int transXYZToPBCbasis(t_box box,t_vecR r,t_vecR *res)
{
	colMatDotVecR(box.trans.c1,box.trans.c2,box.trans.c3,r,res);
	return 0;
}

int transPBCToXYZbasis(t_box box,t_vecR r,t_vecR *res)
{
        colMatDotVecR(box.a,box.b,box.c,r,res);
        return 0;
}

int wrapTrajectoryAtoms(t_atom *atoms,int nAtoms,t_box box)
{
	int i;
	t_vecR tmp;
	real shift;


	for(i=0;i<nAtoms;i++)
	{
		transXYZToPBCbasis(box,atoms[i].crd,&tmp);
		if(tmp.x>0.5 || tmp.x<-0.5)
		{
			shift=tmp.x+0.5;
			shift=floor(shift);
			tmp.x-=shift;
		}
		if(tmp.y>0.5 || tmp.y<-0.5)
                {
                        shift=tmp.y+0.5;
                        shift=floor(shift);
                        tmp.y-=shift;
                }
		if(tmp.z>0.5 || tmp.z<-0.5)
                {
                        shift=tmp.z+0.5;
                        shift=floor(shift);
                        tmp.z-=shift;
                }
		transPBCToXYZbasis(box,tmp,&atoms[i].crd);
	}
	return 0;
}

int wrapTrajectoryAtomsGrp(t_atom *atoms,t_grp grp,t_box box)
{
        int i;
        t_vecR tmp;
        real shift;


        for(i=0;i<grp.nAtoms;i++)
        {
                transXYZToPBCbasis(box,atoms[grp.atoms[i]].crd,&tmp);
                if(tmp.x>0.5 || tmp.x<-0.5)
                {
                        shift=tmp.x+0.5;
                        shift=floor(shift);
                        tmp.x-=shift;
                }
                if(tmp.y>0.5 || tmp.y<-0.5)
                {
                        shift=tmp.y+0.5;
                        shift=floor(shift);
                        tmp.y-=shift;
                }
                if(tmp.z>0.5 || tmp.z<-0.5)
                {
                        shift=tmp.z+0.5;
                        shift=floor(shift);
                        tmp.z-=shift;
                }
                transPBCToXYZbasis(box,tmp,&atoms[grp.atoms[i]].crd);
        }
        return 0;
}

int wrapTrajectoryMol0(t_atom *atoms,t_mol *mols,int nMols,t_box box)
{
	int i,j;
	t_vecR tmp;
	t_vecR shift;

	for(i=0;i<nMols;i++)
	{
		transXYZToPBCbasis(box,atoms[mols[i].atoms[0]].crd,&tmp);
                if(tmp.x>0.5 || tmp.x<-0.5)
                {
                        shift.x=tmp.x+0.5;
                        shift.x=floor(shift.x);
                } else shift.x=0.0;
                if(tmp.y>0.5 || tmp.y<-0.5)
                {
                        shift.y=tmp.y+0.5;
                        shift.y=floor(shift.y);
                } else shift.y=0.0;
                if(tmp.z>0.5 || tmp.z<-0.5)
                {
                        shift.z=tmp.z+0.5;
                        shift.z=floor(shift.z);
                } else shift.z=0.0;
                transPBCToXYZbasis(box,shift,&tmp);

		for(j=0;j<mols[i].nAtoms;j++)
		{
			vecRSub(tmp,atoms[mols[i].atoms[j]].crd,&shift);
			vecRCpy(shift,&atoms[mols[i].atoms[j]].crd);
		}
		for(j=0;j<mols[i].nSites;j++)
                {
                        vecRSub(tmp,mols[i].sites[j].crd,&shift);
                        vecRCpy(shift,&mols[i].sites[j].crd);
                }
	}
	return 0;
}

int wrapTrajectoryMol0Grp(t_atom *atoms,t_mol *mols,t_grp grp,t_box box)
{
        int i,j;
        t_vecR tmp;
        t_vecR shift;

        for(i=0;i<grp.nMols;i++)
        {
                transXYZToPBCbasis(box,atoms[mols[grp.mols[i]].atoms[0]].crd,&tmp);
                if(tmp.x>0.5 || tmp.x<-0.5)
                {
                        shift.x=tmp.x+0.5;
                        shift.x=floor(shift.x);
                } else shift.x=0.0;
                if(tmp.y>0.5 || tmp.y<-0.5)
                {
                        shift.y=tmp.y+0.5;
                        shift.y=floor(shift.y);
                } else shift.y=0.0;
                if(tmp.z>0.5 || tmp.z<-0.5)
                {
                        shift.z=tmp.z+0.5;
                        shift.z=floor(shift.z);
                } else shift.z=0.0;
                transPBCToXYZbasis(box,shift,&tmp);

                for(j=0;j<mols[grp.mols[i]].nAtoms;j++)
                {
                        vecRSub(tmp,atoms[mols[grp.mols[i]].atoms[j]].crd,&shift);
                        vecRCpy(shift,&atoms[mols[grp.mols[i]].atoms[j]].crd);
                }
                for(j=0;j<mols[grp.mols[i]].nSites;j++)
                {
                        vecRSub(tmp,mols[grp.mols[i]].sites[j].crd,&shift);
                        vecRCpy(shift,&mols[grp.mols[i]].sites[j].crd);
                }
        }
        return 0;
}

/* COM expected to be known */
int wrapTrajectoryMolCOM(t_atom *atoms,t_mol *mols,int nMols,t_box box)
{
        int i,j;
        t_vecR tmp;
        t_vecR shift;

        for(i=0;i<nMols;i++)
        {
                transXYZToPBCbasis(box,mols[i].COM,&tmp);
                if(tmp.x>0.5 || tmp.x<-0.5)
                {
                        shift.x=tmp.x+0.5;
                        shift.x=floor(shift.x);
                } else shift.x=0.0;
                if(tmp.y>0.5 || tmp.y<-0.5)
                {
                        shift.y=tmp.y+0.5;
                        shift.y=floor(shift.y);
                } else shift.y=0.0;
                if(tmp.z>0.5 || tmp.z<-0.5)
                {
                        shift.z=tmp.z+0.5;
                        shift.z=floor(shift.z);
                } else shift.z=0.0;
                transPBCToXYZbasis(box,shift,&tmp);

                for(j=0;j<mols[i].nAtoms;j++)
                {
                        vecRSub(tmp,atoms[mols[i].atoms[j]].crd,&shift);
                        vecRCpy(shift,&atoms[mols[i].atoms[j]].crd);
                }
		for(j=0;j<mols[i].nSites;j++)
                {
                        vecRSub(tmp,mols[i].sites[j].crd,&shift);
                        vecRCpy(shift,&mols[i].sites[j].crd);
                }
		vecRSub(tmp,mols[i].COM,&shift);
		vecRCpy(shift,&mols[i].COM);
        }
        return 0;
}

/* COM expected to be known */
int wrapTrajectoryMolCOMGrp(t_atom *atoms,t_mol *mols,t_grp grp,t_box box)
{
        int i,j;
        t_vecR tmp;
        t_vecR shift;

        for(i=0;i<grp.nMols;i++)
        {
                transXYZToPBCbasis(box,mols[grp.mols[i]].COM,&tmp);
                if(tmp.x>0.5 || tmp.x<-0.5)
                {
                        shift.x=tmp.x+0.5;
                        shift.x=floor(shift.x);
                } else shift.x=0.0;
                if(tmp.y>0.5 || tmp.y<-0.5)
                {
                        shift.y=tmp.y+0.5;
                        shift.y=floor(shift.y);
                } else shift.y=0.0;
                if(tmp.z>0.5 || tmp.z<-0.5)
                {
                        shift.z=tmp.z+0.5;
                        shift.z=floor(shift.z);
                } else shift.z=0.0;
                transPBCToXYZbasis(box,shift,&tmp);

                for(j=0;j<mols[grp.mols[i]].nAtoms;j++)
                {
                        vecRSub(tmp,atoms[mols[grp.mols[i]].atoms[j]].crd,&shift);
                        vecRCpy(shift,&atoms[mols[grp.mols[i]].atoms[j]].crd);
                }
                for(j=0;j<mols[grp.mols[i]].nSites;j++)
                {
                        vecRSub(tmp,mols[grp.mols[i]].sites[j].crd,&shift);
                        vecRCpy(shift,&mols[grp.mols[i]].sites[j].crd);
                }
                vecRSub(tmp,mols[grp.mols[i]].COM,&shift);
                vecRCpy(shift,&mols[grp.mols[i]].COM);
        }
        return 0;
}

int linkPBCfull(t_vecR *link,t_box box)
{
	t_vecR trans;

	transXYZToPBCbasis(box,link[0],&trans);
	while(trans.x>0.5) trans.x-=1;
	while(trans.x<-0.5) trans.x+=1;
	while(trans.y>0.5) trans.y-=1;
	while(trans.y<-0.5) trans.y+=1;
	while(trans.z>0.5) trans.z-=1;
	while(trans.z<-0.5) trans.z+=1;
	transPBCToXYZbasis(box,trans,link);

	return 0;
}

int linkPBCfullRaw(t_vecR *link,t_box box)
{
        t_vecR trial,best;
        int a,b,c;
        real cur,min;

        min=9999.9;

        for(a=-1;a<=1;a++)
        {
                for(b=-1;b<=1;b++)
                {
                        for(c=-1;c<=1;c++)
                        {
                                vecRAddScaled(link[0],(real)a,box.a,&trial);
                                vecRAddScaled(link[0],(real)b,box.b,&trial);
                                vecRAddScaled(link[0],(real)c,box.c,&trial);
                                vecRNormSq(trial,&cur);
                                if(cur<min)
                                {
                                        min=cur;
                                        vecRCpy(trial,&best);
                                }
                        }
                }
        }
        vecRCpy(best,link);
        return 0;
}

int linkPBCortho(t_vecR *link,t_box box)
{
        if(link[0].x>0.5*box.a.x) link[0].x-=box.a.x;
        else if (link[0].x<-0.5*box.a.x) link[0].x+=box.a.x;
 
        if(link[0].y>0.5*box.b.y) link[0].y-=box.b.y;
        else if (link[0].y<-0.5*box.b.y) link[0].y+=box.b.y;

        if(link[0].z>0.5*box.c.z) link[0].z-=box.c.z;
        else if (link[0].z<-0.5*box.c.z) link[0].z+=box.c.z;

        return 0;
}

int linkPBC(t_vecR r1,t_vecR r2,t_box box,t_vecR *res)
{
        t_vecR link;

        vecRSub(r1,r2,&link);
 /*       if(box.ortho==1)
        {
		linkPBCortho(&link,box);
	} else {
		linkPBCfull(&link,box);
	}*/
	linkPBCfull(&link,box);
        vecRCpy(link,res);
        return 0;
}

int fixMols(t_atom *atoms,t_mol *mols,int nMols,t_box box)
{
        int i,j;
        t_vecR tmp,tmp2;

        for(i=0;i<nMols;i++)
        {
                for(j=1;j<mols[i].nAtoms;j++)
                {
                        linkPBC(atoms[mols[i].atoms[0]].crd,atoms[mols[i].atoms[j]].crd,box,&tmp);
                        vecRAdd(atoms[mols[i].atoms[0]].crd,tmp,&tmp2);
			vecRCpy(tmp2,&atoms[mols[i].atoms[j]].crd);
		}
	}
	return 0;
}

int distSqPBC(t_vecR r1,t_vecR r2,t_box box,real *dSq)
{
        t_vecR link;

	linkPBC(r1,r2,box,&link);
        vecRNormSq(link,dSq);
        return 0;
}

int distPBC(t_vecR r1,t_vecR r2,t_box box,real *d)
{
	real dSq;
	distSqPBC(r1,r2,box,&dSq);
	d[0]=sqrt(dSq);
	return 0;
}

int distDeriv(t_vecR a,t_vecR da,t_vecR b,t_vecR db,t_box box,real *d)
{
        t_vecR link,vd;
        real tmp,tmp2;

        linkPBC(a,b,box,&link);
        vecRNorm(link,&tmp);
        vecRSub(da,db,&vd);
        vecRProd(link,vd,&tmp2);
        d[0]=tmp2/tmp;
        return 0;
}

int getAnglePBC(t_vecR a,t_vecR b,t_vecR center,t_box box,real *angle)
{
        t_vecR link1,link2;
        real dotproduct;
        real len1,len2;
        real arg;
	char err[100];

        linkPBC(center,a,box,&link1);
        linkPBC(center,b,box,&link2);
        vecRNorm(link1,&len1);
        vecRNorm(link2,&len2);
        vecRProd(link1,link2,&dotproduct);
        arg=dotproduct/(len1*len2);
        if(arg<-1.1 || arg>1.1) {
		sprintf(err,"acos argument out of range {-1,1}: %f -> getAngle \n",arg);
		fatal(err);
	}
        if(arg<-1.0) angle[0]=3.14159; else if(arg>1.0) angle[0]=0.0; else angle[0]=acos(arg);
        return 0;
}

int angleDeriv(t_vecR a,t_vecR da,t_vecR b,t_vecR db,t_vecR c,t_vecR dc,t_box box,real *d)
{
        t_vecR d1,d2,v1,v2;
        real d1d2,d1v1,d1v2,d2v1,d2v2,d1Sq,d2Sq,d1Tr,d2Tr,d1len,d2len,tmp,tmp2;

        linkPBC(c,a,box,&d1);
        linkPBC(c,b,box,&d2);
        vecRSub(dc,da,&v1);
        vecRSub(dc,db,&v2);

        vecRProd(d1,d2,&d1d2);
        vecRProd(d1,v1,&d1v1);
        vecRProd(d1,v2,&d1v2);
        vecRProd(d2,v1,&d2v1);
        vecRProd(d2,v2,&d2v2);
        vecRNormSq(d1,&d1Sq);
        vecRNormSq(d2,&d2Sq);
        vecRNorm(d1,&d1len); d1Tr=d1len*d1len*d1len;
        vecRNorm(d2,&d2len); d2Tr=d2len*d2len*d2len;

        tmp=(d1Sq*d2Sq*(d2v1+d1v2)-d1d2*(d2Sq*d1v1+d1Sq*d2v2))/(d1Tr*d2Tr);
        tmp2=d1d2/(d1len*d2len);
	tmp2=tmp2*tmp2;
	if(tmp2>=1.0) {
		printf("WARNING: trying to evaluate 1/[sqrt(1-x)] with x = %f\n",tmp2);
                printf("  -> angleDeriv\n  ->  geometry.c\n");
	}
        tmp2=-1/(sqrt(1-(tmp2)));
        d[0]=tmp*tmp2;
        return 0;
}

