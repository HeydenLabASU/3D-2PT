#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/dataTypes.h"
#include "../../include/fatal.h"
#include "../../include/alloc.h"
#include "../../include/select.h"
#include "../../include/align.h"
#include "../../include/grps.h"
#include "../../include/geo.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);

void jacobi(double **a,int n,double d[],double **v,int *nrot)
{
  int j,i;
  int iq,ip;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  snew(b,n);
  snew(z,n);
  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
 }
  for (ip=0; ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1; i<=50; i++) {
    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      sfree(z);
      sfree(b);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
            && fabs(d[iq])+g == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if (fabs(h)+g == fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0; j<ip; j++) {
            ROTATE(a,j,ip,j,iq)
	  }
          for (j=ip+1; j<iq; j++) {
            ROTATE(a,ip,j,j,iq)
            }
          for (j=iq+1; j<n; j++) {
            ROTATE(a,ip,j,iq,j)
            }
          for (j=0; j<n; j++) {
            ROTATE(v,j,ip,j,iq)
            }
          ++(*nrot);
        }
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] +=  z[ip];
      d[ip]  =  b[ip];
      z[ip]  =  0.0;
    }
  }
  fatal("Error: Too many iterations!!!\n -> jacobi\n -> align.c\n");
}

void calc_fit_R(int natoms,real *w_rls,rvec *xp,rvec *x,matrix R)
{
  int    c,r,n,j,m,i,irot;
  double **omega,**om;
  double d[2*3],xnr,xpc;
  matrix vh,vk,u;
  real   mn;
  int    index;
  real   max_d;

  snew(omega,2*3);
  snew(om,2*3);
  for(i=0; i<2*3; i++) {
    snew(omega[i],2*3);
    snew(om[i],2*3);
  }
  
  for(i=0; i<2*3; i++) {
    d[i]=0;
    for(j=0; j<2*3; j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms);n++)
    if ((mn = w_rls[n]) != 0.0)
      for(c=0; (c<3); c++) {
	xpc=xp[n][c];
	for(r=0; (r<3); r++) {
	  xnr=x[n][r];
	  u[c][r]+=mn*xnr*xpc;
	}
      }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<2*3; r++)
    for(c=0; c<=r; c++)
      if (r>=3 && c<3) {
        omega[r][c]=u[r-3][c];
        omega[c][r]=u[r-3][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }
  
  /*determine h and k*/
  jacobi(omega,2*3,d,om,&irot);
  /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
   *real       **v = v[0..n-1][0..n-1] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */
  
  if (irot==0) printf("IROT=0\n -> jacobi\n -> align.c");
  
  index=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0; j<2; j++) {
    max_d=-1000;
    for(i=0; i<2*3; i++)
      if (d[i]>max_d) {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0; i<3; i++) {
      vh[j][i]=1.41421356237*om[i][index];
      vk[j][i]=1.41421356237*om[i+3][index];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  
  oprod(vh[0],vh[1],vh[2]);
  oprod(vk[0],vk[1],vk[2]);

  /*determine R*/
  for(r=0; r<3; r++)
    for(c=0; c<3; c++)
      R[r][c] = vk[0][r]*vh[0][c] +
	        vk[1][r]*vh[1][c] +
	        vk[2][r]*vh[2][c];

  for(i=0; i<2*3; i++) {
    sfree(omega[i]);
    sfree(om[i]);
  }
  sfree(omega);
  sfree(om);
}

int getAlignRef(t_grp *grp,t_atom *atoms,t_box box,real **w_rls,rvec **refCrd,t_vecR *refCOM,rvec **x)
{
	int i;

	getGrpCOM(atoms,grp);
	vecRCpy(grp[0].COM,refCOM);
	allocReals(w_rls,grp[0].nAtoms);
	refCrd[0]=(rvec*)save_malloc(grp[0].nAtoms*sizeof(rvec));
	for(i=0;i<grp[0].nAtoms;i++)
	{
		w_rls[0][i]=atoms[grp[0].atoms[i]].mass;
		refCrd[0][i][0]=atoms[grp[0].atoms[i]].crd.x-refCOM[0].x;
		refCrd[0][i][1]=atoms[grp[0].atoms[i]].crd.y-refCOM[0].y;
		refCrd[0][i][2]=atoms[grp[0].atoms[i]].crd.z-refCOM[0].z;
	}
	x[0]=(rvec*)save_malloc(grp[0].nAtoms*sizeof(rvec));

	return 0;
}

int alignGroup(t_grp *grp,t_atom *atoms,int nAtoms,t_box *box,real *w_rls,rvec *refCrd,t_vecR refCOM,rvec *x)
{
	int    i,j,c;
	matrix R;
	rvec   x_old;
	t_vecR tmp;

	getGrpCOM(atoms,grp);
	for(i=0;i<nAtoms;i++) {
		vecRSub(grp[0].COM,atoms[i].crd,&tmp);
		vecRCpy(tmp,&atoms[i].crd);
	}

	for(i=0;i<grp[0].nAtoms;i++)
	{
		x[i][0]=atoms[grp[0].atoms[i]].crd.x;
		x[i][1]=atoms[grp[0].atoms[i]].crd.y;
		x[i][2]=atoms[grp[0].atoms[i]].crd.z;
	}

	/* Calculate the rotation matrix R */
	calc_fit_R(grp[0].nAtoms,w_rls,refCrd,x,R);

	/* rotate all coords */
	for(j=0;j<nAtoms;j++) {
		x_old[0]=atoms[j].crd.x;
		x_old[1]=atoms[j].crd.y;
		x_old[2]=atoms[j].crd.z;
		atoms[j].crd.x=0.0;
		for(c=0; c<3; c++) atoms[j].crd.x+=R[0][c]*x_old[c];
		atoms[j].crd.y=0.0;
                for(c=0; c<3; c++) atoms[j].crd.y+=R[1][c]*x_old[c];
		atoms[j].crd.z=0.0;
                for(c=0; c<3; c++) atoms[j].crd.z+=R[2][c]*x_old[c];
		/*vecRAdd(atoms[j].crd,refCOM,&tmp);
		vecRCpy(tmp,&atoms[j].crd);*/
        }
	/* rotate box */
	x_old[0]=box[0].a.x;
        x_old[1]=box[0].a.y;
        x_old[2]=box[0].a.z;
	box[0].a.x=0.0;
	for(c=0; c<3; c++) box[0].a.x+=R[0][c]*x_old[c];
	box[0].a.y=0.0;
        for(c=0; c<3; c++) box[0].a.y+=R[1][c]*x_old[c];
	box[0].a.z=0.0;
        for(c=0; c<3; c++) box[0].a.z+=R[2][c]*x_old[c];
	x_old[0]=box[0].b.x;
        x_old[1]=box[0].b.y;
        x_old[2]=box[0].b.z;
        box[0].b.x=0.0;
        for(c=0; c<3; c++) box[0].b.x+=R[0][c]*x_old[c];
        box[0].b.y=0.0;
        for(c=0; c<3; c++) box[0].b.y+=R[1][c]*x_old[c];
        box[0].b.z=0.0;
        for(c=0; c<3; c++) box[0].b.z+=R[2][c]*x_old[c];
	x_old[0]=box[0].c.x;
        x_old[1]=box[0].c.y;
        x_old[2]=box[0].c.z;
        box[0].c.x=0.0;
        for(c=0; c<3; c++) box[0].c.x+=R[0][c]*x_old[c];
        box[0].c.y=0.0;
        for(c=0; c<3; c++) box[0].c.y+=R[1][c]*x_old[c];
        box[0].c.z=0.0;
        for(c=0; c<3; c++) box[0].c.z+=R[2][c]*x_old[c];

	makePBCtransMat(box);

	return 0;
}

int alignGroupVel(t_grp *grp,t_atom *atoms,int nAtoms,t_box *box,real *w_rls,rvec *refCrd,t_vecR refCOM,rvec *x)
{
        int    i,j,c;
        matrix R;
        rvec   x_old;
        t_vecR tmp;

        getGrpCOM(atoms,grp);
        for(i=0;i<nAtoms;i++) {
                vecRSub(grp[0].COM,atoms[i].crd,&tmp);
                vecRCpy(tmp,&atoms[i].crd);
        }

        for(i=0;i<grp[0].nAtoms;i++)
        {
                x[i][0]=atoms[grp[0].atoms[i]].crd.x;
                x[i][1]=atoms[grp[0].atoms[i]].crd.y;
                x[i][2]=atoms[grp[0].atoms[i]].crd.z;
        }

        /* Calculate the rotation matrix R */
        calc_fit_R(grp[0].nAtoms,w_rls,refCrd,x,R);

        /* rotate all coords */
        for(j=0;j<nAtoms;j++) {
                x_old[0]=atoms[j].crd.x;
                x_old[1]=atoms[j].crd.y;
                x_old[2]=atoms[j].crd.z;
                atoms[j].crd.x=0.0;
                for(c=0; c<3; c++) atoms[j].crd.x+=R[0][c]*x_old[c];
                atoms[j].crd.y=0.0;
                for(c=0; c<3; c++) atoms[j].crd.y+=R[1][c]*x_old[c];
                atoms[j].crd.z=0.0;
                for(c=0; c<3; c++) atoms[j].crd.z+=R[2][c]*x_old[c];
		x_old[0]=atoms[j].vel.x;
                x_old[1]=atoms[j].vel.y;
                x_old[2]=atoms[j].vel.z;
                atoms[j].vel.x=0.0;
                for(c=0; c<3; c++) atoms[j].vel.x+=R[0][c]*x_old[c];
                atoms[j].vel.y=0.0;
                for(c=0; c<3; c++) atoms[j].vel.y+=R[1][c]*x_old[c];
                atoms[j].vel.z=0.0;
                for(c=0; c<3; c++) atoms[j].vel.z+=R[2][c]*x_old[c];
                /*vecRAdd(atoms[j].crd,refCOM,&tmp);
 *                 vecRCpy(tmp,&atoms[j].crd);*/
        }
        /* rotate box */
        x_old[0]=box[0].a.x;
        x_old[1]=box[0].a.y;
        x_old[2]=box[0].a.z;
        box[0].a.x=0.0;
        for(c=0; c<3; c++) box[0].a.x+=R[0][c]*x_old[c];
        box[0].a.y=0.0;
        for(c=0; c<3; c++) box[0].a.y+=R[1][c]*x_old[c];
        box[0].a.z=0.0;
        for(c=0; c<3; c++) box[0].a.z+=R[2][c]*x_old[c];
        x_old[0]=box[0].b.x;
        x_old[1]=box[0].b.y;
        x_old[2]=box[0].b.z;
        box[0].b.x=0.0;
        for(c=0; c<3; c++) box[0].b.x+=R[0][c]*x_old[c];
        box[0].b.y=0.0;
        for(c=0; c<3; c++) box[0].b.y+=R[1][c]*x_old[c];
        box[0].b.z=0.0;
        for(c=0; c<3; c++) box[0].b.z+=R[2][c]*x_old[c];
        x_old[0]=box[0].c.x;
        x_old[1]=box[0].c.y;
        x_old[2]=box[0].c.z;
        box[0].c.x=0.0;
        for(c=0; c<3; c++) box[0].c.x+=R[0][c]*x_old[c];
        box[0].c.y=0.0;
        for(c=0; c<3; c++) box[0].c.y+=R[1][c]*x_old[c];
        box[0].c.z=0.0;
        for(c=0; c<3; c++) box[0].c.z+=R[2][c]*x_old[c];

        makePBCtransMat(box);

        return 0;
}

int alignGroupNoRot(t_grp *grp,t_atom *atoms,int nAtoms,t_box box,t_vecR refCOM)
{
        int    i,j;
        t_vecR tmp;

	getGrpCOM(atoms,grp);
        for(i=0;i<nAtoms;i++) {
                vecRSub(grp[0].COM,atoms[i].crd,&tmp);
                vecRCpy(tmp,&atoms[i].crd);
        }

	/*for(j=0;j<nAtoms;j++) {
		vecRAdd(atoms[j].crd,refCOM,&tmp);
                vecRCpy(tmp,&atoms[j].crd);
        }*/
	return 0;
}

int alignGroupNoRotVel(t_grp *grp,t_atom *atoms,int nAtoms,t_box box,t_vecR refCOM)
{
        int    i,j;
        t_vecR tmp;

        getGrpCOM(atoms,grp);
	getGrpCOMvel(atoms,grp);
        for(i=0;i<nAtoms;i++) {
                vecRSub(grp[0].COM,atoms[i].crd,&tmp);
                vecRCpy(tmp,&atoms[i].crd);
		vecRSub(grp[0].COMvel,atoms[i].vel,&tmp);
		vecRCpy(tmp,&atoms[i].vel);
        }

        /*for(j=0;j<nAtoms;j++) {
                vecRAdd(atoms[j].crd,refCOM,&tmp);
                vecRCpy(tmp,&atoms[j].crd);
        }*/
        return 0;
}

