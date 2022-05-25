#include <string.h>
#include <math.h>
#include "../../include/GMXtrrio.h"
#include "../../include/fatal.h"
#include "../../include/dataTypes.h"
#include "../../include/geo.h"

int saveOpenRead(FILE **io,char *fn)
{
        char buffer[200];
        if((io[0]=fopen(fn,"r"))==NULL)
        {
                sprintf(buffer,"FILE NOT FOUND!!!: %s\n -> saveOpenRead\n -> io.c\n",fn);
                fatal(buffer);
        }
        return 0;
}

int saveOpenReadBin(FILE **io,char *fn)
{
        char buffer[200];
        if((io[0]=fopen(fn,"rb"))==NULL)
        {
                sprintf(buffer,"FILE NOT FOUND!!!: %s\n -> saveOpenReadBin\n -> io.c\n",fn);
                fatal(buffer);
        }
        return 0;
}

int getFormat(char *coord,int *format)
{
	char buffer[200];
	int namelen;
	char end[4];
	
	namelen=strlen(coord);
	strcpy(end,&coord[namelen-3]);

	if(strcmp(end,"xyz")==0 || strcmp(end,"XYZ")==0) format[0]=1;
	else if (strcmp(end,"crd")==0 || strcmp(end,"CRD")==0) format[0]=2;
	else if (strcmp(end,"gro")==0 || strcmp(end,"GRO")==0) format[0]=3;
	else if (strcmp(end,"trr")==0 || strcmp(end,"TRR")==0) format[0]=4;
	else if (strcmp(end,"dcd")==0 || strcmp(end,"DCD")==0) format[0]=5;
	else {
		sprintf(buffer,"FORMAT COULD NOT BE DETERMINED BY FILE EXTENSION OF: %s\n -> getFormat\n -> io.c\n",coord);
                fatal(buffer);
	}
	return 0;
}

int openInput(int format,int needVelo,FILE **cin,FILE **vin,XDR *xdrin,char *coord,char *velo)
{
	char buffer[200];

	/* 1: xyz (velocities in separate file, if needed)
	   2: AMBER crd with box (velocities in separate file, if needed)
	   3: GRO
	   4: TRR
	   5: DCD (velocities in separate file, if needed)
	*/
        if(format==1)
        {
                saveOpenRead(cin,coord);
                if(needVelo==1)
                {
                        saveOpenRead(vin,velo);
                }
	} else if(format==2)
	{
		saveOpenRead(cin,coord);
                if(needVelo==1)
                {
                        saveOpenRead(vin,velo);
                }
        } else if(format==3)
        {
                saveOpenRead(cin,coord);
        } else if(format==4)
        {
                if(xdropen(xdrin,coord,"rb")==0)
                {
			sprintf(buffer,"FILE NOT FOUND!: %s\n -> openInput\n -> io.c\n",coord);
                        fatal(buffer);
                }
	} else if(format==5)
	{
		saveOpenReadBin(cin,coord);
		if(needVelo==1)
		{
			saveOpenReadBin(vin,velo);
		}
        } else {
		sprintf(buffer,"UNKNOWN INPUT FORMAT DETERMINED BY FILE EXTENSION OF: %s\n -> openInput\n -> io.c\n",coord);
                fatal(buffer);
        }
        return 0;
}

int closeInput(int format,int needVelo,FILE *cin,FILE *vin,XDR *xdrin)
{
        if(format==1 || format==2 || format==5)
        {
                fclose(cin);
                if(needVelo==1)
                {
                        fclose(vin);
                }
        } else if(format==3)
        {
                fclose(cin);
        } else if(format==4)
        {
                xdrclose(xdrin);
        }
        return 0;
}

int checkGROInput(char *buf,t_atom atom,int frame,int element)
{
	int resNr,atomNr;
	char resName[6];
	char atomName[6];

	sscanf(buf,"%d",&resNr);
	if(resNr!=atom.resNr+1) {
		sprintf(buf,"residue number in GRO file and topology don't match!\nread: %d  expected: %d  frame: %d  element: %d\n -> readGRO\n -> io.c\n",resNr,atom.resNr+1,frame,element);
		fatal(buf);
	}
	sscanf(&buf[5],"%s",resName);
	if(strcmp(resName,atom.resName)!=0) {
		sprintf(buf,"residue name in GRO file and topology don't match!\nread: %s  expected: %s  frame: %d  element: %d\n -> readGRO\n -> io.c\n",resName,atom.resName,frame,element);
		fatal(buf);
	}
	sscanf(&buf[10],"%s",atomName);
	if(strcmp(atomName,atom.atomName)!=0) {
		sprintf(buf,"atom name in GRO file and topology don't match!\nread: %s  expected: %s  frame: %d  element: %d\n -> readGRO\n -> io.c\n",atomName,atom.atomName,frame,element);
		fatal(buf);
	}
	sscanf(&buf[15],"%d",&atomNr);
	if(atomNr!=atom.atomNr+1) {
		sprintf(buf,"atom number in GRO file and topology don't match!\nread: %d  expected: %d  frame: %d  element: %d\n -> readGRO\n -> io.c\n",atomNr,atom.atomNr+1,frame,element);
		fatal(buf);
	}
	return 0;
}

int readGRO(FILE *in,t_atom *atoms,int natoms,int needVelo,real *time,int frame,t_box *box)
{
        int i,j,natomstest;
        char buf[200];
        char *key;
        float tmp[6];

        fgets(buf,200,in);
        key=strstr(buf,"t=");
        if(key!=NULL)
        {
                if(sscanf(&key[2],"%f",&tmp[0])!=1) {
			printf("WARNING:   no time found in frame: %d!\n",frame);
			time[0]=(real)frame;
		} else {
			time[0]=tmp[0];
		}
        } else {
                time[0]=(real)frame;
		printf("WARNING:   no time found in frame: %d!\n",frame);
        }

        fgets(buf,200,in);
	sscanf(buf,"%d",&natomstest);
	if(natomstest!=natoms) {
        	strcat(buf," was read as atom number!!!\nUNEXPECTED NUMBER OF PARTICLES IN GRO FILE!!!\n -> readGRO\n -> io.c\n");
        	fatal(buf);
        }

	if(needVelo==1) {
		for(i=0;i<natoms;i++)
	        {
			fgets(buf,200,in);
			checkGROInput(buf,atoms[i],frame,i);
			if(sscanf(&buf[20],"%f %f %f %f %f %f",&tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5])!=6)
			{
				strcat(buf,"WRONG FORMAT IN GRO+VELO FILE!!!\n -> readGRO\n -> io.c\n");
                        	fatal(buf);
			}
			atoms[i].crd.x=(real)10.0*tmp[0];
			atoms[i].crd.y=(real)10.0*tmp[1];
			atoms[i].crd.z=(real)10.0*tmp[2];
			atoms[i].vel.x=(real)10.0*tmp[3];
			atoms[i].vel.y=(real)10.0*tmp[4];
			atoms[i].vel.z=(real)10.0*tmp[5];
		} 
	} else {
		for(i=0;i<natoms;i++)
                {
                        fgets(buf,200,in);
			checkGROInput(buf,atoms[i],frame,i);
                        if(sscanf(&buf[20],"%f %f %f",&tmp[0],&tmp[1],&tmp[2])!=3)
			{
                                strcat(buf,"WRONG FORMAT IN GRO FILE!!!\n -> readGRO\n -> io.c\n");
                                fatal(buf);
                        }
                        atoms[i].crd.x=(real)10.0*tmp[0];
                        atoms[i].crd.y=(real)10.0*tmp[1];
                        atoms[i].crd.z=(real)10.0*tmp[2];
                } 
	}
		
        fgets(buf,200,in);
        if(sscanf(buf,"%f %f %f",&tmp[0],&tmp[1],&tmp[2])!=3)
	{
        	strcat(buf,"WRONG FORMAT FOR BOX in GRO FILE!!!\n -> readGRO\n -> io.c\n");
        	fatal(buf);
       	}
        box[0].a.x=(real)10.0*tmp[0]; box[0].a.y=0.0; box[0].a.z=0.0;
        box[0].b.y=0.0; box[0].b.y=(real)10.0*tmp[1]; box[0].b.z=0.0;
        box[0].c.z=0.0; box[0].c.y=0.0; box[0].c.z=(real)10.0*tmp[2];
	box[0].ortho=1;
	makePBCtransMat(box);
        return 0;
}

int readTRR(XDR *in,t_atom *atoms,int natoms,int needVelo,real *time,int frame,t_box *box)
{
        int i,j;
        char buf[100];
        unsigned int max=100;
        char *ptr;
        int natomstest;
        int dum;
        float dum2;
        float vec[3];
        int velos;
	char err[200];

        if(frame==0)
        {
                printf("*************************************\n");
                printf("*Be aware that the TRR i/o routines *\n");
                printf("*you are using are only suitable for*\n");
                printf("*single precision trajectories!!!   *\n");
                printf("* M. Heyden                         *\n");
                printf("*************************************\n\n");
        }

        ptr=&buf[0];

        if(frame==0)
        {
                xdr_int(in,&dum);
        }
        xdr_string(in,&ptr,max);
        for(i=0;i<10;i++)
        {
                xdr_int(in,&dum);
        }
        xdr_int(in,&natomstest);
        if(natomstest!=natoms)
        {
		sprintf(err,"%d was read as atom number!!!\nUNEXPECTED NUMBER OF PARTICLES IN TRR FILE!!!\n -> readTRR\n -> io.c\n",natomstest);
                fatal(err);
        }
        xdr_int(in,&dum);
        xdr_int(in,&dum);
        xdr_float(in,&dum2);
        time[0]=(real)dum2;
        xdr_float(in,&dum2);

        xdr_vector(in,(char *)&vec,3,(unsigned int)sizeof(float),(xdrproc_t)xdr_float);
        box[0].a.x=10.0*(real)vec[0]; box[0].a.y=10.0*(real)vec[1]; box[0].a.z=10.0*(real)vec[2];
        xdr_vector(in,(char *)&vec,3,(unsigned int)sizeof(float),(xdrproc_t)xdr_float);
	box[0].b.x=10.0*(real)vec[0]; box[0].b.y=10.0*(real)vec[1]; box[0].b.z=10.0*(real)vec[2];
        xdr_vector(in,(char *)&vec,3,(unsigned int)sizeof(float),(xdrproc_t)xdr_float);
	box[0].c.x=10.0*(real)vec[0]; box[0].c.y=10.0*(real)vec[1]; box[0].c.z=10.0*(real)vec[2];
	if(box[0].a.y==0.0 && box[0].a.z==0.0 && box[0].b.x==0.0 && box[0].b.z==0.0 && box[0].c.x==0.0 && box[0].c.y==0.0) box[0].ortho=1;
	else box[0].ortho=0;
	makePBCtransMat(box);

        for(i=0;i<natoms;i++)
        {
                xdr_vector(in,(char *)&vec,3,(unsigned int)sizeof(float),(xdrproc_t)xdr_float);
                atoms[i].crd.x=10.0*vec[0];
                atoms[i].crd.y=10.0*vec[1];
                atoms[i].crd.z=10.0*vec[2];
        }

        xdr_int(in,&velos);             /* expected to be first item of header of the next structure, if there are no velos:  int with value 1993 */

        if(velos!=1993 && needVelo==1)
        {
                /* in this case you just read the x-component of the first velo-vector as an int*/
                /* set the stream back 4 bytes and start reading velocities*/
                fseek((FILE*)in->x_private,-4,SEEK_CUR);

                for(i=0;i<natoms;i++)
                {
                        xdr_vector(in,(char *)&vec,3,(unsigned int)sizeof(float),(xdrproc_t)xdr_float);
                        atoms[i].vel.x=10.0*vec[0];
                        atoms[i].vel.y=10.0*vec[1];
                        atoms[i].vel.z=10.0*vec[2];
                }
                xdr_int(in,&velos);     /* expected to be first item of header of the next structure:  int with value 1993 */
	} else if(velos!=1993) {
		/* in this case you just read the x-component of the first velo-vector as an int*/
                /* for single precision this doesn't matter because it has 4bytes size as well as a float*/
                /* since you don't want to use velocities just read the next 2 components and then natoms-1 times the remaining velo-vectors */
                xdr_int(in,&velos);
                xdr_int(in,&velos);
                if(frame==0)
                {
                        printf("->  velocities found, but ignored!\n");
                }
                for(i=1;i<natomstest;i++)
                {
                        xdr_vector(in,(char *)&vec,3,(unsigned int)sizeof(float),(xdrproc_t)xdr_float);
                }
                xdr_int(in,&velos);  /* expected to be first item of header of the next structure:  int with value 1993 */
        } else if(needVelo==1) {
                sprintf(err,"NO VELOCITIES FOUND IN TRR FILE FOR TIME FRAME: %d!!!\n -> readTRR\n -> io.c\n",frame);
                fatal(err);
        }

        return 0;
}

int checkXYZInput(char *buf,t_atom atom,int frame,int element)
{
	char elem[3];

	sscanf(buf,"%s",elem);
	if(strncmp(elem,atom.atomName,1)!=0) {
		sprintf(buf,"elements in XYZ file and topology don't match!\nread: %s  expected: %s  frame: %d  element: %d\n -> checkXYZInput\n -> io.c\n",elem,atom.atomName,frame,element);
                fatal(buf);
	}
	return 0;
}

int readXYZcrd(FILE *in,t_atom *atoms,int natoms,real *time,int frame)
{
        int i,natomstest;
        char buffer[300];
        char name[5];
        float x,y,z;
	char *key;
	float tmp;

        fgets(buffer,300,in);
        sscanf(buffer,"%d",&natomstest);
        if(natoms!=natomstest)
        {
                strcat(buffer," was read as atom number!!!\nUNEXPECTED NUMBER OF PARTICLES IN CRD XYZ FILE!!!\n -> readXYZcrd\n -> io.c\n");
                fatal(buffer);
        }

	fgets(buffer,300,in);
        key=strstr(buffer,"time =");
        if(key!=NULL)
        {
                if(sscanf(&key[6],"%f",&tmp)!=1) {
                        printf("WARNING:   no time found in frame: %d!\n",frame);
                        time[0]=(real)frame;
                } else {
                        time[0]=tmp;
                }
        } else {
                time[0]=(real)frame;
                if(frame==0) printf("WARNING:   no time found in frame: %d!\n",frame);
        }

        for(i=0;i<natoms;i++)
        {
                fgets(buffer,300,in);
		/*no check when reading into buffer: MH2013*/
		/*checkXYZInput(buffer,atoms[i],frame,i);*/
                if(sscanf(buffer,"%s %f %f %f",name,&x,&y,&z)!=4)
                {
                        strcat(buffer,"WRONG FORMAT IN CRD XYZ FILE!!!\n -> readXYZcrd\n -> io.c\n");
                        fatal(buffer);
                }
                atoms[i].crd.x=(real)x; atoms[i].crd.y=(real)y; atoms[i].crd.z=(real)z;
        }
        return 0;
}

int readXYZcrdWC(FILE *in,t_atom *atoms,int natoms,t_vecR *WCcrd,int nWC,real *time,int frame)
{
        int i,natomstest;
        char buffer[300];
        char name[5];
        float x,y,z;
        char *key;
	float tmp;

        fgets(buffer,300,in);
        sscanf(buffer,"%d",&natomstest);
        if(natoms+nWC!=natomstest)
        {
                strcat(buffer," was read as atom+wc number!!!\nUNEXPECTED NUMBER OF PARTICLES IN CRD+WANNIER XYZ FILE!!!\n -> readXYZcrdWC\n -> io.c\n");
                fatal(buffer);
        }

        fgets(buffer,300,in);
        key=strstr(buffer,"time =");
        if(key!=NULL)
        {
                if(sscanf(&key[6],"%f",&tmp)!=1) {
                        printf("WARNING:   no time found in frame: %d!\n",frame);
                        time[0]=(real)frame;
                } else {
                        time[0]=tmp;
                }
        } else {
                time[0]=(real)frame;
                if(frame==0) printf("WARNING:   no time found in frame: %d!\n",frame);
        }

        for(i=0;i<natoms;i++)
        {
                fgets(buffer,300,in);
		/*no check when reading into buffer: MH2013*/
		/*checkXYZInput(buffer,atoms[i],frame,i);*/
                if(sscanf(buffer,"%s %f %f %f",name,&x,&y,&z)!=4)
                {
                        strcat(buffer,"WRONG FORMAT IN CRD XYZ FILE!!!\n -> readXYZcrdWC\n -> io.c\n");
                        fatal(buffer);
                }
                atoms[i].crd.x=(real)x; atoms[i].crd.y=(real)y; atoms[i].crd.z=(real)z;
        }
	for(i=0;i<nWC;i++)
        {
                fgets(buffer,300,in);
                if(sscanf(buffer,"%s %f %f %f",name,&x,&y,&z)!=4)
                {
                        strcat(buffer,"WRONG FORMAT IN CRD XYZ FILE!!!\n -> readXYZcrdWC\n -> io.c\n");
                        fatal(buffer);
                }
		if(strcmp(name,"X")!=0) {
			sprintf(buffer,"expected Wannier center, but encountered something else\nread: %s  expected: X  frame: %d  element: %d\n -> readXYZcrdWC\n -> io.c\n",name,frame,natoms+i);
			fatal(buffer);
		}
                WCcrd[i].x=(real)x; WCcrd[i].y=(real)y; WCcrd[i].z=(real)z;
        }
	
        return 0;
}

int readXYZvel(FILE *in,t_atom *atoms,int natoms,int frame)
{
        int i,natomstest;
        char buffer[300];
        char name[5];
        float x,y,z;

        fgets(buffer,300,in);
        sscanf(buffer,"%d",&natomstest);
        if(natoms!=natomstest)
        {
                strcat(buffer," was read as atom number!!!\nUNEXPECTED NUMBER OF PARTICLES IN VELO XYZ FILE!!!\n -> readXYZvel\n -> io.c\n");
                fatal(buffer);
        }
	fgets(buffer,300,in);

        for(i=0;i<natoms;i++)
        {
                fgets(buffer,300,in);
		/*no check when reading into buffer: MH2013*/
		/*checkXYZInput(buffer,atoms[i],frame,i);*/
                if(sscanf(buffer,"%s %f %f %f",name,&x,&y,&z)!=4)
                {
                        strcat(buffer,"WRONG FORMAT IN VELO XYZ FILE!!!\n -> readXYZvel\n -> io.c\n");
                        fatal(buffer);
                }
                atoms[i].vel.x=(real)(21876.9*x); atoms[i].vel.y=(real)(21876.9*y); atoms[i].vel.z=(real)(21876.9*z);
        }
        return 0;
}

int readXYZforce(FILE *in,t_atom *atoms,int natoms,int frame)
{
        int i,natomstest;
        char buffer[300];
        char name[5];
        float x,y,z;

        fgets(buffer,300,in);
        sscanf(buffer,"%d",&natomstest);
        if(natoms!=natomstest)
        {
                strcat(buffer," was read as atom number!!!\nUNEXPECTED NUMBER OF PARTICLES IN FORCE XYZ FILE!!!\n -> readXYZforce\n -> io.c\n");
                fatal(buffer);
        }
        fgets(buffer,300,in);

        for(i=0;i<natoms;i++)
        {
                fgets(buffer,300,in);
		/*no check when reading into buffer: MH2013*/
                /*checkXYZInput(buffer,atoms[i],frame,i);*/
                if(sscanf(buffer,"%s %f %f %f",name,&x,&y,&z)!=4)
                {
                        strcat(buffer,"\nWRONG FORMAT IN FORCE XYZ FILE!!!\n -> readXYZforce\n -> io.c\n");
                        fatal(buffer);
                }
                atoms[i].force.x=(real)x; atoms[i].force.y=(real)y; atoms[i].force.z=(real)z;
        }
        return 0;
}

int readXYZboxOrtho(FILE *in,t_box *box,int frame)
{
	float time,x,y,z;
	char buffer[300];

	fgets(buffer,300,in);
	if(sscanf(buffer,"%f %f %f %f",&time,&x,&y,&z)!=4)
        {
        	strcat(buffer,"\nWRONG FORMAT IN BOX XYZ FILE!!!\n -> readXYZbox\n -> io.c\n");
        	fatal(buffer);
        }
	box[0].a.x=(real)x; box[0].a.y=0.0; box[0].a.z=0.0;
	box[0].b.x=0.0; box[0].b.y=(real)y; box[0].b.z=0.0;
	box[0].c.x=0.0; box[0].c.y=0.0; box[0].c.z=(real)z;
	box[0].ortho=1;
	makePBCtransMat(box);
	return 0;
}

int readXYZboxFull(FILE *in,t_box *box,int frame)
{
        float time,ax,ay,az,bx,by,bz,cx,cy,cz;
        char buffer[300];

        fgets(buffer,300,in);
        if(sscanf(buffer,"%f %f %f %f %f %f %f %f %f %f",&time,&ax,&ay,&az,&bx,&by,&bz,&cx,&cy,&cz)!=10)
        {
                strcat(buffer,"\nWRONG FORMAT IN BOX XYZ FILE!!!\n -> readXYZboxFull\n -> io.c\n");
                fatal(buffer);
        }
        box[0].a.x=(real)ax; box[0].a.y=(real)ay; box[0].a.z=(real)az;
        box[0].b.x=(real)bx; box[0].b.y=(real)by; box[0].b.z=(real)bz;
        box[0].c.x=(real)cx; box[0].c.y=(real)cy; box[0].c.z=(real)cz;
	if(box[0].a.y==0.0 && box[0].a.z==0.0 && box[0].b.x==0.0 && box[0].b.z==0.0 && box[0].c.x==0.0 && box[0].c.y==0.0) box[0].ortho=1;
        else box[0].ortho=0;
        makePBCtransMat(box);
        return 0;
}

int readAMBERcrd(FILE *in,t_atom *atoms,int natoms,real *time,int frame,real dt,t_box *box)
{
        int i;
        char buf[200];
        float dum;

	if(frame==0) {
		fgets(buf,200,in);
	}

        time[0]=frame*dt;

        for(i=0;i<natoms;i++)
        {
                fscanf(in,"%f",&dum);
                atoms[i].crd.x=(real)dum;
                fscanf(in,"%f",&dum);
                atoms[i].crd.y=(real)dum;
                fscanf(in,"%f",&dum);
                atoms[i].crd.z=(real)dum;
        }
        fgets(buf,200,in);

        fscanf(in,"%f",&dum);
        box[0].a.x=(real)dum;
	box[0].a.y=0.0;
	box[0].a.z=0.0;
        fscanf(in,"%f",&dum);
	box[0].b.x=0.0;
        box[0].b.y=(real)dum;
	box[0].b.z=0.0;
        fscanf(in,"%f",&dum);
	box[0].c.x=0.0;
	box[0].c.y=0.0;
        box[0].c.z=(real)dum;
	fgets(buf,200,in);
	box[0].ortho=1;
        makePBCtransMat(box);

        return 0;
}

int readAMBERvel(FILE *in,t_atom *atoms,int natoms,int frame)
{
        int i;
        char buf[200];
        float dum;

	if(frame==0) {
                fgets(buf,200,in);
        }

        for(i=0;i<natoms;i++)
        {
                fscanf(in,"%f",&dum);
                atoms[i].vel.x=(real)20.455*dum;
                fscanf(in,"%f",&dum);
                atoms[i].vel.y=(real)20.455*dum;
                fscanf(in,"%f",&dum);
                atoms[i].vel.z=(real)20.455*dum;
        }
        fgets(buf,200,in);

        return 0;
}

int getInt(FILE *in,int *i)
{
	fread(i,sizeof(int),1,in);
	return 0;
}

int getString(FILE *in,char *buf,int bufLen)
{
	fread(buf,sizeof(char),bufLen,in);
	return 0;
}

int getReal(FILE *in,real *r)
{
	fread(r,sizeof(real),1,in);
	return 0;
}

int getDouble(FILE *in,double *d)
{
	fread(d,sizeof(double),1,in);
	return 0;
}

int readDCDheader(FILE *in,int *nFrames,real *time,int nAtoms,real *dt)
{
	char err[300];
        char buffer[300];
        int i,n,start,interval,titleLen;
	real r,s;

        getInt(in,&i);
        if(i!=84) 
        {
                sprintf(err,"EXPECTED INTEGER 84 FROM FIRST 4 BYTES OF DCD FILE: %d FOUND INSTEAD!\n -> readDCDheader\n -> io.c\n",i);
                fatal(err);
        }
        getString(in,buffer,4); 
        buffer[4]=(char)0;
        if(strcmp(buffer,"CORD")!=0) 
        {
                fatal("MISSING KEYWORD 'CORD' IN HEADER OF DCD FILE!\n -> readDCDheader\n -> io.c\n");
        }
	getInt(in,nFrames);
	getInt(in,&start);
	getInt(in,&interval);
	for(n=0;n<6;n++)
	{
		getInt(in,&i);
	}
	getReal(in,&r);
	dt[0]=r/20.455*interval; time[0]=start*r/20.455;

	for(n=0;n<10;n++) getInt(in,&i);
	getInt(in,&i);
	getInt(in,&titleLen);
	if(titleLen>299) fatal("ERROR: TITLE EXCEEDS BUFFER LEN (299)!\n -> readDCDheader\n -> io.c\n");
	getString(in,buffer,titleLen);
	buffer[titleLen]=(char)0;
	getInt(in,&i);
	getInt(in,&i);
	getInt(in,&i);
	if(i!=nAtoms)
	{
		sprintf(err,"NUMBER OF ATOMS BETWEEN DCD (%d) AND TOP (%d) FILE DISAGREE!\n -> readDCDheader\n -> io.c\n",i,nAtoms);
		fatal(err);
	}
	getInt(in,&i);
	return 0;
}

int readDCDcrd(FILE *in,t_atom *atoms,int nAtoms,real *time,int frame,real dt,t_box *box)
{
	char err[300];
	int i;
	double d1,d2;

	/*if(frame!=0) time[0]+=dt;*/
	time[0]=frame*dt;

	getInt(in,&i);
	if(i!=48)
	{
		sprintf(err,"BOX INFORMATION MISSING IN DCD FILE FOR FRAME %d!\n -> readDCDcrd\n -> io.c\n",frame);
		fatal(err);
	}
	getDouble(in,&d1);
	box[0].a.x=(real)d1;
	box[0].a.y=0.0;
	box[0].a.z=0.0;
	getDouble(in,&d1);
	getDouble(in,&d2);
	box[0].b.x=(real)d1*d2;
	d1=sin(acos(d1));
	box[0].b.y=(real)d1*d2;
	box[0].b.z=0.0;
	getDouble(in,&d1);
	if(d1!=0.0) fatal("UNSUPPORTED BOX GEOMETRY: BETA OR GAMMA !=0\n -> readDCDcrd\n -> io.c\n");
	getDouble(in,&d1);
	if(d1!=0.0) fatal("UNSUPPORTED BOX GEOMETRY: BETA OR GAMMA !=0\n -> readDCDcrd\n -> io.c\n");
	getDouble(in,&d1);
	box[0].c.x=0.0;
	box[0].c.y=0.0;
	box[0].c.z=(real)d1;
	if(box[0].a.y==0.0 && box[0].a.z==0.0 && box[0].b.x==0.0 && box[0].b.z==0.0 && box[0].c.x==0.0 && box[0].c.y==0.0) box[0].ortho=1;
        else box[0].ortho=0;
        makePBCtransMat(box);
	getInt(in,&i);
	getInt(in,&i);
	if(i/4!=nAtoms)
	{
		sprintf(err,"%d BYTES FOR X COORDINATES FOUND! EXPECTED %d! FRAME: %d\n -> readDCDcrd\n -> io.c\n",i,nAtoms*4,frame);
		fatal(err);
	}
	for(i=0;i<nAtoms;i++) getReal(in,&atoms[i].crd.x);
	getInt(in,&i);
	getInt(in,&i);
	if(i/4!=nAtoms)
        {
                sprintf(err,"%d BYTES FOR Y COORDINATES FOUND! EXPECTED %d! FRAME: %d\n -> readDCDcrd\n -> io.c\n",i,nAtoms*4,frame);
                fatal(err);
        }
	for(i=0;i<nAtoms;i++) getReal(in,&atoms[i].crd.y);
	getInt(in,&i);
        getInt(in,&i);
	if(i/4!=nAtoms)
        {
                sprintf(err,"%d BYTES FOR Z COORDINATES FOUND! EXPECTED %d! FRAME: %d\n -> readDCDcrd\n -> io.c\n",i,nAtoms*4,frame);
                fatal(err);
        }
	for(i=0;i<nAtoms;i++) getReal(in,&atoms[i].crd.z);
	getInt(in,&i);
	
	return 0;
}

int readDCDvel(FILE *in,t_atom *atoms,int nAtoms,int frame,real dt)
{
        char err[300];
        int i;
        double d1,d2;

        getInt(in,&i);
        if(i/4!=nAtoms)
        {
                sprintf(err,"%d BYTES FOR X VELOCITIES FOUND! EXPECTED %d! FRAME: %d\n -> readDCDcrd\n -> io.c\n",i,nAtoms*4,frame);
                fatal(err);
        }
        for(i=0;i<nAtoms;i++) { 
		getReal(in,&atoms[i].vel.x);
		atoms[i].vel.x*=20.455;
	}
        getInt(in,&i);
        getInt(in,&i);
        if(i/4!=nAtoms)
        {
                sprintf(err,"%d BYTES FOR Y VELOCITIES FOUND! EXPECTED %d! FRAME: %d\n -> readDCDcrd\n -> io.c\n",i,nAtoms*4,frame);
                fatal(err);
        }
        for(i=0;i<nAtoms;i++) {
		getReal(in,&atoms[i].vel.y);
		atoms[i].vel.y*=20.455;
	}
        getInt(in,&i);
        getInt(in,&i);
        if(i/4!=nAtoms)
        {
                sprintf(err,"%d BYTES FOR Z VELOCITIES FOUND! EXPECTED %d! FRAME: %d\n -> readDCDcrd\n -> io.c\n",i,nAtoms*4,frame);
                fatal(err);
        }
	for(i=0;i<nAtoms;i++) {
		getReal(in,&atoms[i].vel.z);
		atoms[i].vel.z*=20.455;
	}
        getInt(in,&i);

        return 0;
}

int prepTRJinput(char *fnCrd,char *fnVel,FILE **cin,FILE **vin,XDR *xdrin,int *format,int *veloFreq,real *time,real *dt,int *nFrames,int needVelo,int nAtoms,int argc)
{
        real dtVel;
        int nFramesVel;
        real timeVel;
        char err[300];
        int tmp;
        float tmp2;
	char test[100];

        getFormat(fnCrd,format);
        if(format[0]==1 || format[0] ==2 || format[0]==5) /*XYZ, AMBER (CRD+BOX) OR DCD*/
        {
                if(needVelo==1 && sscanf(fnVel,"%s",test)!=1) fatal("VELOCITIES NEEDED BUT NO INPUT FILE SUPPLIED!\n -> prepTRJinput\n -> io.c\n");
        }
        openInput(format[0],needVelo,cin,vin,xdrin,fnCrd,fnVel);

        if(format[0]==5)
        {
                if(nFrames[0]==0) readDCDheader(cin[0],nFrames,time,nAtoms,dt);
                else readDCDheader(cin[0],&tmp,time,nAtoms,dt);
        }
        if(format[0]==2)
        {
                printf("AMBER TRAJECTORY TIME STEP:     ");
                fscanf(stdin,"%f",&tmp2);
                dt[0]=tmp2;
        }
        if(format[0]==5 && needVelo==1)
        {
                readDCDheader(vin[0],&nFramesVel,&timeVel,nAtoms,&dtVel);
                if(time[0]!=timeVel) fatal("DIFFERENT TIME STAMPS IN CRD AND VEL FILE FOR FIRST FRAME!\n -> prepTRJinput\n -> io.c\n");
                if(dt[0]!=dtVel)
                {
                        printf("WARNING: DIFFERENT TIME STEPS IN CRD AND VEL FILE!\n");
                        printf("         CRD TIME STEP: %f   VEL TIME STEP: %f\n",dt[0],dtVel);
                        printf("         ENFORCING REDUCED FREQUENCY OF VELOCITY INPUT:\n");
                        veloFreq[0]=(int)dtVel/dt[0];
                        printf("         'veloFreq: dtVel/dt = %d\n",veloFreq[0]);
                }
                if(nFrames[0]>veloFreq[0]*nFramesVel)
                {
                        sprintf(err,"INSUFFICIENT NUMBER OF VELOCITY FRAMES:\nFRAMES IN CRD DCD FILE: %d\nFRAMES IN VEL DCD FILE: %d\nREADING VELOCITIES EVERY %d STEP\n -> prepTRJinput\n -> io.c\n",nFrames[0],nFramesVel,veloFreq[0]);
                        fatal(err);
                }
        }
        return 0;
}

int readTRJ(int frame,int format,FILE *cin,FILE *vin,XDR *xdrin,real dt,real *time,t_box *box,int nAtoms,int nSites,int needVelo,int veloFreq,int wc,t_atom *atoms,t_vecR *siteCrds)
{
	if(format==1)
	{
		if(wc==0) readXYZcrd(cin,atoms,nAtoms,time,frame);
		else readXYZcrdWC(cin,atoms,nAtoms,siteCrds,nSites,time,frame);

		if(needVelo==1 && frame%veloFreq==0) readXYZvel(vin,atoms,nAtoms,frame);
	}

	if(format==2)
	{
		readAMBERcrd(cin,atoms,nAtoms,time,frame,dt,box);

		if(needVelo==1 && frame%veloFreq==0) readAMBERvel(vin,atoms,nAtoms,frame);
	}

	if(format==3) readGRO(cin,atoms,nAtoms,needVelo,time,frame,box);

	if(format==4) readTRR(xdrin,atoms,nAtoms,needVelo,time,frame,box);

	if(format==5)
	{
		readDCDcrd(cin,atoms,nAtoms,time,frame,dt,box);

		if(needVelo==1 && frame%veloFreq==0) readDCDvel(vin,atoms,nAtoms,frame,dt);
	}

	return 0;
}

int printVolumeData(char *fn,char *title,t_vecR ori,t_box vox,int n1,int n2,int n3,real ***data) {
        FILE *cube;
        int a,b,c;
        char name[100];
	t_vecR oriOut;

	/*center of voxel[0][0][0], instead of corner (used internally)*/
	oriOut.x=ori.x+vox.a.x/2.0;
	oriOut.y=ori.y+vox.b.y/2.0;
	oriOut.z=ori.z+vox.c.z/2.0;

        cube=fopen(fn,"w");
        fprintf(cube,"%s\n",title);
        fprintf(cube,"created by M. Heyden\n");
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",1,oriOut.x/0.5292,oriOut.y/0.5292,oriOut.z/0.5292);
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",n1,vox.a.x/0.5292,vox.a.y/0.5292,vox.a.z/0.5292);
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",n2,vox.b.x/0.5292,vox.b.y/0.5292,vox.b.z/0.5292);
        fprintf(cube,"%5d%12.6f%12.6f%12.6f\n",n3,vox.c.x/0.5292,vox.c.y/0.5292,vox.c.z/0.5292);
        fprintf(cube,"%d%12.6f%12.6f%12.6f%12.6f\n",8,0.0,0.0,0.0,0.0);
        for (a=0;a<n1;a++) {
                for (b=0;b<n2;b++) {
                        for (c=0;c<n3;c++) {
                                fprintf(cube," %12.5e",(double)data[a][b][c]);
                                if (c % 6 == 5) fprintf(cube,"\n");
                        }
                        fprintf(cube,"\n");
                }
        }
        fclose(cube);
        sprintf(name,"%-50s",fn);
        printf("   wrote file: %s\n",name);
        return 0;
}
