#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef float real;

int fatal(const char *info)
{
        FILE *out;
        out=fopen("MH_FATAL_ERROR.LOG","w");
        fprintf(out,"%s\n",info);
        fclose(out);
        printf("%s\n",info);
        exit(1);
}

void *save_malloc(size_t size)
{
        void *ptr;
        char buffer[200];

        if((ptr=malloc(size))==NULL)
        {
                sprintf(buffer,"MEMORY ALLOCATION FAILED!!!\ntried to allocate %lu bytes\n -> save_malloc",size);
                fatal(buffer);
        }

        return ptr;
}

int saveOpenRead(FILE **io,char *fn)
{
        char buffer[200];
        if((io[0]=fopen(fn,"r"))==NULL)
        {
                sprintf(buffer,"FILE NOT FOUND!!!: %s\n -> saveOpenRead",fn);
                fatal(buffer);
        }
        return 0;
}

int main()
{
	FILE *in,*out;
	char **incubes,outcube[100];
	real *grid;
	char **header;
	char **atoms;
	int i,j,x,y,z,pos,natoms,nx,ny,nz,nread;
	char buffer[100];
	float tmp;
	real start,end,current;
	int linespercol,count=0;

	printf("number of cube files:       :");
	fscanf(stdin,"%d",&nread);
	incubes=(char**)save_malloc(nread*sizeof(char*));
	for(i=0;i<nread;i++) incubes[i]=(char*)save_malloc(100*sizeof(char));
	for(i=0;i<nread;i++)
	{
		printf("input *.cube file  %d    :",i+1);
		fscanf(stdin,"%s",incubes[i]);
		printf(" %s\n",incubes[i]);
	}
	printf("output *.cube file          :");
	fscanf(stdin,"%s",outcube);
	printf(" %s\n",outcube);

	header=(char**)save_malloc(6*sizeof(char*));
	for(i=0;i<6;i++) header[i]=(char*)save_malloc(100*sizeof(char));

	saveOpenRead(&in,incubes[0]);
	fgets(header[0],100,in);
	fgets(header[1],100,in);
	fgets(header[2],100,in);
	sscanf(header[2],"%d",&natoms);
	
	fgets(header[3],100,in);
	sscanf(header[3],"%d",&nx);
	fgets(header[4],100,in);
	sscanf(header[4],"%d",&ny);
	fgets(header[5],100,in);
	sscanf(header[5],"%d",&nz);

	atoms=(char**)save_malloc(natoms*sizeof(char*));
	for(i=0;i<natoms;i++) atoms[i]=(char*)save_malloc(100*sizeof(char));
	for(i=0;i<natoms;i++)
	{
		fgets(atoms[i],100,in);
	}

	grid=(real*)save_malloc(nx*ny*nz*sizeof(real));
	for(i=0;i<nx*ny*nz;i++) grid[i]=0.0;

	rewind(in);

	for(i=0;i<nread;i++)
	{
		if(i!=0) saveOpenRead(&in,incubes[i]);
		count++;
		for(j=0;j<natoms+6;j++)
		{
			fgets(buffer,100,in);
		}
		for(x=0;x<nx;x++)
		{
			for(y=0;y<ny;y++)
			{
				for(z=0;z<nz;z++)
				{
					pos=z+nz*(y+ny*x);
					fscanf(in,"%f",&tmp);
					grid[pos]+=tmp;
				}
			}
		}
		fgets(buffer,100,in);
		fclose(in);
	}
	
	for(i=0;i<nx*ny*nz;i++)
	{
		grid[i]=grid[i]/count;
	}

	out=fopen(outcube,"w");
	for(i=0;i<6;i++)
	{
		fprintf(out,"%s",header[i]);
	}
	for(i=0;i<natoms;i++)
	{
		fprintf(out,"%s",atoms[i]);
	}
	
	for(x=0;x<nx;x++)
	{
		for(y=0;y<ny;y++)
		{
			for(z=0;z<nz;z++)
			{
				pos=z+nz*(y+ny*x);
				fprintf(out,"%.6e ",grid[pos]);
				if(z % 6 == 5) fprintf(out,"\n");
			}
			fprintf(out,"\n");
		}
	}
	fclose(out);

	return 0;
}

