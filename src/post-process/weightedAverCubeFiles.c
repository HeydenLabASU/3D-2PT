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
	FILE *in,*in2,*out;
	char **incubes,**wcubes,outcube[100];
	real *grid,*grid2;
	char **header;
	char **atoms;
	int i,j,x,y,z,pos,natoms,nx,ny,nz,nread;
	char buffer[100];
	float tmp,tmp2;
	real start,end,current;
	int linespercol,count=0;

	printf("number of cube files:       :");
	fscanf(stdin,"%d",&nread);
	incubes=(char**)save_malloc(nread*sizeof(char*));
	wcubes=(char**)save_malloc(nread*sizeof(char*));
	for(i=0;i<nread;i++) incubes[i]=(char*)save_malloc(100*sizeof(char));
	for(i=0;i<nread;i++) wcubes[i]=(char*)save_malloc(100*sizeof(char));
	for(i=0;i<nread;i++)
	{
		printf("input *.cube file  %d    :",i+1);
		fscanf(stdin,"%s",incubes[i]);
		printf(" %s\n",incubes[i]);
		printf("weight *.cube file %d    :",i+1);
		fscanf(stdin,"%s",wcubes[i]);
		printf(" %s\n",wcubes[i]);
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
	grid2=(real*)save_malloc(nx*ny*nz*sizeof(real));
	for(i=0;i<nx*ny*nz;i++) {
		grid[i]=0.0;
		grid2[i]=0.0;
	}

	rewind(in);

	for(i=0;i<nread;i++)
	{
		if(i!=0) saveOpenRead(&in,incubes[i]);
		saveOpenRead(&in2,wcubes[i]);
		count++;
		for(j=0;j<natoms+6;j++)
		{
			fgets(buffer,100,in);
			fgets(buffer,100,in2);
		}
		for(x=0;x<nx;x++)
		{
			for(y=0;y<ny;y++)
			{
				for(z=0;z<nz;z++)
				{
					pos=z+nz*(y+ny*x);
					fscanf(in,"%f",&tmp);
					fscanf(in2,"%f",&tmp2);
					grid[pos]+=(tmp*tmp2);
					grid2[pos]+=tmp2;
				}
			}
		}
		fgets(buffer,100,in);
		fgets(buffer,100,in2);
		fclose(in);
		fclose(in2);
	}
	
	for(i=0;i<nx*ny*nz;i++)
	{
		if(grid2[i]!=0.0) {
			grid[i]=grid[i]/grid2[i];
		}
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

