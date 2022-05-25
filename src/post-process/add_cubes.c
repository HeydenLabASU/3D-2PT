// add cube files: nothing more, nothing less. Result to stdout.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int i, j, ix, iy, iz, nx, ny, nz, nfiles, ngrid;
	FILE *filelist, *file;
	char line[80];
	double *data, tmp;

	if (argc < 2)
	{
		fprintf(stderr, "SYNTAX: %s {file with files listed}\n", argv[0]);
		return 1;
	}

	filelist	=	fopen(argv[1], "r");
	fgets(line, 80, filelist);
	nfiles		=	atoi(line);

	for (i = 0; i < nfiles; ++i)
	{
		fgets(line, 80, filelist);
		line[strlen(line)-1]	=	0;
		//fprintf(stderr, "%s\n", line);
		file	=	fopen(line, "r");
		if (file == NULL)
		{
			fprintf(stderr, "Error opening '%s'\n", line);
			return 2;
		}
		if (i == 0)
		{
			// read and echo metadata (assumed to be the same for all
			// subsequent cube files)
			// comment lines
			fgets(line, 80, file);
			printf("averaged: %s", line);
			fgets(line, 80, file);
			printf("averaged: %s", line);

			// no. of atoms and position of origin
			fgets(line, 80, file);
			printf("%s", line);

			fgets(line, 80, file);
			sscanf(line, "%i", &nx);
			printf("%s", line);

			fgets(line, 80, file);
			sscanf(line, "%i", &ny);
			printf("%s", line);

			fgets(line, 80, file);
			sscanf(line, "%i", &nz);
			printf("%s", line);

			// we assume there is only one atom
			fgets(line, 80, file);
			printf("%s", line);

			ngrid 	= 	nx*ny*nz;
			data	=	(double *)calloc(ngrid, sizeof(double));
		}
		else
		{
			for (j = 0; j < 7; ++j)
				fgets(line, 80, file);
		}
		for (j = 0; j < ngrid; ++j)
		{
			if (fscanf(file, "%lf", &tmp) != 1)
			{
				fprintf(stderr, "What!? %i %i\n", j, ngrid);	
				return 3;
			}
			data[j]	+=	tmp;
		}

		fclose(file);
	}

	j	=	0;
	for (ix = 0; ix < nx; ++ix)
	{
		for (iy = 0; iy < ny; ++iy)
		{
			for (iz = 0; iz < nz; ++iz)
			{
				printf("%g ", data[j]);
				++j;
				if (iz % 6 == 5)
					printf("\n");
			}
			printf("\n");
		}
	}

	fclose(filelist);

	return 0;
}
