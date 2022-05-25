// Read data file output from water3D and generate entropy cube files
// Assumes translational VDOS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define	PI	(3.14159265)
#define MAX_FREQS (200)

double func(double Delta, double x)
{
	double sqrtDelta, Delta2, sqrtx, x2;

	sqrtDelta = sqrt(Delta);
	Delta2 = Delta*Delta;
	
	sqrtx = sqrt(x);
	x2 = x*x;

	return 2./(Delta2*Delta2*sqrtDelta)*x2*x2*x2*x*sqrtx - 6. / (Delta2*Delta)*x2*x2*x - 1. / (Delta*sqrtDelta)*x2*x*sqrtx + 6. / (Delta*sqrtDelta)*x2*sqrtx + 2.*x - 2.;
}

int main(int argc, char *argv[])
{
	int i, n, nx, ny, nz, ix, iy, iz, inat, ngrid, imax, icell;
	char line[80];
	FILE *numdens, *vdosfile;
	double x0, y0, z0, x1, y1, z1, x2, y2, z2, dgrid, **vdos, *entr, *norm,
			*dens, k, lightspeed, kT, mass, planck, debroglie, x, y, z,
			df, g0, Delta, vdos_solid, vdos_gas, f, zfy, entr_cs, entr_st, 
			w_ho, w_gas, tmp, ox, oy, oz, dgrid3, one_dens = 0.1;

	k = 1.3806488e-23;
	lightspeed = 2.99792458e10; //cm/s
	kT=k*300;
	planck = 6.62606957e-34;

	if (argc < 3)
	{
		fprintf(stderr, "SYNTAX: %s {VDOS data file} {NumDens cube file} [density]\n",
				argv[0]);
		return 1;
	}
	if (argc == 4)
	{
		one_dens = atof(argv[3]);
	}
	else
	{
		one_dens = -1.;
	}

	// Open NumDens file and read the number of grid cells
	numdens = fopen(argv[2], "r");
	for (i = 0; i < 2; ++i)
		fgets(line, 80, numdens);

	printf("Translation entropy of water molecules\n");
	printf("Created by trans-vdos2entropy.x from data by water3D\n");

	fgets(line, 80, numdens);
	sscanf(line, "%i %lg %lg %lg", &inat, &ox, &oy, &oz);
	printf("%s", line);

	fgets(line, 80, numdens);
	sscanf(line, "%i %lg %lg %lg", &nx, &x0, &y0, &z0);
	printf("%s", line);

	fgets(line, 80, numdens);
	sscanf(line, "%i %lg %lg %lg", &ny, &x1, &y1, &z1);
	printf("%s", line);

	fgets(line, 80, numdens);
	sscanf(line, "%i %lg %lg %lg", &nz, &x2, &y2, &z2);
	printf("%s", line);

	dgrid = x0 * 0.529177249e-10;
	dgrid3= dgrid*dgrid*dgrid;

	for (i = 0; i < inat; ++i)
	{
		fgets(line, 80, numdens);
		printf("%s", line);
	}

	ngrid = nx*ny*nz;
	fprintf(stderr, "nx = %i ny = %i nz = %i\n", nx, ny, nz);
	fprintf(stderr, "Allocating VDOS (%i)\n", ngrid*MAX_FREQS);
	vdos = (double **)malloc(ngrid * sizeof(double *));
	for (i = 0; i < ngrid; ++i)
		vdos[i] = (double *)malloc(MAX_FREQS * sizeof(double));
	
	fprintf(stderr, "Allocating norm (%i)\n", ngrid);
	norm = (double *)calloc(ngrid, sizeof(double));
	fprintf(stderr, "Allocating dens (%i)\n", ngrid);
	dens = (double *)malloc(ngrid * sizeof(double));
	fprintf(stderr, "Allocating entr (%i)\n", ngrid);
	entr = (double *)calloc(ngrid, sizeof(double));

	if (one_dens < 0.)
	{
		for (ix = 0; ix < nx; ++ix)
		{
			for (iy = 0; iy < ny; ++iy)
			{
				for (iz = 0; iz < nz; ++iz)
				{
					fscanf(numdens, "%lg ", &tmp);
					tmp = tmp * 1e9 * 1e9 * 1e9;
					dens[iz + iy*nz + ix*nz*ny] = tmp;
				}
			}
		}
	}
	else
	{
		for (ix = 0; ix < nx; ++ix)
		{
			for (iy = 0; iy < ny; ++iy)
			{
				for (iz = 0; iz < nz; ++iz)
				{
					tmp = one_dens * 1e9 * 1e9 * 1e9;
					dens[iz + iy*nz + ix*nz*ny] = tmp;
				}
			}
		}
	}
	fclose(numdens);

	vdosfile = fopen(argv[1], "r");

	icell = -1;
	while (!feof(vdosfile))
	{
		++icell;
		fgets(line, 80, vdosfile);

		if (line[0] != '#')
		{
			fprintf(stderr, "\nFormat error: expected '#' but found '%c'\n",
							line[0]);
			return 1;
		}

		if (sscanf(line, "# %lg %lg %lg df= %lg n= %i m= %lg", &x, &y, &z, &df, &n, &mass) != 6)
		{
			if (sscanf(line, "# %lg %lg %lg df= %lg n= %i", &x, &y, &z, &df, &n) != 5)
			{
				fprintf(stderr, "\nFormat error\n");
				return 1;
			}
			else
			{
				if (icell == 0)
				{
					fprintf(stderr, 
						"\nWARNING: No mass found, resorting to water's\n");
					mass=18.01528;
				}
			}
		}

		df *= lightspeed; //to Hz from cm^-1

		i = -1;
		while (1)
		{
			++i;
			if (i > MAX_FREQS)
			{
				fprintf(stderr, "Insufficient memory allocation: freqs\n");
				return 1;
			}
			if (fscanf(vdosfile, "%lg", &tmp) != 1)
			{
				fscanf(vdosfile, "\n");
				break;
			}
			vdos[icell][i] = tmp;
			norm[icell] += tmp*df;
		}
		if (icell % 1000 == 0)
		{
			fprintf(stderr, "Reading progress: %.2f%% %c", 100.*(double)icell / (double)ngrid, 13);
		}
	}
	fclose(vdosfile);
	fprintf(stderr, "\n");

	mass *= 1.66053892e-27;
	debroglie = planck / sqrt(2.*PI*mass*kT);

	imax = i;
	for (icell = 0; icell < ngrid; ++icell)
	{
		if (norm[icell] == 0)
			continue;

		for (i = 0; i < imax; ++i)
		{
			vdos[icell][i] = vdos[icell][i] / norm[icell] * 3.;
		}
		g0 = vdos[icell][0];
		if (g0 <= 0)
			continue;

		Delta = 2.*g0 / (9.) * sqrt(PI*kT/mass) * pow(dens[icell], 1./3.)*pow(6./PI, 2./3.);

		// find f by numerical search here
		f = 0.;
		while (func(Delta, f) < 0)
			f += 0.001;
		f -= 0.001;
		while (func(Delta, f) < 0)
			f += 0.0001;
		f -= 0.0001;
		while (func(Delta, f) < 0)
			f += 0.00001;
		f -= 0.00001;
		while (func(Delta, f) < 0)
			f += 0.000001;

		if (f > 1.)
		{
			fprintf(stderr, "f out of bounds\n");
			return 1;
		}

		y = pow(f / Delta, 1.5);
		zfy = (1. + f*y + (f*y)*(f*y) - (f*y)*(f*y)*(f*y)) / ( (1. - f*y)*(1. - f*y)*(1. - f*y) );

		entr_cs = (log(zfy) + f*y*(3.*f*y - 4.) / ( (1. - f*y)*(1. - f*y) ) )*k;
		entr_st = (-log(f*dens[icell]*(debroglie*debroglie*debroglie)) + 2.5)*k;
		w_gas = (entr_cs + entr_st) / (3.*k);
		for (i = 0; i < imax; ++i)
		{
			w_ho = planck*(i+1)*df / (kT*(exp(planck*(i+1)*df/kT) - 1.)) - log(1. - exp(-planck*(i+1)*df / kT));

			vdos_gas = g0 / (1. + pow((g0*(i+1)*df*PI / (6.*f)), 2));
			vdos_solid = vdos[icell][i] - vdos_gas;

			entr[icell] = entr[icell] + w_ho*vdos_solid*df; 
			entr[icell] = entr[icell] + w_gas*vdos_gas*df; 
		}

		entr[icell] = entr[icell]*k;
		if (icell % 1000 == 0)
		{
			fprintf(stderr, "Calculation progress: %.2f%% %c",
				(double)icell / (double)ngrid * 100., 13);
		}
	}
	fprintf(stderr, "\n");

	for (ix = 0; ix < nx; ++ix)
	{
		for (iy = 0; iy < ny; ++iy)
		{
			for (iz = 0; iz < nz; ++iz)
			{
				// The units are J / (K mol)
				printf("%g ", entr[iz + iy*nz + ix*nz*ny]*6.0224e23);
				if (iz % 6 == 5)
					printf("\n");
			}
			if (iz % 6 != 0 ) printf("\n");
		}
	}

	return 0;
}
