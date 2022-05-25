#include <stdio.h>
#include <stdlib.h>

int fatal(const char *info)
{
	FILE *out;
	out=fopen("MH_FATAL_ERROR.LOG","w");
	fprintf(out,"%s\n",info);
	fclose(out);
	printf("%s\n",info);
	exit(1);
}

