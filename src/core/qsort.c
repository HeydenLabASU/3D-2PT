#include <stdlib.h>
#include <string.h>

/* qsort int comparison function */
int int_cmp(const void *a, const void *b)
{
        const int *ia = (const int *)a; /* casting pointer types  */
        const int *ib = (const int *)b;
        return *ia  - *ib;
        /* integer comparison: returns negative if b > a 
        and positive if a > b */
}

/* qsort real comparison function */
int real_cmp(const void *a, const void *b)
{
	int res;
        const float *ia = (const float *)a; /* casting pointer types  */
        const float *ib = (const float *)b;
        if(*ia-*ib<0.0)
	{
		res=-1;
	} else res=1;
	return res;
        /* real comparison: returns negative if b > a 
        and positive if a > b */
}

/* qsort C-string comparison function */
int cstring_cmp(const void *a, const void *b)
{
        const char **ia = (const char **)a;
        const char **ib = (const char **)b;
        return strcmp(*ia, *ib);
        /* strcmp functions works exactly as expected from
        comparison function */
}

/*USAGE OF QSORT*/
/* qsort(array,nelem,sizeof(arrayelem),oneOfTheFunctionsAbove); */
