// Romberg integration utilities

#include <stdlib.h>
#include <stdio.h>

#include "Utils.h"


//**********************************************************************************


#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
// Numerical Recipes standard error handler 
{
	fprintf(stderr,"Run-time error...\n");
	fprintf(stderr,"%s\n", error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


//************************************************************************************


double *vector(long nl, long nh)
// allocate a double vector with subscript range v[nl..nh] 
{
	double *v;

	v = (double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}


//************************************************************************************


void free_vector(double *v, long nl, long nh)
// free a double vector allocated with vector() 
{
        free((FREE_ARG) (v+nl-NR_END));
}


#undef NR_END
#undef FREE_ARG


//************************************************************************************
