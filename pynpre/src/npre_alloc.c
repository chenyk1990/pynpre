#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

void *np_alloc (size_t n    /* number of elements */, 
			  size_t size /* size of one element */)
	  /*< output-checking allocation >*/
{
    void *ptr; 
    
    size *= n;
    
    ptr = malloc (size);

    if (NULL == ptr)
	{
	printf("cannot allocate %lu bytes:", size);
	return NULL;
	}

    return ptr;
}

float *np_floatalloc (size_t n /* number of elements */)
	  /*< float allocation >*/ 
{
    float *ptr;
    ptr = (float*) np_alloc (n,sizeof(float));
    return ptr;
}

int *np_intalloc (size_t n /* number of elements */)
	  /*< int allocation >*/  
{
    int *ptr;
    ptr = (int*) np_alloc (n,sizeof(int));
    return ptr;
}

typedef float complex np_complex;

np_complex *np_complexalloc (size_t n /* number of elements */) 
	  /*< complex allocation >*/
{
    np_complex *ptr;
    ptr = (np_complex*) np_alloc (n,sizeof(np_complex));
    return ptr;
}


