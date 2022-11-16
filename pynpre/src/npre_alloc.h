

#ifndef KISS_FFT_H
#include "npre_kissfft.h"
#endif

#include "npre_dtype.h"

void *np_alloc (size_t n    /* number of elements */, 
			  size_t size /* size of one element */);

float *np_floatalloc (size_t n /* number of elements */);

float **np_floatalloc2 (size_t n1 /* fast dimension */, 
				  size_t n2 /* slow dimension */);

float ***np_floatalloc3 (size_t n1 /* fast dimension */, 
				   size_t n2 /* slower dimension */, 
				   size_t n3 /* slowest dimension */);

int *np_intalloc (size_t n /* number of elements */);

int **np_intalloc2 (size_t n1 /* fast dimension */, 
			      size_t n2 /* slow dimension */);

int ***np_intalloc3 (size_t n1 /* fast dimension */, 
			       size_t n2 /* slower dimension */, 
			       size_t n3 /* slowest dimension */);
			       
// np_complex *np_complexalloc (size_t n /* number of elements */);

np_complex *np_complexalloc (size_t n /* number of elements */);

np_complex **np_complexalloc2 (size_t n1 /* fast dimension */, 
					 size_t n2 /* slow dimension */);

np_complex ***np_complexalloc3 (size_t n1 /* fast dimension */, 
					  size_t n2 /* slower dimension */, 
					  size_t n3 /* slowest dimension */);