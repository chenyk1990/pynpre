/* 1D Fourier transform */
/*
  Copyright (C) 2016 Yangkang Chen
*/

#include<math.h>
#include <stdio.h>
#include <complex.h>
#include <stdbool.h>

#ifdef USE_SIMD
# include <xmmintrin.h>
# define kiss_fft_scalar __m128
#define KISS_FFT_MALLOC(nbytes) _mm_malloc(nbytes,16)
#define KISS_FFT_FREE _mm_free
#else	
#define KISS_FFT_MALLOC malloc
#define KISS_FFT_FREE free
#endif	


#ifdef FIXED_POINT
#include <sys/types.h>	
# if (FIXED_POINT == 32)
#  define kiss_fft_scalar int32_t
# else	
#  define kiss_fft_scalar int16_t
# endif
#else
# ifndef kiss_fft_scalar
/*  default is float */
#   define kiss_fft_scalar float
# endif
#endif

typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
}kiss_fft_cpx;


typedef struct kiss_fftr_state *kiss_fftr_cfg;

kiss_fftr_cfg kiss_fftr_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
void kiss_fftr(kiss_fftr_cfg cfg,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata);
void kiss_fftri(kiss_fftr_cfg cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);

#define SF_PI (3.14159265358979323846264338328)

int omp_init(void)
/*< init OMP parameters >*/
{
    int ompnth=1;
    
    return ompnth;
}

typedef float complex np_complex;
np_complex *np_complexalloc (size_t n /* number of elements */);
void *np_alloc (size_t n    /* number of elements */, size_t size /* size of one element */);
			  
kiss_fft_cpx np_cmul(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex multiplication >*/
{
    kiss_fft_cpx c;
    c.r = a.r*b.r-a.i*b.i;
    c.i = a.i*b.r+a.r*b.i;
    return c;
}

void fft1(float **ompP /*input data*/, 
		kiss_fft_cpx **ompQ,
		int n2,
		int n1,
		float d1,
		float o1, 		
		int nt,
		int nw,
		float dw,
		bool sym, 
		bool opt,
		bool verb, 		
		bool inv)	
/*<1D forward and inverse Fourier transform>*/    
{

    float shift;
    float wght;
    int   ompnth=1; /* number of threads */
    int   ompith=0; /* thread index */
    int i1,i2;
    kiss_fft_cpx  *ompE;
            
    ompnth=omp_init();
    
    kiss_fftr_cfg *ompcfg;

    /*------------------------------------------------------------*/
    /* Hermitian weight */
    wght = sym ? 1.0f/sqrtf((float) nt): 1.0f/nt;

    /*------------------------------------------------------------*/
    if(verb) printf("allocate arrays %d %d\n",n2,ompnth);

    ompE = (kiss_fft_cpx* ) np_complexalloc (ompnth);

    /*------------------------------------------------------------*/
    if(verb) printf("init FFT\n");


    ompcfg  = (kiss_fftr_cfg*) np_alloc(ompnth,sizeof(kiss_fftr_cfg));
    for(ompith=0; ompith<ompnth; ompith++)
	ompcfg[ompith] = kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
    ompith=0;


	if (!inv) { /* FORWARD TRANSFORM */

	    for(i2=0; i2<n2; i2++) {

		if (sym) for(i1=0;  i1<n1; i1++) ompP[i2][i1] *= wght;
		;        for(i1=n1; i1<nt; i1++) ompP[i2][i1]  = 0.0;

		kiss_fftr(ompcfg[ompith],ompP[i2],ompQ[i2]);


		if (0. != o1) { shift = -2.0*SF_PI*dw*o1;
		    for (i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[i2][i1]=np_cmul(ompQ[i2][i1],ompE[ompith]);
		    }
		}

	    }

	} else { /* INVERSE TRANSFORM */

	    for(i2=0; i2<n2; i2++) {

		if (0. != o1) { shift = +2.0*SF_PI*dw*o1;
		    for(i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[i2][i1]=np_cmul(ompQ[i2][i1],ompE[ompith]);
		    }
		}
		kiss_fftri(ompcfg[ompith],ompQ[i2],ompP[i2]);

		for(i1=0; i1<n1; i1++) ompP[i2][i1] *= wght;
	    }

	}
}

