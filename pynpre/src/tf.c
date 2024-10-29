#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "ntfa_divnnsc.h"

#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))
#define SF_MAX_DIM 9
#define SF_PI (3.14159265358979323846264338328)

#include "ntfa_alloc.h"

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

static PyObject *stft1d(PyObject *self, PyObject *args){
    
	/**initialize data input**/
    int nd, nd2;
    
    PyObject *f1=NULL;
    PyObject *arrf1=NULL;

	int ndata;	/*integer parameter*/
	float fpar; /*float parameter*/
    int ndim, i;
    float *data;
    
    int niter,verb,rect0,n1,ntw,opt,sym,window;
    float dt,alpha;
    int ifb,inv;
	PyArg_ParseTuple(args, "Oiiiiiiff", &f1,&n1,&verb,&window,&inv,&sym,&opt,&ntw,&dt,&ot);
	
    if (ntw%2 == 0)
        ntw = (ntw+1);
    
	ndata=n1;


    int i1, iw, nt, nw, i2, n2, n12, n1w;
    int m[SF_MAX_DIM], *rect;
    float t, d1, w, w0, dw, mean=0.0f;
    float *inp, *kbsc, *mkbsc, *sscc, *mm, *ww;
    
    d1=dt;
    if(opt)
	    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
    else
        nt=ntw;
    if (nt%2) nt++;
    nw = nt/2+1;
    dw = 1./(nt*d1);
	ow = 0.;


	for(i2=0; i2 < SF_MAX_DIM; i2 ++) {
	    m[i2] = 1;
	}
	m[0] = n1;

    n1w = n1*nw;
    n12 = 2*n1w;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;


	printf("niter=%d,n1=%d,nw=%d,nt=%d,n12=%d,rect0=%d\n",niter,n1,nw,nt,n12,rect0);
	printf("dw=%g,w0=%g,alpha=%g,dt=%g\n",dw,w0,alpha,dt);


    inp = tf_floatalloc(n1);

	    
	if(!inv)
	{
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    
    nd2=PyArray_NDIM(arrf1);
    
    npy_intp *sp=PyArray_SHAPE(arrf1);

	data  = (float*)malloc(ndata * sizeof(float));
	
    if (*sp != ndata)
    {
    	printf("Dimension mismatch, N_input = %d, N_data = %d\n", *sp, ndata);
    	return NULL;
    }
    
    /*reading data*/
    for (i=0; i<ndata; i++)
    {
        inp[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }
	printf("ndata=%d\n",ndata);
	    
//        sf_floatread (inp,n1,in);
        for (i=0; i < n1; i++)  {
        for (j=0; j < ntw; j++) {
            if (i+j-m < 0 || i+j-m >= n1) {
            p[j] = 0.;
            } else {
            p[j] = inp[i+j-m];
            }
        }

        if (wind) {
            for (i1=0; i1 < ntw; i1++) {
                p[i1] *= (0.5 - 0.5*cosf(2*SF_PI*i1 / (ntw - 1)));
            }
        }
        
        if (sym) {
            for (i1=0; i1 < ntw; i1++) {
            p[i1] *= wt;
            }
        }
        
        for (i1=ntw; i1 < nt; i1++) {
            p[i1]=0.0;
        }
        
        kiss_fftr (cfg,p,pp);
        
        if (0. != o1) {
            for (i1=0; i1 < nw; i1++) {
            shift = -2.0*SF_PI*i1*dw*o1;
            ce.r = cosf(shift);
            ce.i = sinf(shift);
            pp[i1]=sf_cmul(pp[i1],ce);
            }
        }
        for (i1=0; i1 < nw; i1++) {
            outp[i1*n1+i] = sf_cmplx(pp[i1].r,pp[i1].i);
        }
        }
        sf_complexwrite(outp,n1*nw,out);
        
	}else{
	/*This part is to reconstruct the data given the basis functions and their weights (i.e., TF spectrum)*/
	
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    
    nd2=PyArray_NDIM(arrf1);
    
    npy_intp *sp=PyArray_SHAPE(arrf1);
    	
    for (i=0; i<n1w*2; i++)
    {
        sscc[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }

        sf_complexread(outp,n1*nw,in);

        for (i=0; i < n1; i++)  {
        inp[i] = 0.;
        }
        for (i=0; i < n1; i++)  {
        for (i1=0; i1 < nw; i1++) {
            pp[i1].r = crealf(outp[i1*n1+i]);
            pp[i1].i = cimagf(outp[i1*n1+i]);
        }
        if (0. != o1) {
            for (i1=0; i1 < nw; i1++) {
            shift = +2.0*SF_PI*i1*dw*ow;
            ce.r = cosf(shift);
            ce.i = sinf(shift);
            pp[i1]=sf_cmul(pp[i1],ce);
            }
        }
        
        kiss_fftri(cfg,pp,p);
        
        for (i1=0; i1 < ntw; i1++) {
            p[i1] *= wt;
            if (i+i1-m >= 0 && i+i1-m < n1) {
            inp[i+i1-m] += p[i1];
            }
        }
        }

        for (i=0; i < m; i++) {
        inp[i] = inp[i]/((i+m+1)*1.);
        }
        for (i=m; i < n1-m; i++)  {
        inp[i] = inp[i]/(ntw*1.);
        }
        for (i=n1-m; i < n1; i++) {
        inp[i] = inp[i]/((n1-i+m)*1.);
        }
        sf_floatwrite(inp,n1,out);
	/*sub-function goes here*/
	
    /*Below is the output part*/
    PyArrayObject *vecout;
    npy_intp dims[2];
    
    if(!inv)
    {

	if(ifb) /*if output the basis functions, e.g., Fourier bases in this case*/
	{
	dims[0]=ndata*nw*4+3;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = sscc[i];
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i+ndata*nw*2))) = kbsc[i];
		
	(*((float*)PyArray_GETPTR1(vecout,0+ndata*nw*4))) = w0;
	(*((float*)PyArray_GETPTR1(vecout,1+ndata*nw*4))) = dw;
	(*((float*)PyArray_GETPTR1(vecout,2+ndata*nw*4))) = nw;
	}else{
	dims[0]=ndata*nw*2+3;dims[1]=1;
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = sscc[i];
	printf("w0=%g,dw=%g,nw=%d\n",w0,dw,nw);
	(*((float*)PyArray_GETPTR1(vecout,0+ndata*nw*2))) = w0;
	(*((float*)PyArray_GETPTR1(vecout,1+ndata*nw*2))) = dw;
	(*((float*)PyArray_GETPTR1(vecout,2+ndata*nw*2))) = nw;
	}
	
	
	}else{
	
	dims[0]=n1;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = inp[i];
	}
	
	return PyArray_Return(vecout);
	
}

/*documentation for each functions.*/
static char ntfacfun_document[] = "Document stuff for this C module...";

/*defining our functions like below:
  function_name, function, METH_VARARGS flag, function documents*/
static PyMethodDef functions[] = {
  {"stft1d", stft1d, METH_VARARGS, ntfacfun_document},
  {NULL, NULL, 0, NULL}
};

/*initializing our module informations and settings in this structure
for more informations, check head part of this file. there are some important links out there.*/
static struct PyModuleDef ntfacfunModule = {
  PyModuleDef_HEAD_INIT, /*head informations for Python C API. It is needed to be first member in this struct !!*/
  "ftfacfun",  /*module name (FFT-based time-frequency analysis)*/
  NULL, /*means that the module does not support sub-interpreters, because it has global state.*/
  -1,
  functions  /*our functions list*/
};

/*runs while initializing and calls module creation function.*/
PyMODINIT_FUNC PyInit_ntfacfun(void){
  
    PyObject *module = PyModule_Create(&ntfacfunModule);
    import_array();
    return module;
}
