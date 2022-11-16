#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <numpy/arrayobject.h>

#include "npre_win.h"
#include "npre_memcpy.h"
#include "npre_fxynpre.h"
#include "npre_alloc.h"

int kiss_fft_next_fast_size(int n);

/*from user/chenyk/vecoper.c */
static PyObject *Cfxynpre(PyObject *self, PyObject *args){
    
	/**initialize data input**/
    int nd, nd2;
    
    PyObject *f1=NULL;
    PyObject *arrf1=NULL;

	int ndata;	/*integer parameter*/
	float fpar; /*float parameter*/
    int ndim, i;
    float *data;

	/*sub-function goes here*/
    int  verb=0, sym=0, opt=0, fs=0; /*fs: if apply frequency-dependent smoothing*/
    int   ix,iy,nx,ny;
    int rect[9], rect1, rect2, rect3;    
    int niter,mode;    
    int nsx,nsy,j;
    int n1win,n2win,n3win,n1pad,n2pad,n3pad;/*window sizes,and padding sizes*/
    int nw1,nw2,nw3,iw1,iw2,iw3;/*number of windows in first and second axis,indices of nw1,nw2*/
    int s1,s2,s3; /*indices in selecting windows*/
    int nov1,nov2,nov3,ov1,ov2,ov3;/*overlapping or non-overlapping points*/
    float o1,o11, d1,r1,r2,r3;
    float ***din, ***dout,***dtmp; 
	int nt,n1;
	
    float Nfrac; /*starting varying from nw/Nfrac*F_nyquist*/
    float Ntimes; /*largest radius is Ntimes of the ref radius*/
    float pow;   /*sharp rate of the varying curve, the lower the sharper*/

	printf("verb=%d,sym=%d,opt=%d,fs=%d\n",verb,sym,opt,fs);
	PyArg_ParseTuple(args, "Oiiiiiiiiiiiiiiiiiffffffff",&f1,&n1,&nx,&ny,&n1win,&n2win,&n3win,&nsx,&nsy,&niter,&mode,&rect1,&rect2,&rect3,&opt,&sym,&fs,&verb,&d1,&o1,&Nfrac,&Ntimes,&pow,&r1,&r2,&r3);
	printf("verb=%d,sym=%d,opt=%d,fs=%d\n",verb,sym,opt,fs);
	printf("n1=%d,nx=%d,ny=%d\n",n1,nx,ny);
	printf("n1win=%d,n2win=%d,n3win=%d\n",n1win,n2win,n3win);
	printf("nsx=%d,nsy=%d\n",nsx,nsy);
	printf("rect=%d,%d,%d\n",rect1,rect2,rect3);
	printf("r=%g,%g,%g\n",r1,r2,r3);


    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    
    nd2=PyArray_NDIM(arrf1);
    
    npy_intp *sp=PyArray_SHAPE(arrf1);

	ndata=n1*nx*ny;
	
	data  = (float*)malloc(ndata * sizeof(float));
	
    if (*sp != ndata)
    {
    	printf("Dimension mismatch, N_input = %d, N_data = %d\n", *sp, ndata);
    	return NULL;
    }
    
    /*reading data*/
    for (i=0; i<ndata; i++)
    {
        data[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }


	printf("verb=%d,sym=%d,opt=%d,fs=%d\n",verb,sym,opt,fs);
	
	printf("Nfrac=%g,Ntimes=%g,pow=%g\n",Nfrac,Ntimes,pow);
	 
//     verb=false; /* Verbosity flag */ 
//     sym=false; /* y, symmetric scaling for Hermitian FFT */
//     opt=true;  /* y, determine optimal size for efficiency */   
//     fs=false;  /* y, determine frequency-dependent smoothing */

    
// 	if(fs)
// 	{
// //     if (!np_getfloat("Nfrac",&Nfrac))    
//     Nfrac=6.0;     /*frequency-dependent smoothing starts from 1/Nfrac * (Nyquist frequency) */
// //     if (!np_getfloat("Ntimes",&Ntimes))    
//     Ntimes=5.0;  /*Maximum smoothing radius is Ntimes*(reference smoothing radius) */
// //     if (!np_getfloat("pow",&pow))    
//     pow=0.5;  /*fraction parameter*/
// 	}

//     for (j=0; j < 3; j++) {
// // 	snprintf(key,6,"rect%d",j+1);
// // 	if (!np_getint(key,rect+j)) 
// 	rect[j]=1; 
// 	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
//     }
    rect[0]=rect1;
    rect[1]=rect2;
    rect[2]=rect3;
    
//  	if (!np_histint  (Fin,"n1",&n1)) 
//  	n1=100;
// // 	if (!np_histfloat(Fin,"d1",&d1)) 
// 	d1=1.;
// // 	if (!np_histfloat(Fin,"o1",&o1)) 
// 	o1=0.;
// //  	if (!np_histint  (Fin,"n2",&nx)) 
//  	nx=100;
// //  	if (!np_histint  (Fin,"n3",&ny)) 
//  	ny=100;
	printf("n1=%d,nx=%d,ny=%d\n",n1,nx,ny);

//     if (!np_getint("niter",&niter)) 
//     niter=100;
    /* number of iterations */

	/*windowing parameters*/
//     if (!np_getint("n1win",&n1win)) 
//     n1win=n1; /*first window length*/
// //     if (!np_getint("n2win",&n2win)) 
//     n2win=nx; /*second window length*/	
// //     if (!np_getint("n3win",&n3win)) 
//     n3win=ny; /*second window length*/	
//     
// // 	if (!np_getfloat("r1",&r1)) 
// 	r1=0.5;		  /*first overlapping ratio*/
// // 	if (!np_getfloat("r2",&r2)) 
// 	r2=0.5;		  /*second overlapping ratio*/
// // 	if (!np_getfloat("r3",&r3)) 
// 	r3=0.5;		  /*third overlapping ratio*/
// 
// //     if (!np_getint("mode",&mode)) 
//     mode=0; /*predictive filtering mode; default: non-stationary*/	
//     
// // 	if(! np_getint("nsx",&nsx)) 
// 	nsx=2; /* number of shifts in non-causal prediction filtering */
// // 	if(! np_getint("nsy",&nsy)) 
// 	nsy=2; /* number of shifts in non-causal prediction filtering */

	nov1=(1-r1)*n1win;nov2=(1-r2)*n2win;nov3=(1-r3)*n3win;/*non-overlapping size 1,2,3*/
	ov1=r1*n1win;ov2=r2*n2win;ov3=r3*n3win;		/*overlapping size 1,2,3*/
	
	printf("nsx=%d,nsy=%d\n",nsx,nsy);

	/*padding */	
	n1pad=n1win;nw1=1;while(n1pad<n1){n1pad=n1pad+nov1;nw1++;	 }   
	n2pad=n2win;nw2=1;while(n2pad<nx){n2pad=n2pad+nov2;nw2++;	 }   	
	n3pad=n3win;nw3=1;while(n3pad<ny){n3pad=n3pad+nov3;nw3++;	 }   	
			    
    din =                  np_floatalloc3  (n1pad,n2pad,n3pad); 
    dout =                  np_floatalloc3  (n1pad,n2pad,n3pad);           
    dtmp =                   np_floatalloc3  (n1win,n2win,n3win);    

//  	for(iy=0; iy<ny; iy++)			    
//  	for(ix=0; ix<nx; ix++)
// 		np_floatread (din[iy][ix],n1,Fin);

	mcp3d1d(din,data,0,0,0,0,0,0,n1,nx,ny,n1,nx,ny);
	
	for(iw3=0;iw3<nw3;iw3++)
	for(iw2=0;iw2<nw2;iw2++)
		for(iw1=0;iw1<nw1;iw1++)
			{
			s1=iw1*nov1;s2=iw2*nov2;s3=iw3*nov3;
/*				printf("Hello,n1win=%d,n2win=%d",n1win,n2win);
			printf("nw2=%d,nw1=%d,d1=%g,o11=%g,niter=%d,ns=%d",nw2,nw1,d1,o11,niter,ns);*/
			mcp3d(dtmp,din,0,0,0,s1,s2,s3,n1win,n2win,n3win);
			o11=o1+s1*d1;
			if(mode==0)
    		{printf("Non-stationary Predictive Filtering\n");
    		if(fs)
    		{fxynpre(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,Nfrac,Ntimes,pow,sym,opt,verb);}
    		else
    		{fxynpre0(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,sym,opt,verb);}}
    		else{
    		if(mode==1)
    		{printf("Regularized Non-stationary Autoregression\n");
    		fxynpre1(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,sym,opt,verb);}
    		else
    		{if(mode==2){printf("Stationary Predictive Filtering\n");
    		fxypre(dtmp[0],n2win*n3win,n2win,n3win,n1win,d1,o11,niter,nsx,nsy,rect,sym,opt,verb);}else{printf("wrong mode\n\n\n");}}
    		}    		
    		win_weight3d(dtmp,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3);
    		mcp_ad3d(dout,dtmp,s1,s2,s3,0,0,0,n1win,n2win,n3win);
						
			}
	
	mcp3d3d1d(data,dout,0,0,0,0,0,0,n1,nx,ny,n1,nx,ny);

	
    /*Below is the output part*/
    PyArrayObject *vecout;
	npy_intp dims[2];
	dims[0]=ndata;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = data[i];
	
	return PyArray_Return(vecout);
	
}


/*documentation for each functions.*/
static char npre3dcfun_document[] = "Document stuff for this C module...";

/*defining our functions like below:
  function_name, function, METH_VARARGS flag, function documents*/
static PyMethodDef functions[] = {
  {"Cfxynpre", Cfxynpre, METH_VARARGS, npre3dcfun_document},
  {NULL, NULL, 0, NULL}
};

/*initializing our module informations and settings in this structure
for more informations, check head part of this file. there are some important links out there.*/
static struct PyModuleDef npre3dcfunModule = {
  PyModuleDef_HEAD_INIT, /*head informations for Python C API. It is needed to be first member in this struct !!*/
  "npre3dcfun",  /*module name*/
  NULL, /*means that the module does not support sub-interpreters, because it has global state.*/
  -1,
  functions  /*our functions list*/
};

/*runs while initializing and calls module creation function.*/
PyMODINIT_FUNC PyInit_npre3dcfun(void){
  
    PyObject *module = PyModule_Create(&npre3dcfunModule);
    import_array();
    return module;
}
