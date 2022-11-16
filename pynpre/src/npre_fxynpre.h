#ifndef _fxynpre_h
#define _fxynpre_h


// #include <rsf.h>
// #include "fft1.h"
// #include "cdivn.h"
// #include "cdivn_fs.h"
// #include "cdivn_rnar.h"
// #include "cdivnn.h"
// #include "cdivns.h" /*stationary division in traditional FXY*/
void fxynpre0(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		bool sym, 
		bool opt,
		bool verb);
/*<Non-stationary predictive filtering in 3D (without frequency-dependent smoothing)>*/


void fxynpre(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		float Nfrac, /*starting varying from 1/Nfrac*F_nyquist*/
		float Ntimes, /*largest radius is Ntimes of the ref radius*/
		float pow,   /*sharp rate of the varying curve, the lower the sharper*/
		bool sym, 
		bool opt,
		bool verb);
/*<Non-stationary predictive filtering in 3D with frequency-dependent smoothing>*/


void fxynpre1(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		bool sym, 
		bool opt,
		bool verb);
/*<Non-stationary predictive filtering in 3D>*/


void fxynpre3(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		float **rct /* triangle lengths [ndim][nd] */,
        int **sft 	/* triangle shifts [ndim][nd] */,
		bool sym, 
		bool opt,
		bool verb);
/*<Non-stationary predictive filtering in 3D with arbitrary smoothing strategies>*/


void fxypre(float **dtime /*input and output data*/, 
		int n2,
		int nx,
		int ny,
		int n1,
		float d1,
		float o1, 		
		int niter,
		int nsx,
		int nsy,
		int *rect,
		bool sym, 
		bool opt,
		bool verb);
/*<Stationary predictive filtering in 3D>*/

#endif
