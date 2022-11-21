from npre3dcfun import *
import numpy as np

def npre3d(din,d1=0.004,o1=0,n1win=None,n2win=None,n3win=None,nsx=2,nsy=2,niter=100,mode=0,rect=[5,5,5],opt=1,sym=0,fs=0,verb=1,Nfrac=6.0,Ntimes=5.0,pow=0.5,r1=0.5,r2=0.5,r3=0.5):
	"""
	Non-stationary predictive filtering (3D)
	
	By Yangkang Chen
	Nov 6, 2022
	
	INPUT
	din:	input data
	n1win:  first window length
	n2win:  second window length
	n3win:  third window length
	nsx:	number of shifts in non-causal prediction filtering 
	nsy: 	number of shifts in non-causal prediction filtering 
	niter:  number of iterations 
	mode: 	predictive filtering mode; default: non-stationary
			mode=0: 	Non-stationary Predictive Filtering
			mode=0, fs=1: Non-stationary Predictive Filtering with frequency-dependent smoothing
			mode=1:		Regularized Non-stationary Autoregression
			mode=2: 	Stationary Predictive Filtering	
	
	OUTPUT
	dout:	denoised data
	
	"""
	if din.ndim==2:	#for 2D problems
		din=np.expand_dims(din, axis=2)
		
	[n1,n2,n3]=din.shape
	
	din=np.float32(din.flatten(order='F'))
	
	rect1=rect[0]
	rect2=rect[1]
	rect3=rect[2]
	
	if n1win is None:
		n1win=n1
	if n2win is None:
		n2win=n2
	if n3win is None:
		n3win=n3
	
	dout=Cfxynpre(din,n1,n2,n3,n1win,n2win,n3win,nsx,nsy,niter,mode,rect1,rect2,rect3,opt,sym,fs,verb,d1,o1,Nfrac,Ntimes,pow,r1,r2,r3);
	
	dout=dout.reshape(n1,n2,n3,order='F')
	
	if n3==1:	#for 2D problems
		dout=np.squeeze(dout)

	return dout

