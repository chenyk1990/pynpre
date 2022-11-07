from npre3dcfun import *

def npre3d(din):
	"""
	Non-stationary predictive filtering (3D)
	
	By Yangkang Chen
	Nov 6, 2022
	"""
	import numpy as np
	
	din=np.float32(din)
	ndata=din.size;
	dout=fxynpre(din,ndata,3.0)
	
	return dout

