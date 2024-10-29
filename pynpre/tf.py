from ftfacfun import *

def st1d(din,dt=0.004,inv=0,opt=1,sym=0,ntw=7,ot=0,wind=0,verb=1):
    """
    st1d: S transform for 1D signals
    
    INPUT
    din: input trace (1D, n1/nt)
    inv: flag of inverse transform (0: forward; 1: inverse; default: 0)
    opt:
    sym:
    
    
    OUTPUT
    dout: output Time-frequency spectrum (first axis: Time; second axis: Frequency; third axis: imag/real)
    NOTE: 1) dout is of size n1*nw*2 (e.g., dout.reshape([n1,nw,2],order='F'))
          2) The first component in 3rd axis is Imaginary; and the second component is Real (Remember it when constructing a complex number)
    
    EXAMPLE
    demos/test_pyntfa_syn1d_st.py
    
    HISTORY
    Original version by Yangkang Chen, Oct 28, 2024
    
    """
    import numpy as np
    
    din=np.float32(din)
    n1=din.shape[0];
    
    print(n1)
    
    dtmp=stft1d(din.flatten(order='F'),n1,verb,wind,inv,sym,opt,ntw,dt,ot)
    
    
    return dtmp
    
    

def stft1d(din,dt=0.004,niter=100,rect=10,ifb=0,ifn=0,inv=0,verb=1,alpha=0,rect1=None,rect2=None):
    """
    stft1d: Short-Time Fourier transform for 1D signals
    
    INPUT
    din: input trace (1D, n1/nt)
    fhi: High frequency in band, default is Nyquist
    flo: Low frequency in band, default is 0
    inv: flag of inverse transform (0: forward; 1: inverse; default: 0)
    
    OUTPUT
    dout: output Time-frequency spectrum (first axis: Time; second axis: Frequency; third axis: imag/real)
    NOTE: 1) dout is of size n1*nw*2 (e.g., dout.reshape([n1,nw,2],order='F'))
          2) The first component in 3rd axis is Imaginary; and the second component is Real (Remember it when constructing a complex number)
    
    EXAMPLE
    demos/test_pyntfa_syn1d_st.py
    
    HISTORY
    Original version by Yangkang Chen, Oct 28, 2024
    
    """
    import numpy as np
    
    
    
    
    
    
    
    return
    


