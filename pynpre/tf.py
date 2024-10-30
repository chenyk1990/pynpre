from ftfacfun import *

def st1d(din,dt=0.004,inv=0,flo=0,fhi=0.5,verb=1):
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
    
    dtmp=st1dc(din.flatten(order='F'),n1,verb,inv,flo,fhi,dt)
    
    if inv==0:
        nw=np.int32((dtmp.size-3)/n1/2);
        dout=dtmp[0:n1*nw*2]
        w0=dtmp[n1*nw*2]    #different from PYntfa
        dw=dtmp[n1*nw*2+1]  #different from PYntfa
        nw2=np.int32(dtmp[n1*nw*2+2])
        print(dtmp.size,w0,dw,nw2)
        if nw2 != nw:
            print('nw=',nw,'nw2=',nw2,'dimension discrepancy')
        else:
            print('dimension consistent')
        return dout,w0,dw,nw
    else:
        return dtmp
    
    

def stft1d(din,dt=0.004,inv=0,opt=1,sym=0,ntw=7,ot=0,wind=0,verb=1):
    """
    stft1d: Short-Time Fourier transform for 1D signals
    
    INPUT
    din: input trace (1D, n1/nt)

    inv: flag of inverse transform (0: forward; 1: inverse; default: 0)
    
    OUTPUT
    dout: output Time-frequency spectrum (first axis: Time; second axis: Frequency; third axis: imag/real)
    NOTE: 1) dout is of size n1*nw*2 (e.g., dout.reshape([n1,nw,2],order='F'))
          2) The first component in 3rd axis is Real; and the second component is Imaginary (Remember it when constructing a complex number)
          3) The third axis real/imag sequence of stft1d is opposite to the case in ntfa1d (pyntfa)
    
    EXAMPLE
    demos/test_pyntfa_syn1d_st.py
    
    HISTORY
    Original version by Yangkang Chen, Oct 28, 2024
    
    """
    import numpy as np

    n1=din.shape[0];
    din=np.float32(din)

    
    print(n1)
    
    dtmp=stft1dc(din.flatten(order='F'),n1,verb,wind,inv,sym,opt,ntw,dt,ot)
    
    
    if inv==0:
        nw=np.int32((dtmp.size-3)/n1/2);
        dout=dtmp[0:n1*nw*2]
        w0=dtmp[n1*nw*2]    #different from PYntfa
        dw=dtmp[n1*nw*2+1]  #different from PYntfa
        nw2=np.int32(dtmp[n1*nw*2+2])
        print(dtmp.size,w0,dw,nw2)
        if nw2 != nw:
            print('nw=',nw,'nw2=',nw2,'dimension discrepancy')
        else:
            print('dimension consistent')
        return dout,w0,dw,nw
    else:
        return dtmp
    


