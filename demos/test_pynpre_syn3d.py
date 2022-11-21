
#  DEMO script (python version) for NPRE (benchmarked with Madagascar version, should be correct)
#  This synthetic data is biased towards the DRR method
#  Please keep an eye on other examples that show superior performance of NPRE 
#
#  Copyright (C) 2022 Yangkang Chen
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details: http://www.gnu.org/licenses/
#  
## generate synthetic data
#This synthetic data was used in Huang et al., 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
import numpy as np
import matplotlib.pyplot as plt
import pydrr as pd #pd: DRR
import pynpre as npre

## generate the synthetic data
a1=np.zeros([300,20])
[n,m]=a1.shape
a3=np.zeros([300,20])
a4=np.zeros([300,20])

k=-1;
a=0.1;
b=1;
pi=np.pi

ts=np.arange(-0.055,0.055+0.002,0.002)
b1=np.zeros([len(ts)])
b2=np.zeros([len(ts)])
b3=np.zeros([len(ts)])
b4=np.zeros([len(ts)])

for t in ts:
    k=k+1;
    b1[k]=(1-2*(pi*30*t)*(pi*30*t))*np.exp(-(pi*30*t)*(pi*30*t));
    b2[k]=(1-2*(pi*40*t)*(pi*40*t))*np.exp(-(pi*40*t)*(pi*40*t));
    b3[k]=(1-2*(pi*40*t)*(pi*40*t))*np.exp(-(pi*40*t)*(pi*40*t));
    b4[k]=(1-2*(pi*30*t)*(pi*30*t))*np.exp(-(pi*30*t)*(pi*30*t));

t1=np.zeros([m],dtype='int')
t3=np.zeros([m],dtype='int')
t4=np.zeros([m],dtype='int')
for i in range(m):
  t1[i]=np.round(140);
  t3[i]=np.round(-6*i+180);
  t4[i]=np.round(6*i+10);
  a1[t1[i]:t1[i]+k+1,i]=b1; 
  a3[t3[i]:t3[i]+k+1,i]=b1; 
  a4[t4[i]:t4[i]+k+1,i]=b1; 

temp=a1[0:300,:]+a3[0:300,:]+a4[0:300,:];

shot=np.zeros([300,20,20])
for j in range(20):
    a4=np.zeros([300,20]);
    for i in range(m):
    	t4[i]=np.round(6*i+10+3*j); 
    	a4[t4[i]:t4[i]+k+1,i]=b1;
  
    	t1[i]=np.round(140-2*j);
    	a1[t1[i]:t1[i]+k+1,i]=b1;

    shot[:,:,j]=a1[0:300,:]+a3[0:300,:]+a4[0:300,:];

d0=shot

## add noise
[n1,n2,n3]=d0.shape
np.random.seed(201415)
n=0.2*np.random.randn(n1,n2,n3);
dn=d0+n;
print(np.std(dn))


# d1=pd.drr3d(dn,0,120,0.004,3,100);	#RR
# d1=dn;
# d1=npre.npre3d(dn,n1win=50,n2win=20,n3win=20,rect=[3,3,3],niter=50,mode=1)
d1=npre.npre3d(dn,rect=[2,2,2],niter=50,mode=0) #mode=2 is better
noi1=dn-d1;

d2=pd.drr3d(dn,0,120,0.004,3,3);	#DRR
noi2=dn-d2;

## compare SNR
print('SNR of RR is %g'%pd.snr(d0,d1,2));
print('SNR of DRR is %g'%pd.snr(d0,d2,2));

## plotting
fig = plt.figure(figsize=(8, 7))
ax=fig.add_subplot(3, 2, 1)
plt.imshow(dn.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noisy data');
ax=fig.add_subplot(3, 2, 3)
plt.imshow(d1.reshape(n1,n2*n3,order='F'),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Denoised (NPRE, SNR=%.4g dB)'%pd.snr(d0,d1,2));
ax=fig.add_subplot(3, 2, 4)
plt.imshow(noi1.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noise (NPRE)');
ax=fig.add_subplot(3, 2, 5)
plt.imshow(d2.reshape(n1,n2*n3,order='F'),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Denoised (DRR, SNR=%.4g dB)'%pd.snr(d0,d2,2));
ax=fig.add_subplot(3, 2, 6)
plt.imshow(noi2.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.1, 0.1),aspect='auto');ax.set_xticks([]);ax.set_yticks([]);
plt.title('Noise (DRR)');
plt.savefig('test_pynpre_syn3d.png',format='png',dpi=300);
plt.show()

## Save as binary files
dn=np.float32(dn)
fid = open ("syn3d_dn.bin", "wb") #binary file format, int
fid.write(dn.flatten(order='F'))

d0=np.float32(d0)
fid = open ("syn3d_dc.bin", "wb") #binary file format, int
fid.write(d0.flatten(order='F'))

d1=np.float32(d1)
fid = open ("syn3d_d1.bin", "wb") #binary file format, int
fid.write(d1.flatten(order='F'))

#sfbin2rsf <syn3d_dn.bin n1=300 n2=400 d1=1 d2=1 o1=0 o2=0 |sfgrey |sfpen

