  #T1 MAP and M0 MAP FROM MULTIANGLE SPGR DATA
                        

import nifti as nt

import numpy as np

import scipy.stats as st

np.set_printoptions(threshold = 'nan')

tr = 50

rd = nt.NiftiImage('/opt/data/MTransfer/rev-agarMT/FA-Nifti/mt/4Dmt.nii')

rd = np.array(rd.data,dtype=float)

rd[rd < 30] = 'nan'

sz = rd.shape

print sz

x = np.zeros((1, sz[0]),dtype = float)

y = np.zeros((1, sz[0]),dtype = float) 

t1map = np.zeros((sz[1],sz[2],sz[3]),dtype = float)

mmap  = np.zeros((sz[1],sz[2],sz[3]),dtype = float)

for sl in range(0, sz[1]):
    
    for cl in range(0, sz[3]):
                
        for rw in range(0, sz[2]):
            
            fa = 10
            
            for an in range(0, 8):
                                
                ra = np.deg2rad(fa)
                                                    
                x[0, an] = rd[an,sl,rw,cl] / np.tan(ra)
                                                                 
                y[0, an] = rd[an,sl,rw,cl] / np.sin(ra)
                               
                fa = fa + 10
                    
            slope, intercept, r_value, p_value, std_err = st.linregress(x,y)
                    
            t1map[sl,rw,cl] =  -tr / np.log(slope)
            
            mmap[sl,rw,cl] = intercept / (1 - np.exp(-tr / t1map[sl,rw,cl]))
            
            
mmap[np.isnan(mmap)] = 0

mmap[mmap < 0] = 0

t1map[np.isnan(t1map)] = 0

t1map[t1map < 0] = 0
            
img_wrt = nt.NiftiImage(t1map)
            
img_wrt.save('/opt/data/MTransfer/rev-agarMT/results/mt/t1map.nii')
            
img_wrt = nt.NiftiImage(mmap)
            
img_wrt.save('/opt/data/MTransfer/rev-agarMT/results/mt/mmap.nii')
            
            
