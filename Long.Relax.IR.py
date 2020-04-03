    #T1 AND M0 MAP FROM INVERSION-RECOVERY DATA.
                               
import nifti as nt

import numpy as np

import scipy as sp

from scipy.optimize import curve_fit

np.set_printoptions(threshold = 'nan')

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

rd = nt.NiftiImage('/opt/Stroke Project/rat/IR-HighRes_fatSat/niifti/edited.nii')

rd = np.array(rd.data,dtype = float)

sz = rd.shape

print sz

x = np.array([500, 750, 1000, 1500, 2000, 3000, 5000, 7000])

tone = np.zeros((1,sz[2],sz[3]),dtype = float) 

mzero = np.zeros((1,sz[2],sz[3]),dtype = float) 

for sl in range(0,sz[1]):

    for cl in range(0,sz[2]):
    
        for rw in range(0,sz[3]):
           
            y = rd[:,sl,cl,rw]   
            
            mi =  min(y)
            
            ix = np.where(y==mi)
            
            ix = np.array(ix)
            
            ix =ix + 1
                                
            print rd[:,sl,cl,rw]
            
            print sl,cl,rw
            
            y[0:ix] = -y[0:ix]
            
            print y
            
            def func(x,p1,p2):
            
                return p1*(1-2*(np.exp(-x/p2)))
            
            popt,pcov = curve_fit(func,x, y, p0 = (1500,500)) 
            
            tone[:,cl,rw] = popt[1]
            
            print popt[1]
            
            print popt[0]
            
            mzero[:,cl,rw] = popt[0]

img_wrt = nt.NiftiImage(tone)
            
img_wrt.save('/opt/data/pdata/t1map_rat.nii')
            
img_wrt = nt.NiftiImage(mzero)
            
img_wrt.save('/opt/data/pdata/m0map_rat.nii')


