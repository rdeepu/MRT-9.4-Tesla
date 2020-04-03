import nifti as nt

import numpy as np

import scipy as sp

from scipy.optimize import curve_fit

np.set_printoptions(threshold = 'nan')

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

rd = nt.NiftiImage('/opt/ime-157/data/Stroke-Project/rat/T2-SE/niifti/edited.nii')

rd = np.array(rd.data,dtype = float)

sz = rd.shape

print sz

x = np.array([10.6, 21.2, 31.8, 42.4, 53,63.6,74.2,84.8,95.4])

ttwo = np.zeros((1,sz[2],sz[3]),dtype = float) 

mzero = np.zeros((1,sz[2],sz[3]),dtype = float) 

for sl in range(4,5):

    for cl in range(0,sz[2]):
    
        for rw in range(0,sz[3]):
           
            y = rd[:,sl,cl,rw]   
                                   
            def func(x,p1,p2):
            
                return p1*(np.exp(-x/p2))
            
            popt,pcov = curve_fit(func,x, y, p0 = (3000,50)) 
            
            ttwo[:,cl,rw] = popt[1]
            
            #print popt[1]
            
            #print popt[0]
            
            mzero[:,cl,rw] = popt[0]

img_wrt = nt.NiftiImage(ttwo)
            
img_wrt.save('/opt/ime-157/data/Stroke-Project/rat/T2-SE/results/t2map_rat.nii')
            
img_wrt = nt.NiftiImage(mzero)
            
img_wrt.save('/opt/ime-157/data/Stroke-Project/rat/T2-SE/results/m0map_rat.nii')

