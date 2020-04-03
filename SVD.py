                                       #SINGULAR VALUE DECOMPOSITION
import nifti as nt

import numpy as np

import scipy as sp

np.set_printoptions(threshold = 'nan')

rd = nt.NiftiImage('/opt/data/cort_map/Animal/Odd/2016-05-18/IR-TSE/SL-19/results/python_test/edtd.nii')

rd = np.array(rd.data,dtype = int)

sz = rd.shape

print sz

nvols = sz[0]

nslices = sz[1]

ncols = sz[2]

nrows = sz[3]

rd = np.reshape(rd,(nrows*ncols*nslices,nvols))

U,S,V = np.linalg.svd(rd,full_matrices=0)

print U.shape

print S.shape

print V.shape

rd = np.reshape(U,(nvols,nslices,ncols,nrows))

print rd.shape

wr = nt.NiftiImage(rd)

wr.save('/opt/data/cort_map/Animal/Odd/2016-05-18/IR-TSE/SL-19/results/python_test/srvd.nii')

