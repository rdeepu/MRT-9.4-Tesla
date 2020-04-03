   #PRINCIPAL COMPONENT ANALYSIS

from matplotlib.mlab import PCA

import nifti as nt

import numpy as np

import scipy as sp

np.set_printoptions(threshold = 'nan')

rd = nt.NiftiImage('/opt/ime-157/data/Lekshmi-training/edtd.nii')

rd = np.array(rd.data)

sz = rd.shape

print sz

nvols = sz[0]

nslices = sz[1]

ncols = sz[2]

nrows = sz[3]

rd = np.reshape(rd,(nrows*ncols*nslices,nvols))

res =  PCA(rd)

rd = np.reshape(res.a,(nvols,nslices,ncols,nrows))

print rd.shape

wr = nt.NiftiImage(rd)

wr.save('/opt/ime-157/data/Lekshmi-training/pca_trial.nii')
