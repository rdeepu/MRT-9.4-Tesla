                #DATA CLUSTERING
    
import nifti as nt

import numpy as np

import scipy as sp

from scipy.cluster.vq import kmeans,vq

from sklearn.cluster import KMeans

np.set_printoptions(threshold = 'nan')

rd = nt.NiftiImage('/opt/ime-157/data/Lekshmi-training/pca_crtx_trial.nii')

rd = np.array(rd.data,dtype=float)

rd[rd==0] = 'nan'

sz = rd.shape

ncols = sz[2]

nrows = sz[3]

print nrows

print ncols

rd = np.reshape(rd,(nrows*ncols),1)

print rd.shape

kmeans = KMeans(n_clusters=7, random_state=0).fit(rd)

kmeans.labels

#codebook, distortion = kmeans(rd,7)

#labels,_ = vq(rd,codebook)

#labels

#rd = np.reshape(labels,(ncols,nrows))

#img_wrt = nt.NiftiImage(rd)
            
#img_wrt.save('/opt/ime-157/data/Lekshmi-training/pca_cluster.nii')
            

