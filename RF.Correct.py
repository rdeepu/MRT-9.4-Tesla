 #B1+ CORRECTION FOR PARAMETRIC MAPS
                          

import nifti as nt

import numpy as np

np.set_printoptions(threshold = 'nan')

rdt1 = nt.NiftiImage('/opt/data/MTransfer/rev-agarMT/results/mt/t1map.nii')

rdt1 = np.array(rdt1.data,dtype=float)

rdfm = nt.NiftiImage('/opt/data/MTransfer/rev-agarMT/results/flmap/flipmap.nii')

rdfm = np.array(rdfm.data,dtype=float)

rdfm[rdfm == 0] = 0.01

rdcor = rdt1 / rdfm

img_wrt = nt.NiftiImage(rdcor)
            
img_wrt.save('/opt/data/MTransfer/rev-agarMT/results/mt/fmcor_t1map.nii')
