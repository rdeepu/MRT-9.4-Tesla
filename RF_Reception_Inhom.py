import nifti as nt

import numpy as np

np.set_printoptions(threshold = 'nan')

s1 = nt.NiftiImage('/opt/data/Scanner-Calibration/16-02-12/b.coil/01.nii')

s1 = np.array(s1.data, dtype = float)

s1[s1 < 40] = 'nan'

s2 = nt.NiftiImage('/opt/data/Scanner-Calibration/16-02-12/r.coil/01.nii')

s2 = np.array(s2.data, dtype = float)

s2[s2 < 40] = 'nan'

c = s1/s2

wr = nt.NiftiImage(c)

wr.save('/opt/data/Scanner-Calibration/16-02-12/AFI/3_db_afi_seql/niifti/flipmap-II.nii')

