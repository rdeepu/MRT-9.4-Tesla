                            #Transmit RF Inhomogeneity Map Calculation employing AFI
                                              
import nifti as nt

import numpy as np

np.set_printoptions(threshold = 'nan')

nu = 5

fa = 60

s1 = nt.NiftiImage('/opt/data/Scanner-Calibration/16-02-12/AFI/3_db_afi_seql/niifti/01.nii')

s1 = np.array(s1.data, dtype = float)

s1[s1 < 40] = 'nan'

s2 = nt.NiftiImage('/opt/data/Scanner-Calibration/16-02-12/AFI/3_db_afi_seql/niifti/02.nii')

s2 = np.array(s2.data, dtype = float)

s2[s2 < 40] = 'nan'

r = s2/s1

r = ((r*nu)-1)/(nu-r)

r = np.rad2deg(np.arccos(r))
            
r = r/fa

wr = nt.NiftiImage(r)

wr.save('/opt/data/Scanner-Calibration/16-02-12/AFI/3_db_afi_seql/niifti/flipmap-I.nii')
